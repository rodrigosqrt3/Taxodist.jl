"""
    TaxodistTree

A rooted classification tree in the layout of an ape `phylo` object:
`edge` is an `(n_edges, 2)` matrix of `(parent, child)` node numbers with tips
numbered `1:n_tips` and internal nodes from `n_tips + 1` (the root) in
preorder; `edge_length` holds branch lengths; `tip_label` and `node_label`
hold input labels and clade names. Use `newick(tree)` to export it and
`cophenetic(tree)` for its patristic distances.
"""
struct TaxodistTree
    edge::Matrix{Int}
    edge_length::Vector{Float64}
    tip_label::Vector{String}
    node_label::Vector{String}
end

n_tips(tree::TaxodistTree) = length(tree.tip_label)
n_nodes(tree::TaxodistTree) = length(tree.node_label)

function Base.show(io::IO, ::MIME"text/plain", tree::TaxodistTree)
    print(
        io,
        "TaxodistTree with $(n_tips(tree)) tips and $(n_nodes(tree)) internal nodes " *
        "(root: $(repr(first(tree.node_label))))",
    )
end

function _tree_children(tree::TaxodistTree)
    children = Dict{Int,Vector{Tuple{Int,Float64}}}()
    for row in axes(tree.edge, 1)
        push!(
            get!(children, tree.edge[row, 1], Tuple{Int,Float64}[]),
            (tree.edge[row, 2], tree.edge_length[row]),
        )
    end
    return children
end

function _newick_label(label::AbstractString)
    isempty(label) && return ""
    occursin(r"[\s(),:;'\[\]]", label) || return String(label)
    return "'" * replace(label, "'" => "''") * "'"
end

"""
    newick(tree::TaxodistTree)

Return the tree in Newick format with node labels and branch lengths.
"""
function newick(tree::TaxodistTree)
    children = _tree_children(tree)
    tips = n_tips(tree)
    function render(node)
        node <= tips && return _newick_label(tree.tip_label[node])
        inner = join(
            ["$(render(child)):$(branch)" for (child, branch) in children[node]],
            ",",
        )
        return "($(inner))$(_newick_label(tree.node_label[node - tips]))"
    end
    return render(tips + 1) * ";"
end

"""
    cophenetic(tree::TaxodistTree)

Return the patristic distances between the tips of a tree as a
`TaxonomicDistanceMatrix`.
"""
function cophenetic(tree::TaxodistTree)
    children = _tree_children(tree)
    tips = n_tips(tree)
    depth = Dict(tips + 1 => 0.0)
    parent = Dict{Int,Int}()
    queue = [tips + 1]
    while !isempty(queue)
        node = popfirst!(queue)
        for (child, branch) in get(children, node, Tuple{Int,Float64}[])
            depth[child] = depth[node] + branch
            parent[child] = node
            push!(queue, child)
        end
    end
    function ancestors(node)
        path = [node]
        while haskey(parent, node)
            node = parent[node]
            push!(path, node)
        end
        return path
    end
    paths = [ancestors(tip) for tip in 1:tips]
    values = zeros(tips, tips)
    for i in 1:tips, j in (i + 1):tips
        seen = Set(paths[i])
        common = first(node for node in paths[j] if node in seen)
        values[i, j] = values[j, i] = depth[i] + depth[j] - 2 * depth[common]
    end
    return TaxonomicDistanceMatrix(values, tree.tip_label)
end

"""
    taxo_tree(taxa; verbose=false, progress=true)

Build a tree whose internal nodes are the named clades of the source
classification, directly from the lineages of the supplied taxa (a vector of
names, a `TaxodistResolution`, or a `TaxodistBundle`). Clades with a single
descendant among the supplied taxa are collapsed.

Branch lengths make the tree ultrametric with patristic distances equal to
`distance_matrix`: a clade at depth `h` sits at height `1 / (2h)` and tips at
height zero; taxa with identical lineages are joined by zero-length branches.
Taxa without a lineage are dropped with a warning; all lineages must share the
same root node.
"""
function taxo_tree(taxa; verbose::Bool=false, progress::Bool=true)
    labels, lineages = _matrix_input(taxa; verbose=verbose, progress=progress)
    missing_rows = findall(isnothing, lineages)
    if !isempty(missing_rows)
        @warn "$(length(missing_rows)) taxa without a lineage were dropped from the tree: " *
              join(labels[missing_rows], ", ")
        keep = findall(!isnothing, lineages)
        labels = labels[keep]
        lineages = lineages[keep]
    end
    length(lineages) >= 2 || throw(ArgumentError(
        "At least two taxa with lineages are required to build a tree.",
    ))
    if any(isempty, lineages) || length(unique(first.(lineages))) > 1
        throw(ArgumentError(
            "All lineages must share the same root node; taxa without a common " *
            "ancestor have infinite distance and cannot be placed in one tree.",
        ))
    end

    # A node is (tip=index) or (label, height, children).
    function build(indices, depth)
        clade = lineages[first(indices)][depth]
        ended = [i for i in indices if length(lineages[i]) == depth]
        continuing = [i for i in indices if length(lineages[i]) > depth]

        children = Any[]
        if length(ended) == 1
            push!(children, (tip=only(ended),))
        elseif length(ended) > 1
            # Identical lineages: join at height zero so their distance is 0.
            push!(children, (
                label=isempty(continuing) ? clade : "",
                height=0.0,
                children=Any[(tip=i,) for i in ended],
            ))
        end
        groups = Dict{String,Vector{Int}}()
        order = String[]
        for i in continuing
            node = lineages[i][depth + 1]
            haskey(groups, node) || push!(order, node)
            push!(get!(groups, node, Int[]), i)
        end
        for node in order
            push!(children, build(groups[node], depth + 1))
        end

        length(children) == 1 && return only(children)
        return (label=clade, height=1 / (2 * depth), children=children)
    end
    root = build(collect(eachindex(lineages)), 1)

    tips = length(lineages)
    tip_label = Vector{String}(undef, tips)
    node_label = String[]
    edges = Tuple{Int,Int}[]
    lengths = Float64[]
    counters = [0, tips]

    function visit(node, parent=nothing, parent_height=0.0)
        if haskey(node, :tip)
            counters[1] += 1
            number = counters[1]
            tip_label[number] = labels[node.tip]
            height = 0.0
        else
            counters[2] += 1
            number = counters[2]
            push!(node_label, node.label)
            height = node.height
        end
        if parent !== nothing
            push!(edges, (parent, number))
            push!(lengths, parent_height - height)
        end
        if !haskey(node, :tip)
            for child in node.children
                visit(child, number, height)
            end
        end
    end
    visit(root)

    edge = Matrix{Int}(undef, length(edges), 2)
    for (row, (parent, child)) in enumerate(edges)
        edge[row, 1] = parent
        edge[row, 2] = child
    end
    return TaxodistTree(edge, lengths, tip_label, node_label)
end
