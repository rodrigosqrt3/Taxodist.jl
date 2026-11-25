# -- Distance metrics --------------------------------------------------------

"""
    DISTANCE_METRICS

Names of the available distance metrics: `"reciprocal"` (default), `"path"`,
`"wu_palmer"`, and `"ranked"`. See `taxo_distance`.
"""
const DISTANCE_METRICS = ("reciprocal", "path", "wu_palmer", "ranked")

"""
    DEFAULT_RANK_LEVELS

Standard Linnaean ranks used by the `"ranked"` metric by default.
"""
const DEFAULT_RANK_LEVELS = (
    "kingdom", "phylum", "class", "order", "family", "genus", "species",
)

"""
    RANK_LEVELS

Rank levels used by the `"ranked"` distance metric (default
`DEFAULT_RANK_LEVELS`: kingdom, phylum, class, order, family, genus, species).
Levels are compared, in lower case, with the last word of each node's rank, so
`"Anamorph genus"` counts as a genus. Change them with, for example,
`Taxodist.RANK_LEVELS[] = ["kingdom", "phylum", "class", "order"]`.
"""
const RANK_LEVELS = Ref{Vector{String}}(collect(DEFAULT_RANK_LEVELS))

"""
    _match_metric(metric)

Validate a metric name, which must be a single string (a `Symbol` such as
`:path` is also accepted, as in 1.0), and return it as a `String`. Vectors
and other values are rejected. Messages follow R's `.match_metric()`.
"""
function _match_metric(metric)
    metric isa Union{AbstractString,Symbol} || throw(ArgumentError(
        "metric must be a single string, one of " *
        "$(join(repr.(DISTANCE_METRICS), ", ", ", or ")).",
    ))
    name = String(metric)
    name in DISTANCE_METRICS || throw(ArgumentError(
        "Unknown distance metric $(repr(name)). Use one of " *
        "$(join(repr.(DISTANCE_METRICS), ", ", ", or ")).",
    ))
    return name
end

"""
    _rank_levels()

Validated rank levels of the `"ranked"` metric from `RANK_LEVELS`, stripped
and in lower case.
"""
function _rank_levels()
    levels = RANK_LEVELS[]
    valid = levels isa AbstractVector && !isempty(levels) &&
            all(level -> level isa AbstractString, levels)
    valid || throw(ArgumentError(
        "Taxodist.RANK_LEVELS must be a non-empty vector of strings.",
    ))
    return String[lowercase(String(_trimws(level))) for level in levels]
end

# A lineage may carry the ranks of its nodes as `(nodes=..., ranks=...)`.
_lineage_nodes(lineage::NamedTuple) = lineage.nodes
_lineage_nodes(lineage) = lineage

function _has_ranks(lineage)
    lineage isa NamedTuple || return false
    ranks = lineage.ranks
    return ranks !== nothing && length(ranks) == length(lineage.nodes) &&
           !all(ismissing, ranks)
end

"""
    _project_lineage(nodes, ranks, levels)

Reduce a lineage to the nodes whose rank is one of `levels` (compared on the
lower-case last word of the rank), always keeping the root node, so that taxa
of different kingdoms remain connected, and the terminal node.
"""
function _project_lineage(nodes::AbstractVector, ranks::AbstractVector, levels)
    n_nodes = length(nodes)
    n_nodes == 0 && return String[]
    keep = [
        !ismissing(rank) && lowercase(replace(rank, r"^.*\s" => ""; count=1)) in levels
        for rank in ranks
    ]
    keep[1] = true
    keep[n_nodes] = true
    return String.(nodes[keep])
end

"""
    _metric_lineages(lineages, labels, metric, levels=_rank_levels())

Lineages on which a metric is computed: the plain nodes for depth-based
metrics; for `"ranked"`, the lineages reduced with `_project_lineage`.
Lineages without ranks cause an error for `"ranked"`.
"""
function _metric_lineages(lineages, labels, metric::AbstractString, levels=nothing)
    if metric != "ranked"
        return Union{Nothing,Vector{String}}[
            lineage === nothing ? nothing : String.(_lineage_nodes(lineage))
            for lineage in lineages
        ]
    end
    levels === nothing && (levels = _rank_levels())
    unranked = [lineage !== nothing && !_has_ranks(lineage) for lineage in lineages]
    if any(unranked)
        throw(ArgumentError(
            "The \"ranked\" metric requires the rank of every lineage node. " *
            "No ranks are available for $(join(repr.(String.(labels[unranked])), ", ")). " *
            "Supply `ranks` to taxo_from_lineages() or use another metric.",
        ))
    end
    return Union{Nothing,Vector{String}}[
        lineage === nothing ? nothing : _project_lineage(lineage.nodes, lineage.ranks, levels)
        for lineage in lineages
    ]
end

"""
    _metric_value(h, depth_a, depth_b, metric)

Distance from the MRCA depth `h` and the lineage depths. Identical lineages
have distance 0; lineages without a shared root have infinite distance, except
under `"wu_palmer"`, where it is 1.
"""
function _metric_value(h::Integer, depth_a::Integer, depth_b::Integer, metric::Symbol)
    h > 0 && h == depth_a && h == depth_b && return 0.0
    if metric === :reciprocal
        return 1.0 / h
    elseif metric === :path
        return h == 0 ? Inf : Float64(depth_a + depth_b - 2h)
    elseif metric === :wu_palmer
        return 1.0 - (2.0 * h) / (depth_a + depth_b)
    elseif metric === :ranked
        return h == 0 ? Inf : (depth_a + depth_b - 2h) / 2
    end
    throw(ArgumentError("Unknown distance metric $(repr(String(metric)))."))
end

_metric_value(h, depth_a, depth_b, metric::AbstractString) =
    _metric_value(h, depth_a, depth_b, Symbol(metric))

# -- Pairwise distances ------------------------------------------------------

function _common_prefix_depth(lineage_a, lineage_b)
    depth = 0
    for (node_a, node_b) in zip(lineage_a, lineage_b)
        node_a == node_b || break
        depth += 1
    end
    return depth
end

function _distance_value(lineage_a, lineage_b; metric="reciprocal")
    return _metric_value(
        _common_prefix_depth(lineage_a, lineage_b),
        length(lineage_a),
        length(lineage_b),
        Symbol(_match_metric(metric)),
    )
end

"""
    _compute_distance(lineage_a, lineage_b; taxon_a="A", taxon_b="B", metric="reciprocal")

Compute the hierarchy distance between two resolved root-to-node lineages.
Shared ancestry is restricted to their continuous common prefix. For the
`"ranked"` metric the lineages must carry their ranks, as
`(nodes=..., ranks=...)`; the MRCA and depths always refer to the complete
lineages.
"""
function _compute_distance(
    lineage_a,
    lineage_b;
    taxon_a::AbstractString="A",
    taxon_b::AbstractString="B",
    metric="reciprocal",
)
    metric = _match_metric(metric)
    nodes_a = _lineage_nodes(lineage_a)
    nodes_b = _lineage_nodes(lineage_b)
    depth_a = length(nodes_a)
    depth_b = length(nodes_b)
    mrca_depth = _common_prefix_depth(nodes_a, nodes_b)

    distance = if metric == "ranked"
        reduced = _metric_lineages([lineage_a, lineage_b], [taxon_a, taxon_b], metric)
        _metric_value(
            _common_prefix_depth(reduced[1], reduced[2]),
            length(reduced[1]),
            length(reduced[2]),
            :ranked,
        )
    else
        _metric_value(mrca_depth, depth_a, depth_b, Symbol(metric))
    end

    return (
        distance=distance,
        mrca=mrca_depth == 0 ? nothing : String(nodes_a[mrca_depth]),
        mrca_depth=mrca_depth,
        depth_a=depth_a,
        depth_b=depth_b,
        taxon_a=_taxon_label(taxon_a),
        taxon_b=_taxon_label(taxon_b),
        metric=metric,
    )
end

# Ranked lineage `(nodes, ranks)` of a taxon, from lineage_ranks().
function _ranked_lineage(taxon; verbose::Bool=false)
    ranked = lineage_ranks(taxon; verbose=verbose)
    ranked === nothing && return nothing
    return (nodes=String.(ranked.node), ranks=_rank_vector(ranked.rank))
end

# Lineage retrieval of the pairwise functions: plain lineages for depth-based
# metrics, ranked lineages for "ranked".
function _lineage_for_metric(taxon, metric::AbstractString; verbose::Bool=false)
    metric == "ranked" && return _ranked_lineage(taxon; verbose=verbose)
    return get_lineage(taxon; verbose=verbose)
end

"""
    taxo_distance(taxon_a, taxon_b; verbose=false, metric="reciprocal")

Retrieve two lineages from The Taxonomicon and compute their taxonomic
hierarchy distance. Returns a named tuple with `distance`, `mrca`,
`mrca_depth`, `depth_a`, `depth_b`, `taxon_a`, `taxon_b`, and `metric`, or
`nothing` when either taxon cannot be found. For the Open Tree or NCBI
taxonomies, resolve the taxa with `taxo_resolve(taxa; source="ott")` (or
`"ncbi"`) and pass the resolution to [`distance_matrix`](@ref).

# Metrics

With `h` the MRCA depth and `d_A`, `d_B` the lineage depths (numbers of nodes,
root included), every metric is zero for identical lineages:

- `"reciprocal"` (default): `1/h`. An ultrametric that weights shallow
  divergences heavily and compresses deep ones.
- `"path"`: `d_A + d_B - 2h`, the number of edges between the two nodes of the
  classification tree.
- `"wu_palmer"`: `1 - 2h/(d_A + d_B)`, the complement of the Wu-Palmer
  similarity, in `[0, 1)`; it normalizes for the different resolution of
  different parts of the classification.
- `"ranked"`: half the path length after reducing both lineages to the root,
  the nodes of the standard Linnaean ranks (see `RANK_LEVELS`), and the
  terminal node; for taxa of the same rank, the step weights of Clarke and
  Warwick's taxonomic distinctness (1 for congeneric species, 2 for
  confamilial species, ...). Requires lineage ranks (see `lineage_ranks`).

Lineages without a shared root have infinite distance, except under
`"wu_palmer"`, where it is 1.
"""
function taxo_distance(taxon_a, taxon_b; verbose::Bool=false, metric="reciprocal")
    metric = _match_metric(metric)
    lineage_a = _lineage_for_metric(taxon_a, metric; verbose=verbose)
    lineage_b = _lineage_for_metric(taxon_b, metric; verbose=verbose)

    if lineage_a === nothing
        @warn "Could not retrieve lineage for $(taxon_a)"
        return nothing
    end
    if lineage_b === nothing
        @warn "Could not retrieve lineage for $(taxon_b)"
        return nothing
    end

    return _compute_distance(
        lineage_a,
        lineage_b;
        taxon_a=_taxon_label(taxon_a),
        taxon_b=_taxon_label(taxon_b),
        metric=metric,
    )
end

"""
    mrca(taxon_a, taxon_b; verbose=false)

Return the most recent common ancestor of two taxa, or `nothing` when either
lineage cannot be resolved or the lineages have no common root.
"""
function mrca(taxon_a, taxon_b; verbose::Bool=false)
    result = taxo_distance(taxon_a, taxon_b; verbose=verbose)
    result === nothing && return nothing
    return result.mrca
end

"""
    TaxonomicDistanceMatrix(values, taxa; metric=nothing, method=nothing)

A square taxonomic-distance matrix that preserves taxon labels and supports
both integer indexing (`matrix[1, 2]`) and label indexing
(`matrix["Carnotaurus", "Tyrannosaurus"]`). The `metric` field records the
distance metric (`nothing` when unknown), like the `"method"` attribute of
R's `dist` objects. The `method` field is `nothing` for distances between
taxa; for dissimilarities between samples (`taxo_dissimilarity`), whose
labels are sample names, it is the index (`"gamma_plus"` or `"theta_plus"`)
and `metric` is the metric of the underlying taxon distances.
"""
struct TaxonomicDistanceMatrix <: AbstractMatrix{Float64}
    values::Matrix{Float64}
    taxa::Vector{String}
    metric::Union{Nothing,String}
    method::Union{Nothing,String}

    function TaxonomicDistanceMatrix(
        values::AbstractMatrix{<:Real},
        taxa::AbstractVector{<:AbstractString};
        metric::Union{Nothing,AbstractString,Symbol}=nothing,
        method::Union{Nothing,AbstractString}=nothing,
    )
        size(values, 1) == size(values, 2) || throw(ArgumentError(
            "A taxonomic distance matrix must be square.",
        ))
        size(values, 1) == length(taxa) || throw(ArgumentError(
            "The number of taxon labels must match the matrix dimensions.",
        ))
        return new(
            Matrix{Float64}(values),
            String.(taxa),
            metric === nothing ? nothing : String(metric),
            method === nothing ? nothing : String(method),
        )
    end

    # Internal: wrap freshly computed values without copying them.
    TaxonomicDistanceMatrix(
        ::Val{:nocopy}, values::Matrix{Float64}, taxa::Vector{String}, metric::String,
    ) = new(values, taxa, metric, nothing)
end

Base.size(matrix::TaxonomicDistanceMatrix) = size(matrix.values)
Base.IndexStyle(::Type{TaxonomicDistanceMatrix}) = IndexCartesian()
Base.getindex(matrix::TaxonomicDistanceMatrix, i::Int, j::Int) = matrix.values[i, j]
Base.Matrix(matrix::TaxonomicDistanceMatrix) = copy(matrix.values)

function Base.getindex(
    matrix::TaxonomicDistanceMatrix,
    taxon_a::AbstractString,
    taxon_b::AbstractString,
)
    i = findfirst(==(String(taxon_a)), matrix.taxa)
    j = findfirst(==(String(taxon_b)), matrix.taxa)
    i === nothing && throw(KeyError(taxon_a))
    j === nothing && throw(KeyError(taxon_b))
    return matrix.values[i, j]
end

function Base.show(io::IO, ::MIME"text/plain", matrix::TaxonomicDistanceMatrix)
    println(io, "$(size(matrix, 1))×$(size(matrix, 2)) TaxonomicDistanceMatrix")
    matrix.method === nothing || println(io, "method: ", matrix.method)
    matrix.metric === nothing || println(io, "metric: ", matrix.metric)
    isempty(matrix.taxa) || println(io, "taxa: ", join(matrix.taxa, ", "))
    show(io, MIME"text/plain"(), matrix.values)
end

"""
    distance_matrix(taxa; verbose=false, progress=true, metric="reciprocal")

Retrieve the requested lineages and calculate a symmetric labeled distance
matrix with the given metric (see `taxo_distance`), recorded in its `metric`
field. `taxa` may be a vector of names (retrieved from The Taxonomicon), a
`TaxodistResolution` from any source (`taxo_resolve(taxa; source="ott")` or
`"ncbi"`) whose stored lineages (and ranks) are reused without retrieval, or a
`TaxodistBundle`
whose validated stored matrix is returned (`metric` is then ignored).
Unresolved taxa receive `NaN` in their off-diagonal comparisons.
"""
function distance_matrix(
    taxa::AbstractVector;
    verbose::Bool=false,
    progress::Bool=true,
    metric="reciprocal",
)
    return _distance_matrix_impl(
        taxa, _match_metric(metric), _rank_levels(); verbose=verbose, progress=progress,
    )
end

function distance_matrix(
    resolution::TaxodistResolution;
    verbose::Bool=false,
    progress::Bool=true,
    metric="reciprocal",
)
    return _distance_matrix_impl(
        resolution, _match_metric(metric), _rank_levels(); verbose=verbose, progress=progress,
    )
end

function _distance_matrix_impl(taxa, metric::AbstractString, levels; verbose=false, progress=true)
    labels, lineages = _matrix_input(
        taxa; verbose=verbose, progress=progress, ranks=metric == "ranked",
    )
    reduced = _metric_lineages(lineages, labels, metric, levels)
    return _distance_matrix_from_lineages(labels, reduced; metric=metric)
end

"""
    _prefix_node_ids(lineages)

Replace every lineage node by an integer identifying its complete path from
the root. Two lineages share the node at depth `d` if and only if their
identifiers at `d` are equal, so MRCA depths can be computed by comparing
integers instead of strings. Returns the identifier vectors and a mask of
missing (`nothing`) lineages, whose identifier vector is empty.
"""
function _prefix_node_ids(lineages)
    codes = Dict{Tuple{Int,String},Int}()
    ids = Vector{Vector{Int}}(undef, length(lineages))
    missing_mask = falses(length(lineages))
    for (index, lineage) in enumerate(lineages)
        if lineage === nothing
            missing_mask[index] = true
            ids[index] = Int[]
            continue
        end
        path = Vector{Int}(undef, length(lineage))
        parent = 0
        for (depth, node) in enumerate(lineage)
            parent = get!(codes, (parent, String(node)), length(codes) + 1)
            path[depth] = parent
        end
        ids[index] = path
    end
    return ids, missing_mask
end

"""
    _MrcaIndex

Index for MRCA-depth queries, mirroring R's `.mrca_index()`. Lineages are
sorted lexicographically by their node identifiers (a missing depth sorts
before every node, so a lineage precedes its descendants); in that order the
common prefix of two lineages is the minimum of the common prefixes of all
adjacent lineages between them, as in the LCP table of a suffix array. The
depths from one lineage to all others therefore follow from two running
minima, in O(n) per lineage instead of a node-by-node comparison of every
pair.

- `order`: lineage at each sorted position;
- `pos`: sorted position of each lineage;
- `adj`: common prefix of sorted neighbours (length `n - 1`);
- `lens`: number of nodes of each lineage.
"""
struct _MrcaIndex
    order::Vector{Int}
    pos::Vector{Int}
    adj::Vector{Int}
    lens::Vector{Int}
end

function _mrca_index(ids::AbstractVector{Vector{Int}})
    n = length(ids)
    lens = Int[length(path) for path in ids]
    # `isless` on integer vectors is lexicographic with a proper prefix first,
    # which equals R's order on identifiers padded with 0 (identifiers are
    # positive); the stable sort keeps duplicates in input order.
    order = sortperm(ids; alg=Base.Sort.DEFAULT_STABLE)
    adj = Vector{Int}(undef, max(n - 1, 0))
    @inbounds for k in 1:(n - 1)
        adj[k] = _common_prefix_depth(ids[order[k]], ids[order[k + 1]])
    end
    pos = Vector{Int}(undef, n)
    pos[order] = 1:n
    return _MrcaIndex(order, pos, adj, lens)
end

"""
    _mrca_depths_from(index, j)

MRCA depths from lineage `j` to every lineage of `index`, in input order. The
depth of a lineage with itself is its length.
"""
function _mrca_depths_from(index::_MrcaIndex, j::Integer)
    n = length(index.order)
    p = index.pos[j]
    sorted = Vector{Int}(undef, n)
    sorted[p] = index.lens[j]
    running = typemax(Int)
    @inbounds for q in (p + 1):n
        running = min(running, index.adj[q - 1])
        sorted[q] = running
    end
    running = typemax(Int)
    @inbounds for q in (p - 1):-1:1
        running = min(running, index.adj[q])
        sorted[q] = running
    end
    return sorted[index.pos]
end

function _distance_matrix_from_lineages(labels, lineages; metric::AbstractString="reciprocal")
    n_taxa = length(labels)
    length(lineages) == n_taxa || throw(ArgumentError(
        "The number of lineages must match the number of labels.",
    ))
    ids, missing_mask = _prefix_node_ids(lineages)
    values = _pairwise_values(_mrca_index(ids), missing_mask, Symbol(metric))
    return TaxonomicDistanceMatrix(Val(:nocopy), values, String.(labels), String(metric))
end

# Square matrix of distances from an MRCA index. Column `a` is filled from
# the running minima of `adj` to the right and to the left of the sorted
# position of lineage `a` (column-wise writes keep the fill cache-friendly).
# Unresolved lineages give NaN off the diagonal.
function _pairwise_values(index::_MrcaIndex, missing_mask::AbstractVector{Bool}, metric::Symbol)
    order, adj, lens = index.order, index.adj, index.lens
    n_taxa = length(order)
    values = Matrix{Float64}(undef, n_taxa, n_taxa)
    @inbounds for p in 1:n_taxa
        a = order[p]
        length_a = lens[a]
        missing_a = missing_mask[a]
        values[a, a] = 0.0
        running = typemax(Int)
        for q in (p + 1):n_taxa
            running = min(running, adj[q - 1])
            b = order[q]
            values[b, a] = missing_a || missing_mask[b] ? NaN :
                _metric_value(running, lens[b], length_a, metric)
        end
        running = typemax(Int)
        for q in (p - 1):-1:1
            running = min(running, adj[q])
            b = order[q]
            values[b, a] = missing_a || missing_mask[b] ? NaN :
                _metric_value(running, lens[b], length_a, metric)
        end
    end
    return values
end

"""
    TaxonomicCrossDistances(values, rows, columns)

A rectangular matrix of taxonomic distances between two sets of taxa,
returned by `cross_distances`. Supports integer indexing and label indexing
(`matrix["Alpha", "Beta"]`) with row labels from the first set and column
labels from the second.
"""
struct TaxonomicCrossDistances <: AbstractMatrix{Float64}
    values::Matrix{Float64}
    rows::Vector{String}
    columns::Vector{String}
end

Base.size(matrix::TaxonomicCrossDistances) = size(matrix.values)
Base.IndexStyle(::Type{TaxonomicCrossDistances}) = IndexCartesian()
Base.getindex(matrix::TaxonomicCrossDistances, i::Int, j::Int) = matrix.values[i, j]
Base.Matrix(matrix::TaxonomicCrossDistances) = copy(matrix.values)

function Base.getindex(
    matrix::TaxonomicCrossDistances,
    row::AbstractString,
    column::AbstractString,
)
    i = findfirst(==(String(row)), matrix.rows)
    j = findfirst(==(String(column)), matrix.columns)
    i === nothing && throw(KeyError(row))
    j === nothing && throw(KeyError(column))
    return matrix.values[i, j]
end

function Base.show(io::IO, ::MIME"text/plain", matrix::TaxonomicCrossDistances)
    println(io, "$(size(matrix, 1))×$(size(matrix, 2)) TaxonomicCrossDistances")
    println(io, "rows: ", join(matrix.rows, ", "))
    println(io, "columns: ", join(matrix.columns, ", "))
    show(io, MIME"text/plain"(), matrix.values)
end

"""
    _matrix_input(taxa; verbose=false, progress=true, ranks=false)

Labels and lineages of the supported inputs (a vector of names or IDs, a
`TaxodistResolution`, or a `TaxodistBundle`). Unresolved entries are
`nothing`. With `ranks=true`, lineages carry their ranks as
`(nodes=..., ranks=...)` when available.
"""
function _matrix_input(
    taxa::AbstractVector;
    verbose::Bool=false,
    progress::Bool=true,
    ranks::Bool=false,
)
    # Numbers are taxon IDs; `missing` (and `nothing`) entries, like R's NA,
    # keep their row with NaN distances.
    labels = String[_taxon_label(taxon) for taxon in taxa]
    n_taxa = length(labels)
    lineages = Vector{Any}(undef, n_taxa)
    for (i, taxon) in enumerate(taxa)
        progress && verbose && println("Retrieving lineage $(i)/$(n_taxa): $(labels[i])")
        if ranks
            lineages[i] = (taxon === nothing || ismissing(taxon)) ? nothing :
                _ranked_lineage(taxon; verbose=verbose)
        else
            lineage = get_lineage(taxon; verbose=verbose)
            lineages[i] = lineage === nothing ? nothing : String.(lineage)
        end
    end
    ranks && return labels, lineages
    return labels, Union{Nothing,Vector{String}}[lineage for lineage in lineages]
end

function _matrix_input(resolution::TaxodistResolution; ranks::Bool=false, kwargs...)
    issubset(Set(RESOLUTION_COLUMNS), Set(names(resolution.data))) || throw(ArgumentError(
        "Invalid TaxodistResolution object.",
    ))
    labels = String.(resolution.input)
    lineages = Union{Nothing,Vector{String}}[
        lineage === nothing ? nothing : String.(lineage)
        for lineage in resolution.lineage
    ]
    (ranks && "ranks" in names(resolution.data)) || return labels, lineages
    ranked = Vector{Any}(undef, length(lineages))
    for index in eachindex(lineages)
        lineage = lineages[index]
        rank_vector = resolution.ranks[index]
        ranked[index] = if lineage !== nothing && rank_vector !== nothing &&
                           length(rank_vector) == length(lineage)
            (nodes=lineage, ranks=_rank_vector(rank_vector))
        else
            lineage
        end
    end
    return labels, ranked
end

"""
    cross_distances(x, y; verbose=false, progress=true, metric="reciprocal")

Compute the rectangular matrix of taxonomic distances from every taxon in `x`
to every taxon in `y` (for example, community members to reference taxa)
without computing the full square matrix. `x` and `y` may be vectors of names
(retrieved from The Taxonomicon), `TaxodistResolution`s from any source, or
`TaxodistBundle`s; distances are only meaningful when both come from the same
source. `metric` is one of the metrics
of `taxo_distance`. Unresolved taxa receive `NaN`; taxa without a common
ancestor receive `Inf` (1 under `"wu_palmer"`).
"""
function cross_distances(x, y; verbose::Bool=false, progress::Bool=true, metric="reciprocal")
    metric = _match_metric(metric)
    need_ranks = metric == "ranked"
    rows, lineages_x = _matrix_input(x; verbose=verbose, progress=progress, ranks=need_ranks)
    columns, lineages_y = _matrix_input(y; verbose=verbose, progress=progress, ranks=need_ranks)
    reduced_x = _metric_lineages(lineages_x, rows, metric)
    reduced_y = _metric_lineages(lineages_y, columns, metric)
    ids, missing_mask = _prefix_node_ids(vcat(reduced_x, reduced_y))
    values = _cross_values(
        _mrca_index(ids), missing_mask, length(rows), length(columns), Symbol(metric),
    )
    return TaxonomicCrossDistances(values, rows, columns)
end

# Distances from the first `n_x` lineages of `index` (rows) to the next `n_y`
# (columns), one column per query of the index; unresolved taxa give NaN.
function _cross_values(index::_MrcaIndex, missing_mask::AbstractVector{Bool},
                       n_x::Int, n_y::Int, metric::Symbol)
    values = fill(NaN, n_x, n_y)
    for j in 1:n_y
        column = n_x + j
        missing_mask[column] && continue
        depths = _mrca_depths_from(index, column)
        length_column = index.lens[column]
        @inbounds for i in 1:n_x
            missing_mask[i] && continue
            values[i, j] = _metric_value(depths[i], index.lens[i], length_column, metric)
        end
    end
    return values
end

"""
    _distance_order(distances)

Stable ascending order of distances with unresolved (`NaN`) values after
`Inf`, like R's `order(distance, na.last = TRUE)`.
"""
_distance_order(distances) = sortperm(distances; alg=Base.Sort.DEFAULT_STABLE, lt=isless)

"""
    closest_relative(taxon, candidates; verbose=false, metric="reciprocal")

Return candidate taxa ordered from the smallest hierarchy distance to the
largest under `metric` (see `taxo_distance`). Unresolved candidates are
retained at the end (after taxa without a common ancestor, whose distance is
`Inf`) with distance `NaN`.
"""
function closest_relative(
    taxon,
    candidates::AbstractVector;
    verbose::Bool=false,
    metric="reciprocal",
)
    metric = _match_metric(metric)
    focal = _taxon_label(taxon)
    lineage = _lineage_for_metric(taxon, metric; verbose=verbose)
    if lineage === nothing
        @warn "Could not retrieve lineage for $(focal)"
        return nothing
    end

    result = DataFrame(taxon=String[], distance=Float64[])
    for candidate in candidates
        label = _taxon_label(candidate)
        candidate_lineage = _lineage_for_metric(candidate, metric; verbose=verbose)
        distance = candidate_lineage === nothing ? NaN : _compute_distance(
            lineage,
            candidate_lineage;
            taxon_a=focal,
            taxon_b=label,
            metric=metric,
        ).distance
        push!(result, (taxon=label, distance=distance))
    end

    return result[_distance_order(result.distance), :]
end

"""
    focal_distances(focal, taxa; verbose=false, progress=true, metric="reciprocal")

Calculate the distance (under `metric`, see `taxo_distance`) and MRCA from one
focal taxon to every requested taxon. The return value is a named tuple
containing the focal label and a `DataFrame` in its `data` field, sorted by
distance with unresolved taxa last; these have distance `NaN`, `mrca`
`nothing`, and `mrca_depth` `missing` (R's `NA`). With `progress` (and `verbose`),
one line is printed per taxon.
"""
function focal_distances(
    focal,
    taxa::AbstractVector;
    verbose::Bool=false,
    progress::Bool=true,
    metric="reciprocal",
)
    metric = _match_metric(metric)
    focal_name = _taxon_label(focal)
    focal_lineage = _lineage_for_metric(focal, metric; verbose=verbose)
    if focal_lineage === nothing
        @warn "Could not retrieve lineage for $(focal_name)"
        return nothing
    end

    result = DataFrame(
        taxon=String[],
        distance=Float64[],
        mrca=Union{Nothing,String}[],
        mrca_depth=Union{Missing,Int}[],
    )

    n_taxa = length(taxa)
    for (index, original) in enumerate(taxa)
        taxon = _taxon_label(original)
        progress && verbose && println("Computing focal distance $(index)/$(n_taxa): $(taxon)")
        lineage = isequal(original, focal) ? focal_lineage :
            _lineage_for_metric(original, metric; verbose=verbose)
        if lineage === nothing
            push!(result, (
                taxon=taxon,
                distance=NaN,
                mrca=nothing,
                mrca_depth=missing,
            ))
            continue
        end

        comparison = _compute_distance(
            focal_lineage,
            lineage;
            taxon_a=focal_name,
            taxon_b=taxon,
            metric=metric,
        )
        push!(result, (
            taxon=taxon,
            distance=comparison.distance,
            mrca=comparison.mrca,
            mrca_depth=comparison.mrca_depth,
        ))
    end

    return (focal=focal_name, data=result[_distance_order(result.distance), :])
end

"""
    lineage_depth(taxon; verbose=false)

Return the number of cleaned hierarchy nodes in a resolved lineage.
"""
function lineage_depth(taxon; verbose::Bool=false)
    lineage = get_lineage(taxon; verbose=verbose)
    return lineage === nothing ? nothing : length(lineage)
end

"""
    check_coverage(taxa; verbose=false)

Return a table indicating whether each taxon can be resolved to a biological
Taxonomicon entry. `covered` is `true` when found, `false` when absent, and
`missing` when the answer is unknown because the server could not be reached.
"""
function check_coverage(taxa::AbstractVector; verbose::Bool=false)
    result = DataFrame(taxon=String[], covered=Union{Missing,Bool}[])
    for taxon in _taxon_label.(taxa)
        status = _lookup_taxonomicon_id(taxon; verbose=verbose).status
        covered = status == "found" ? true : status == "not_found" ? false : missing
        push!(result, (taxon=taxon, covered=covered))
    end
    return result
end

function _analysis_matrix(input; kwargs...)
    if input isa TaxonomicDistanceMatrix
        return input
    end
    input isa Union{AbstractVector,TaxodistResolution} || throw(ArgumentError(
        "Expected a vector of taxon names, a TaxodistResolution, or a TaxonomicDistanceMatrix.",
    ))
    return distance_matrix(input; kwargs...)
end

function _invalid_matrix_reason(matrix::TaxonomicDistanceMatrix)
    values = matrix.values
    any(isnan, values) && return :nan
    any(value -> !isfinite(value), values) && return :infinite
    return nothing
end

# Clustering methods: the name accepted by taxo_cluster() (case-insensitive),
# the linkage symbol recorded in the result, and the method code of R's
# Fortran HCLUST routine (ward.D = 1, single = 2, complete = 3, average = 4,
# mcquitty = 5, median = 6, centroid = 7, ward.D2 = 8).
const _linkage_aliases = Dict(
    "ward.d" => (:ward_presquared, 1),
    "ward_presquared" => (:ward_presquared, 1),
    "single" => (:single, 2),
    "complete" => (:complete, 3),
    "average" => (:average, 4),
    "mcquitty" => (:mcquitty, 5),
    "median" => (:median, 6),
    "centroid" => (:centroid, 7),
    "ward.d2" => (:ward, 8),
    "ward" => (:ward, 8),
)

# Position of the pair (i, j), i < j, in R's `dist` storage (Fortran IOFFST).
@inline _ioffst(n::Int, i::Int, j::Int) = j + (i - 1) * n - (i * (i + 1)) ÷ 2

"""
    _hclust_merges(diss, n, iopt)

Line-by-line port of the Fortran routine `HCLUST` of R's `stats` package
(F. Murtagh; R Core Team), which `stats::hclust()` calls: the nearest-
neighbour agglomeration with the Lance-Williams update of method `iopt`.
`diss` holds the dissimilarities in `dist` storage and is overwritten.
Returns the agglomeration history `ia`, `ib` and the heights `crit`. Ties are
broken exactly as in R, so that merges and leaf order are identical.
"""
function _hclust_merges(diss::Vector{Float64}, n::Int, iopt::Int)
    inf = 1.0e300
    ia = zeros(Int, n)
    ib = zeros(Int, n)
    crit = zeros(Float64, n)
    membr = ones(Float64, n)
    nn = zeros(Int, n)
    disnn = zeros(Float64, n)
    flag = trues(n)
    im = 0
    jj = 0
    jm = 0
    ncl = n
    if iopt == 8  # Ward "D2": squared dissimilarities
        for i in eachindex(diss)
            diss[i] = diss[i] * diss[i]
        end
    end

    # Nearest neighbour of each point to its right.
    for i in 1:(n - 1)
        dmin = inf
        for j in (i + 1):n
            ind = _ioffst(n, i, j)
            if dmin > diss[ind]
                dmin = diss[ind]
                jm = j
            end
        end
        nn[i] = jm
        disnn[i] = dmin
    end

    is_ward = iopt == 1 || iopt == 8
    while true
        # Least dissimilarity in the list of nearest neighbours.
        dmin = inf
        for i in 1:(n - 1)
            if flag[i] && disnn[i] < dmin
                dmin = disnn[i]
                im = i
                jm = nn[i]
            end
        end
        ncl -= 1

        i2 = min(im, jm)
        j2 = max(im, jm)
        ia[n - ncl] = i2
        ib[n - ncl] = j2
        iopt == 8 && (dmin = sqrt(dmin))
        crit[n - ncl] = dmin
        flag[j2] = false

        # Dissimilarities from the new cluster.
        dmin = inf
        for k in 1:n
            (flag[k] && k != i2) || continue
            ind1 = i2 < k ? _ioffst(n, i2, k) : _ioffst(n, k, i2)
            ind2 = j2 < k ? _ioffst(n, j2, k) : _ioffst(n, k, j2)
            d12 = diss[_ioffst(n, i2, j2)]
            if is_ward
                diss[ind1] = (membr[i2] + membr[k]) * diss[ind1] +
                             (membr[j2] + membr[k]) * diss[ind2] - membr[k] * d12
                diss[ind1] = diss[ind1] / (membr[i2] + membr[j2] + membr[k])
            elseif iopt == 2
                diss[ind1] = min(diss[ind1], diss[ind2])
            elseif iopt == 3
                diss[ind1] = max(diss[ind1], diss[ind2])
            elseif iopt == 4
                diss[ind1] = (membr[i2] * diss[ind1] + membr[j2] * diss[ind2]) /
                             (membr[i2] + membr[j2])
            elseif iopt == 5
                diss[ind1] = (diss[ind1] + diss[ind2]) / 2
            elseif iopt == 6
                diss[ind1] = ((diss[ind1] + diss[ind2]) - d12 / 2) / 2
            else  # iopt == 7, unweighted centroid
                diss[ind1] = (membr[i2] * diss[ind1] + membr[j2] * diss[ind2] -
                              membr[i2] * membr[j2] * d12 / (membr[i2] + membr[j2])) /
                             (membr[i2] + membr[j2])
            end

            if i2 < k
                if diss[ind1] < dmin
                    dmin = diss[ind1]
                    jj = k
                end
            elseif diss[ind1] < disnn[k]
                # Fix for non-monotone methods (centroid, median): k may now
                # have i2 as its nearest neighbour.
                disnn[k] = diss[ind1]
                nn[k] = i2
            end
        end
        membr[i2] = membr[i2] + membr[j2]
        disnn[i2] = dmin
        nn[i2] = jj

        # Update the nearest neighbours that pointed to i2 or j2.
        for i in 1:(n - 1)
            if flag[i] && (nn[i] == i2 || nn[i] == j2)
                dmin = inf
                for j in (i + 1):n
                    if flag[j]
                        ind = _ioffst(n, i, j)
                        if diss[ind] < dmin
                            dmin = diss[ind]
                            jj = j
                        end
                    end
                end
                nn[i] = jj
                disnn[i] = dmin
            end
        end
        ncl > 1 || break
    end
    return ia, ib, crit
end

"""
    _hclust_order(ia, ib, n)

Port of the Fortran routine `HCASS2` of R's `stats` package: converts the
agglomeration history into R's `merge` matrix (singletons negative, clusters
by merge step) and computes the leaf `order` of the dendrogram.
"""
function _hclust_order(ia::Vector{Int}, ib::Vector{Int}, n::Int)
    iia = copy(ia)
    iib = copy(ib)
    for i in 1:(n - 2)
        k = min(ia[i], ib[i])
        for j in (i + 1):(n - 1)
            ia[j] == k && (iia[j] = -i)
            ib[j] == k && (iib[j] = -i)
        end
    end
    for i in 1:(n - 1)
        iia[i] = -iia[i]
        iib[i] = -iib[i]
    end
    for i in 1:(n - 1)
        if iia[i] > 0 && iib[i] < 0
            iia[i], iib[i] = iib[i], iia[i]
        end
        if iia[i] > 0 && iib[i] > 0
            iia[i], iib[i] = min(iia[i], iib[i]), max(iia[i], iib[i])
        end
    end

    iorder = zeros(Int, n)
    iorder[1] = iia[n - 1]
    iorder[2] = iib[n - 1]
    loc = 2
    for i in (n - 2):-1:1
        for j in 1:loc
            iorder[j] == i || continue
            # Replace cluster i by its two members.
            iorder[j] = iia[i]
            loc += 1
            for k in loc:-1:(j + 2)
                iorder[k] = iorder[k - 1]
            end
            iorder[j + 1] = iib[i]
            break
        end
    end
    merges = hcat(iia[1:(n - 1)], iib[1:(n - 1)])
    return merges, -iorder
end

"""
    _hclust(values, method)

Hierarchical clustering of the square dissimilarity matrix `values` (its
lower triangle is used, as `as.dist()` does in R) with the same algorithm as
R's `stats::hclust()`, so that heights, merges, and leaf order are identical
to R's for all eight methods, ties included. Returns a `Clustering.Hclust`.
"""
function _hclust(values::AbstractMatrix{<:Real}, method::AbstractString="complete")
    linkage, iopt = _clustering_method(method)
    n = size(values, 1)
    n >= 2 || throw(ArgumentError("At least two objects are required for clustering."))
    diss = Vector{Float64}(undef, n * (n - 1) ÷ 2)
    for i in 1:(n - 1), j in (i + 1):n
        diss[_ioffst(n, i, j)] = values[j, i]
    end
    ia, ib, crit = _hclust_merges(diss, n, iopt)
    merges, order = _hclust_order(ia, ib, n)
    return Clustering.Hclust(merges, crit[1:(n - 1)], order, linkage)
end

function _clustering_method(method)
    entry = method isa Union{AbstractString,Symbol} ?
        get(_linkage_aliases, lowercase(String(method)), nothing) : nothing
    entry === nothing && throw(ArgumentError(
        "Unsupported clustering method: $(repr(method)). Valid methods are R's " *
        "\"ward.D\", \"single\", \"complete\", \"average\", \"mcquitty\", \"median\", " *
        "\"centroid\", and \"ward.D2\"; \"ward\" is an alias of \"ward.D2\" and " *
        "\"ward_presquared\" of \"ward.D\".",
    ))
    return entry
end

"""
    taxo_cluster(taxa; method="average", kwargs...)

Compute a taxonomic distance matrix and perform hierarchical clustering.
`taxa` may also be an existing `TaxonomicDistanceMatrix`. Other keyword
arguments (`verbose`, `progress`, `metric`) are passed to `distance_matrix`;
the metric used is recorded in `dist.metric`. `method` is one of R's
`stats::hclust()` methods: `"ward.D"`, `"single"`, `"complete"`, `"average"`
(UPGMA, the default), `"mcquitty"`, `"median"`, `"centroid"`, or
`"ward.D2"`; `"ward"` is accepted for `"ward.D2"` and `"ward_presquared"` for
`"ward.D"`. The clustering uses a port of R's algorithm, so heights, merges,
and leaf order equal R's, also when distances are tied. `hclust` is a
`Clustering.Hclust` (fields `merges`, `heights`, `order`, `linkage`).
"""
function taxo_cluster(taxa; method="average", kwargs...)
    linkage, _ = _clustering_method(method)
    matrix = _analysis_matrix(taxa; kwargs...)

    reason = _invalid_matrix_reason(matrix)
    if reason === :nan
        @warn "Distance matrix contains NaN values (taxa not found or server offline). Clustering skipped."
        return (hclust=nothing, dist=matrix, method=linkage)
    elseif reason === :infinite
        @warn "Distance matrix contains infinite values (no shared ancestor). Clustering skipped."
        return (hclust=nothing, dist=matrix, method=linkage)
    elseif size(matrix, 1) < 2
        @warn "At least two taxa are required for clustering. Clustering skipped."
        return (hclust=nothing, dist=matrix, method=linkage)
    end

    return (hclust=_hclust(matrix.values, String(method)), dist=matrix, method=linkage)
end

function _pcoa(matrix::TaxonomicDistanceMatrix, k::Int)
    values = matrix.values
    n_taxa = size(values, 1)
    centering = Matrix{Float64}(I, n_taxa, n_taxa) .- (1.0 / n_taxa)
    gram = centering * (-0.5 .* values .^ 2) * centering
    decomposition = eigen(Symmetric(gram))
    order = sortperm(decomposition.values; rev=true)
    eigenvalues = decomposition.values[order]
    eigenvectors = decomposition.vectors[:, order]

    selected = eigenvalues[1:k]
    positive = findall(>(0.0), selected)
    if length(positive) < k
        @warn "Only $(length(positive)) of the first $(k) eigenvalues are positive."
    end

    retained_values = selected[positive]
    retained_vectors = eigenvectors[:, positive]
    coordinates = retained_vectors .* reshape(sqrt.(retained_values), 1, :)

    numerator = sum(retained_values)
    absolute_total = sum(abs, eigenvalues)
    positive_total = sum(max(value, 0.0) for value in eigenvalues)
    gof = [
        absolute_total > 0 ? numerator / absolute_total : NaN,
        positive_total > 0 ? numerator / positive_total : NaN,
    ]

    return coordinates, eigenvalues, gof
end

"""
    taxo_ordinate(taxa; k=2, kwargs...)

Apply principal coordinates analysis (classical multidimensional scaling) to
a taxonomic distance matrix. Other keyword arguments (`verbose`, `progress`,
`metric`) are passed to `distance_matrix`.
"""
function taxo_ordinate(taxa; k=2, kwargs...)
    matrix = _analysis_matrix(taxa; kwargs...)
    reason = _invalid_matrix_reason(matrix)
    if reason === :nan
        @warn "Distance matrix contains NaN values. Ordination skipped."
        return (points=nothing, dist=matrix, GOF=nothing, eig=nothing)
    elseif reason === :infinite
        @warn "Distance matrix contains infinite values (no shared ancestor). Ordination skipped."
        return (points=nothing, dist=matrix, GOF=nothing, eig=nothing)
    end

    n_taxa = size(matrix, 1)
    if n_taxa < 2
        @warn "At least two taxa are required for ordination. Ordination skipped."
        return (points=nothing, dist=matrix, GOF=nothing, eig=nothing)
    end

    valid_k = k isa Real && !(k isa Bool) && isfinite(k) && k >= 1 && k == floor(k)
    valid_k || throw(ArgumentError("`k` must be a single positive integer."))
    dimensions = Int(k)
    maximum_dimensions = n_taxa - 1
    if dimensions > maximum_dimensions
        @warn "`k` reduced from $(dimensions) to $(maximum_dimensions), the maximum for $(n_taxa) taxa."
        dimensions = maximum_dimensions
    end

    coordinates, eigenvalues, gof = _pcoa(matrix, dimensions)
    points = DataFrame(taxon=copy(matrix.taxa))
    for axis in axes(coordinates, 2)
        points[!, Symbol("PC$(axis)")] = coordinates[:, axis]
    end

    return (points=points, dist=matrix, GOF=gof, eig=eigenvalues)
end
