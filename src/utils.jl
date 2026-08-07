"""
    compare_lineages(taxon_a, taxon_b; verbose=false)

Print the shared trunk and the two lineage branches after their MRCA. Returns
the complete lineages and MRCA depth for programmatic use.
"""
function compare_lineages(taxon_a, taxon_b; verbose::Bool=false)
    lineage_a = get_lineage(taxon_a; verbose=verbose)
    lineage_b = get_lineage(taxon_b; verbose=verbose)

    if lineage_a === nothing || lineage_b === nothing
        println("Error: Could not retrieve one or both lineages")
        return nothing
    end

    result = _compute_distance(
        lineage_a,
        lineage_b;
        taxon_a=String(taxon_a),
        taxon_b=String(taxon_b),
    )
    mrca_depth = result.mrca_depth

    println("Lineage Comparison")
    println("MRCA: $(result.mrca) at depth $(mrca_depth)")

    if mrca_depth > 0
        println("\nShared lineage ($(mrca_depth) nodes):")
        for node in lineage_a[1:mrca_depth]
            println("  $(node)")
        end
    end

    println("\n$(taxon_a) only ($(length(lineage_a) - mrca_depth) nodes):")
    if length(lineage_a) > mrca_depth
        for node in lineage_a[(mrca_depth + 1):end]
            println("> $(node)")
        end
    end

    println("\n$(taxon_b) only ($(length(lineage_b) - mrca_depth) nodes):")
    if length(lineage_b) > mrca_depth
        for node in lineage_b[(mrca_depth + 1):end]
            println("> $(node)")
        end
    end

    return (
        lineage_a=String.(lineage_a),
        lineage_b=String.(lineage_b),
        mrca_depth=mrca_depth,
    )
end

"""
    shared_clades(taxon_a, taxon_b; verbose=false)

Return the continuous shared lineage from the root through the MRCA. Returns
an empty vector for disconnected hierarchies and `nothing` when either taxon
cannot be resolved.
"""
function shared_clades(taxon_a, taxon_b; verbose::Bool=false)
    lineage_a = get_lineage(taxon_a; verbose=verbose)
    lineage_b = get_lineage(taxon_b; verbose=verbose)
    (lineage_a === nothing || lineage_b === nothing) && return nothing

    result = _compute_distance(
        lineage_a,
        lineage_b;
        taxon_a=String(taxon_a),
        taxon_b=String(taxon_b),
    )
    result.mrca_depth == 0 && return String[]
    return String.(lineage_a[1:result.mrca_depth])
end

"""
    is_member(taxon, clade; verbose=false)

Test whether the complete clade name occurs in a taxon's lineage. Matching is
case-insensitive and ignores outer whitespace, but does not perform partial or
regular-expression matching.
"""
function is_member(taxon, clade; verbose::Bool=false)
    lineage = get_lineage(taxon; verbose=verbose)
    lineage === nothing && return nothing

    normalized_clade = lowercase(strip(String(clade)))
    normalized_lineage = lowercase.(strip.(String.(lineage)))
    return normalized_clade in normalized_lineage
end

"""
    filter_clade(taxa, clade; verbose=false)

Return only the taxa whose resolved lineage contains the requested clade.
Unresolved taxa are omitted.
"""
function filter_clade(
    taxa::AbstractVector,
    clade;
    verbose::Bool=false,
)
    kept = String[]
    for taxon in String.(taxa)
        membership = is_member(taxon, clade; verbose=verbose)
        membership === true && push!(kept, taxon)
    end
    return kept
end

"""
    taxo_path(taxon_a, taxon_b; verbose=false)

Return the node-by-node path ascending from taxon A to the MRCA and descending
to taxon B. The result contains `taxon_a`, `taxon_b`, and a `data` DataFrame
with columns `node`, `depth`, and `direction`.
"""
function taxo_path(taxon_a, taxon_b; verbose::Bool=false)
    name_a = String(taxon_a)
    name_b = String(taxon_b)
    lineage_a = get_lineage(name_a; verbose=verbose)
    lineage_b = get_lineage(name_b; verbose=verbose)

    if lineage_a === nothing
        @warn "Could not retrieve lineage for $(name_a)"
        return nothing
    end
    if lineage_b === nothing
        @warn "Could not retrieve lineage for $(name_b)"
        return nothing
    end

    result = _compute_distance(
        lineage_a,
        lineage_b;
        taxon_a=name_a,
        taxon_b=name_b,
    )
    mrca_depth = result.mrca_depth
    if mrca_depth == 0
        @warn "No common ancestor found for $(repr(name_a)) and $(repr(name_b))."
        return nothing
    end

    path = DataFrame(node=String[], depth=Int[], direction=String[])

    if mrca_depth < length(lineage_a)
        for index in length(lineage_a):-1:(mrca_depth + 1)
            push!(path, (
                node=String(lineage_a[index]),
                depth=index,
                direction="a",
            ))
        end
    end

    push!(path, (
        node=String(result.mrca),
        depth=mrca_depth,
        direction="mrca",
    ))

    if mrca_depth < length(lineage_b)
        for index in (mrca_depth + 1):length(lineage_b)
            push!(path, (
                node=String(lineage_b[index]),
                depth=index,
                direction="b",
            ))
        end
    end

    return (taxon_a=name_a, taxon_b=name_b, data=path)
end

"""
    summary_taxodist_ord(ordination)

Return a table of retained PCoA eigenvalues, variance percentages, and
cumulative percentages.
"""
function summary_taxodist_ord(ordination)
    ordination.points === nothing && return nothing
    ordination.eig === nothing && return nothing

    n_axes = max(ncol(ordination.points) - 1, 0)
    n_axes == 0 && return DataFrame(
        Axis=String[],
        Eigenvalue=Float64[],
        Variance_Pct=Float64[],
        Cumulative_Pct=Float64[],
    )

    positive_eigenvalues = ordination.eig[ordination.eig .> 0]
    total = sum(positive_eigenvalues)
    retained = positive_eigenvalues[1:min(n_axes, length(positive_eigenvalues))]
    variance = total > 0 ? 100 .* retained ./ total : fill(NaN, length(retained))

    return DataFrame(
        Axis=["PC$(axis)" for axis in eachindex(retained)],
        Eigenvalue=retained,
        Variance_Pct=variance,
        Cumulative_Pct=cumsum(variance),
    )
end

"""
    plot_taxodist_ord(ordination; kwargs...)

Create a labeled scatter plot of the first two PCoA axes.
"""
function plot_taxodist_ord(ordination; kwargs...)
    ordination.points === nothing && return nothing
    points = ordination.points
    n_axes = ncol(points) - 1
    n_axes == 0 && return nothing

    x = points.PC1
    y = n_axes >= 2 ? points.PC2 : zeros(length(x))
    goodness = ordination.GOF === nothing ? NaN : ordination.GOF[1]
    title = isfinite(goodness) ?
        "Taxonomic Ordination (PCoA)  (GOF = $(round(goodness; digits=3)))" :
        "Taxonomic Ordination (PCoA)"

    plot = Plots.scatter(
        x,
        y;
        series_annotations=points.taxon,
        legend=false,
        xlabel="PC1",
        ylabel="PC2",
        title=title,
        kwargs...,
    )
    return plot
end

function _dendrogram_coordinates(hclust_result)
    leaf_positions = Dict(
        leaf => Float64(position)
        for (position, leaf) in enumerate(hclust_result.order)
    )
    node_positions = Dict{Int,Float64}()
    node_heights = Dict{Int,Float64}()
    xs = Float64[]
    ys = Float64[]

    child_position(id) = id < 0 ? leaf_positions[-id] : node_positions[id]
    child_height(id) = id < 0 ? 0.0 : node_heights[id]

    for row in axes(hclust_result.merges, 1)
        left = hclust_result.merges[row, 1]
        right = hclust_result.merges[row, 2]
        height = Float64(hclust_result.heights[row])
        left_x = child_position(left)
        right_x = child_position(right)
        left_height = child_height(left)
        right_height = child_height(right)

        append!(xs, [left_x, left_x, right_x, right_x, NaN])
        append!(ys, [left_height, height, height, right_height, NaN])
        node_positions[row] = (left_x + right_x) / 2
        node_heights[row] = height
    end

    return xs, ys
end

"""
    plot_taxodist_cluster(cluster; kwargs...)

Create a labeled dendrogram from a result returned by `taxo_cluster`.
"""
function plot_taxodist_cluster(cluster; kwargs...)
    cluster.hclust === nothing && return nothing
    xs, ys = _dendrogram_coordinates(cluster.hclust)
    labels = cluster.dist.taxa[cluster.hclust.order]
    return Plots.plot(
        xs,
        ys;
        label=false,
        xticks=(1:length(labels), labels),
        xrotation=45,
        xlabel="",
        ylabel="Distance",
        title="Taxonomic Clustering",
        kwargs...,
    )
end

"""
    taxo_heatmap(taxa; display_plot=true, kwargs...)

Plot a hierarchically ordered taxonomic-distance heatmap and return the
underlying `TaxonomicDistanceMatrix`.
"""
function taxo_heatmap(taxa; display_plot::Bool=true, kwargs...)
    matrix = _analysis_matrix(taxa; kwargs...)
    reason = _invalid_matrix_reason(matrix)
    if reason === :nan
        @warn "Distance matrix contains NaN values. Heatmap skipped."
        return matrix
    elseif reason === :infinite
        @warn "Distance matrix contains infinite values (no shared ancestor). Heatmap skipped."
        return matrix
    elseif size(matrix, 1) < 2
        @warn "At least two taxa are required for a heatmap. Heatmap skipped."
        return matrix
    end

    cluster = Clustering.hclust(Matrix(matrix); linkage=:average, branchorder=:r)
    order = cluster.order
    labels = matrix.taxa[order]
    values = matrix.values[order, order]
    plot = Plots.heatmap(
        labels,
        labels,
        values;
        aspect_ratio=1,
        color=:viridis,
        xrotation=45,
        title="Taxonomic Distance Heatmap",
    )
    display_plot && display(plot)
    return matrix
end