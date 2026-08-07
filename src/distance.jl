"""
    _compute_distance(lineage_a, lineage_b; taxon_a="A", taxon_b="B")

Compute the hierarchy distance between two resolved root-to-node lineages.
Shared ancestry is restricted to their continuous common prefix.
"""
function _compute_distance(
    lineage_a::AbstractVector{<:AbstractString},
    lineage_b::AbstractVector{<:AbstractString};
    taxon_a::AbstractString="A",
    taxon_b::AbstractString="B",
)
    depth_a = length(lineage_a)
    depth_b = length(lineage_b)
    mrca_depth = 0

    for (node_a, node_b) in zip(lineage_a, lineage_b)
        node_a == node_b || break
        mrca_depth += 1
    end

    if mrca_depth == 0
        return (
            distance=Inf,
            mrca=nothing,
            mrca_depth=0,
            depth_a=depth_a,
            depth_b=depth_b,
            taxon_a=String(taxon_a),
            taxon_b=String(taxon_b),
        )
    end

    return (
        distance=lineage_a == lineage_b ? 0.0 : 1.0 / mrca_depth,
        mrca=String(lineage_a[mrca_depth]),
        mrca_depth=mrca_depth,
        depth_a=depth_a,
        depth_b=depth_b,
        taxon_a=String(taxon_a),
        taxon_b=String(taxon_b),
    )
end

"""
    taxo_distance(taxon_a, taxon_b; verbose=false)

Retrieve two lineages and compute their taxonomic hierarchy distance.
"""
function taxo_distance(taxon_a, taxon_b; verbose::Bool=false)
    lineage_a = get_lineage(taxon_a; verbose=verbose)
    lineage_b = get_lineage(taxon_b; verbose=verbose)

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
        taxon_a=String(taxon_a),
        taxon_b=String(taxon_b),
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
    TaxonomicDistanceMatrix(values, taxa)

A square taxonomic-distance matrix that preserves taxon labels and supports
both integer indexing (`matrix[1, 2]`) and label indexing
(`matrix["Carnotaurus", "Tyrannosaurus"]`).
"""
struct TaxonomicDistanceMatrix <: AbstractMatrix{Float64}
    values::Matrix{Float64}
    taxa::Vector{String}

    function TaxonomicDistanceMatrix(
        values::AbstractMatrix{<:Real},
        taxa::AbstractVector{<:AbstractString},
    )
        size(values, 1) == size(values, 2) || throw(ArgumentError(
            "A taxonomic distance matrix must be square.",
        ))
        size(values, 1) == length(taxa) || throw(ArgumentError(
            "The number of taxon labels must match the matrix dimensions.",
        ))
        labels = String.(taxa)
        length(unique(labels)) == length(labels) || throw(ArgumentError(
            "Taxon labels must be unique.",
        ))
        return new(Matrix{Float64}(values), labels)
    end
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
    isempty(matrix.taxa) || println(io, "taxa: ", join(matrix.taxa, ", "))
    show(io, MIME"text/plain"(), matrix.values)
end

"""
    distance_matrix(taxa; verbose=false, progress=true)

Retrieve the requested lineages and calculate a symmetric labeled distance
matrix. Unresolved taxa receive `NaN` in their off-diagonal comparisons.
"""
function distance_matrix(
    taxa::AbstractVector;
    verbose::Bool=false,
    progress::Bool=true,
)
    labels = String.(taxa)
    length(unique(labels)) == length(labels) || throw(ArgumentError(
        "Taxon names must be unique.",
    ))

    n_taxa = length(labels)
    values = fill(NaN, n_taxa, n_taxa)
    for i in 1:n_taxa
        values[i, i] = 0.0
    end

    lineages = Vector{Union{Nothing,Vector{String}}}(undef, n_taxa)
    for (i, taxon) in pairs(labels)
        progress && verbose && println("Retrieving lineage $(i)/$(n_taxa): $(taxon)")
        lineage = get_lineage(taxon; verbose=verbose)
        lineages[i] = lineage === nothing ? nothing : String.(lineage)
    end

    if n_taxa >= 2
        for i in 1:(n_taxa - 1), j in (i + 1):n_taxa
            lineage_i = lineages[i]
            lineage_j = lineages[j]
            if lineage_i !== nothing && lineage_j !== nothing
                value = _compute_distance(
                    lineage_i,
                    lineage_j;
                    taxon_a=labels[i],
                    taxon_b=labels[j],
                ).distance
                values[i, j] = value
                values[j, i] = value
            end
        end
    end

    return TaxonomicDistanceMatrix(values, labels)
end

_distance_sort_key(value::Real) = isnan(value) ? Inf : Float64(value)

"""
    closest_relative(taxon, candidates; verbose=false)

Return candidate taxa ordered from the smallest hierarchy distance to the
largest. Unresolved candidates are retained at the end with distance `NaN`.
"""
function closest_relative(
    taxon,
    candidates::AbstractVector;
    verbose::Bool=false,
)
    focal = String(taxon)
    lineage = get_lineage(focal; verbose=verbose)
    if lineage === nothing
        @warn "Could not retrieve lineage for $(focal)"
        return nothing
    end

    result = DataFrame(taxon=String[], distance=Float64[])
    for candidate in String.(candidates)
        candidate_lineage = get_lineage(candidate; verbose=verbose)
        distance = candidate_lineage === nothing ? NaN : _compute_distance(
            lineage,
            candidate_lineage;
            taxon_a=focal,
            taxon_b=candidate,
        ).distance
        push!(result, (taxon=candidate, distance=distance))
    end

    order = sortperm(result.distance; by=_distance_sort_key)
    return result[order, :]
end

"""
    focal_distances(focal, taxa; verbose=false)

Calculate the distance and MRCA from one focal taxon to every requested taxon.
The return value is a named tuple containing the focal label and a sorted
`DataFrame` in its `data` field.
"""
function focal_distances(
    focal,
    taxa::AbstractVector;
    verbose::Bool=false,
)
    focal_name = String(focal)
    focal_lineage = get_lineage(focal_name; verbose=verbose)
    if focal_lineage === nothing
        @warn "Could not retrieve lineage for $(focal_name)"
        return nothing
    end

    result = DataFrame(
        taxon=String[],
        distance=Float64[],
        mrca=Union{Nothing,String}[],
        mrca_depth=Int[],
    )

    for taxon in String.(taxa)
        lineage = taxon == focal_name ? focal_lineage : get_lineage(taxon; verbose=verbose)
        if lineage === nothing
            push!(result, (
                taxon=taxon,
                distance=NaN,
                mrca=nothing,
                mrca_depth=0,
            ))
            continue
        end

        comparison = _compute_distance(
            focal_lineage,
            lineage;
            taxon_a=focal_name,
            taxon_b=taxon,
        )
        push!(result, (
            taxon=taxon,
            distance=comparison.distance,
            mrca=comparison.mrca,
            mrca_depth=comparison.mrca_depth,
        ))
    end

    order = sortperm(result.distance; by=_distance_sort_key)
    return (focal=focal_name, data=result[order, :])
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
Taxonomicon entry.
"""
function check_coverage(taxa::AbstractVector; verbose::Bool=false)
    result = DataFrame(taxon=String[], covered=Bool[])
    for taxon in String.(taxa)
        id = get_taxonomicon_id(taxon; verbose=verbose)
        push!(result, (taxon=taxon, covered=id !== nothing))
    end
    return result
end