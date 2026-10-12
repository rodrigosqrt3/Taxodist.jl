function _common_prefix_depth(lineage_a, lineage_b)
    depth = 0
    for (node_a, node_b) in zip(lineage_a, lineage_b)
        node_a == node_b || break
        depth += 1
    end
    return depth
end

function _distance_value(lineage_a, lineage_b)
    mrca_depth = _common_prefix_depth(lineage_a, lineage_b)
    mrca_depth == 0 && return Inf
    length(lineage_a) == length(lineage_b) == mrca_depth && return 0.0
    return 1.0 / mrca_depth
end

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
    mrca_depth = _common_prefix_depth(lineage_a, lineage_b)

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
        distance=(depth_a == depth_b == mrca_depth) ? 0.0 : 1.0 / mrca_depth,
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
        return new(Matrix{Float64}(values), String.(taxa))
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
matrix. `taxa` may be a vector of names, a `TaxodistResolution` whose stored
lineages are reused without retrieval, or a `TaxodistBundle` whose validated
stored matrix is returned. Unresolved taxa receive `NaN` in their
off-diagonal comparisons.
"""
function distance_matrix(
    taxa::AbstractVector;
    verbose::Bool=false,
    progress::Bool=true,
)
    labels, lineages = _matrix_input(taxa; verbose=verbose, progress=progress)
    return _distance_matrix_from_lineages(labels, lineages)
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

function _prefix_distance(ids_a::Vector{Int}, ids_b::Vector{Int})
    common = min(length(ids_a), length(ids_b))
    depth = 0
    @inbounds while depth < common && ids_a[depth + 1] == ids_b[depth + 1]
        depth += 1
    end
    depth == 0 && return Inf
    length(ids_a) == length(ids_b) == depth && return 0.0
    return 1.0 / depth
end

function _distance_matrix_from_lineages(labels, lineages)
    n_taxa = length(labels)
    length(lineages) == n_taxa || throw(ArgumentError(
        "The number of lineages must match the number of labels.",
    ))
    ids, missing_mask = _prefix_node_ids(lineages)
    values = Matrix{Float64}(undef, n_taxa, n_taxa)
    @inbounds for j in 1:n_taxa
        values[j, j] = 0.0
        for i in (j + 1):n_taxa
            value = missing_mask[i] || missing_mask[j] ? NaN : _prefix_distance(ids[i], ids[j])
            values[i, j] = value
            values[j, i] = value
        end
    end

    return TaxonomicDistanceMatrix(values, labels)
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

function _matrix_input(taxa::AbstractVector; verbose::Bool=false, progress::Bool=true)
    labels = String.(taxa)
    n_taxa = length(labels)
    lineages = Vector{Union{Nothing,Vector{String}}}(undef, n_taxa)
    for (i, taxon) in pairs(labels)
        progress && verbose && println("Retrieving lineage $(i)/$(n_taxa): $(taxon)")
        lineage = get_lineage(taxon; verbose=verbose)
        lineages[i] = lineage === nothing ? nothing : String.(lineage)
    end
    return labels, lineages
end

function _matrix_input(resolution::TaxodistResolution; kwargs...)
    issubset(Set(RESOLUTION_COLUMNS), Set(names(resolution.data))) || throw(ArgumentError(
        "Invalid TaxodistResolution object.",
    ))
    labels = String.(resolution.input)
    lineages = Union{Nothing,Vector{String}}[
        lineage === nothing ? nothing : String.(lineage)
        for lineage in resolution.lineage
    ]
    return labels, lineages
end

"""
    cross_distances(x, y; verbose=false, progress=true)

Compute the rectangular matrix of taxonomic distances from every taxon in `x`
to every taxon in `y` (for example, community members to reference taxa)
without computing the full square matrix. `x` and `y` may be vectors of names,
`TaxodistResolution`s, or `TaxodistBundle`s. Unresolved taxa receive `NaN`;
taxa without a common ancestor receive `Inf`.
"""
function cross_distances(x, y; verbose::Bool=false, progress::Bool=true)
    rows, lineages_x = _matrix_input(x; verbose=verbose, progress=progress)
    columns, lineages_y = _matrix_input(y; verbose=verbose, progress=progress)
    ids, missing_mask = _prefix_node_ids(vcat(lineages_x, lineages_y))
    n_x = length(rows)
    values = Matrix{Float64}(undef, n_x, length(columns))
    for j in eachindex(columns), i in 1:n_x
        values[i, j] = missing_mask[i] || missing_mask[n_x + j] ? NaN :
            _prefix_distance(ids[i], ids[n_x + j])
    end
    return TaxonomicCrossDistances(values, rows, columns)
end

function distance_matrix(
    resolution::TaxodistResolution;
    verbose::Bool=false,
    progress::Bool=true,
)
    labels, lineages = _matrix_input(resolution)
    return _distance_matrix_from_lineages(labels, lineages)
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
Taxonomicon entry. `covered` is `true` when found, `false` when absent, and
`missing` when the answer is unknown because the server could not be reached.
"""
function check_coverage(taxa::AbstractVector; verbose::Bool=false)
    result = DataFrame(taxon=String[], covered=Union{Missing,Bool}[])
    for taxon in String.(taxa)
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

"""
    taxo_cluster(taxa; method="average", kwargs...)

Compute a taxonomic distance matrix and perform hierarchical clustering.
`taxa` may also be an existing `TaxonomicDistanceMatrix`.
"""
function taxo_cluster(taxa; method="average", kwargs...)
    matrix = _analysis_matrix(taxa; kwargs...)
    linkage = Symbol(lowercase(String(method)))
    valid_linkages = (:single, :average, :complete, :ward, :ward_presquared)
    linkage in valid_linkages || throw(ArgumentError(
        "Unsupported clustering method: $(repr(method)).",
    ))

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

    result = Clustering.hclust(
        Matrix(matrix);
        linkage=linkage,
        branchorder=:r,
    )
    return (hclust=result, dist=matrix, method=linkage)
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
a taxonomic distance matrix.
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