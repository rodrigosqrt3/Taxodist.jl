"""
    taxo_distinctness(comm, d=nothing; metric="ranked", verbose=false, progress=true,
                      taxa=nothing, samples=nothing)

Compute Warwick and Clarke's taxonomic diversity and distinctness indices for
one or more samples: taxonomic diversity Δ, taxonomic distinctness Δ*, average
taxonomic distinctness Δ+, and variation in taxonomic distinctness Λ+.

`comm` is one of:

- a vector of taxon names (one sample, presence only);
- a dictionary, named tuple, or vector of `taxon => abundance` pairs (one
  sample);
- a `DataFrame` of abundances or presences with one row per sample and one
  column per taxon (an optional string column `sample` names the samples);
- a numeric or `Bool` matrix with one row per sample and one column per taxon,
  with the column (taxon) names given by `taxa` and, optionally, the sample
  names by `samples`.

`d` is a `TaxonomicDistanceMatrix` whose labels include every taxon present
in `comm`. When it is `nothing`, distances are computed with
`distance_matrix(taxa; metric, verbose, progress)`. The default `"ranked"`
metric gives the step weights of Clarke and Warwick (1 for congeneric taxa,
2 for confamilial, ...); any metric can be used.

Returns a named tuple `(metric, data)`: `metric` is the metric recorded in
`d` (`nothing` when unknown) and `data` a `DataFrame` with one row per sample
and columns `sample`, `S` (number of taxa present), `Delta`, `Delta_star`,
`Delta_plus`, `Lambda_plus`, and `SDelta_plus` (S Δ+). The indices are `NaN`
for samples with fewer than two taxa, or when a distance between present
taxa is missing or infinite.

With `x_i` the abundance of taxon `i`, `n = Σ x_i`, `S` the number of taxa
present, and `ω_ij` the distance between taxa `i` and `j`, sums running over
pairs `i < j`:

- `Δ = ΣΣ ω_ij x_i x_j / [n(n-1)/2]`
- `Δ* = ΣΣ ω_ij x_i x_j / ΣΣ x_i x_j`
- `Δ+ = ΣΣ ω_ij / [S(S-1)/2]`
- `Λ+ = ΣΣ (ω_ij - Δ+)² / [S(S-1)/2]`

With presence-only data, Δ and Δ* equal Δ+. Results agree with
`vegan::taxondive()` in R.

References: Warwick & Clarke (1995) Mar. Ecol. Prog. Ser. 129, 301-305;
Clarke & Warwick (1998) J. Appl. Ecol. 35, 523-531; Clarke & Warwick (2001)
Mar. Ecol. Prog. Ser. 216, 265-278.
"""
function taxo_distinctness(
    comm,
    d=nothing;
    metric="ranked",
    verbose::Bool=false,
    progress::Bool=true,
    taxa=nothing,
    samples=nothing,
)
    sample_names, taxon_names, abundances = _distinctness_samples(comm; taxa=taxa, samples=samples)
    if d === nothing
        d = distance_matrix(
            taxon_names; verbose=verbose, progress=progress, metric=_match_metric(metric),
        )
    end
    full = _distinctness_matrix(d, taxon_names)

    data = DataFrame(
        sample=String[],
        S=Int[],
        Delta=Float64[],
        Delta_star=Float64[],
        Delta_plus=Float64[],
        Lambda_plus=Float64[],
        SDelta_plus=Float64[],
    )
    for (row, name) in enumerate(sample_names)
        push!(data, merge((sample=name,), _distinctness_indices(abundances[row, :], full)))
    end
    return (metric=d.metric, data=data)
end

"""
    taxo_distinctness_test(sample, d; n_iter=999, rng=Random.default_rng())

Compare the average taxonomic distinctness Δ+ and its variation Λ+ of a
sample (a vector of taxon names) with their distribution in `n_iter` random
samples of the same number of taxa drawn without replacement from the master
list (species pool) of `d`, a `TaxonomicDistanceMatrix` whose labels include
every sampled taxon (Clarke & Warwick 1998, 2001).

Returns a named tuple with:

- `statistics`: a `DataFrame` with rows `Delta_plus` and `Lambda_plus` and
  columns `statistic`, `observed`, `expected`, `lower` and `upper` (the 2.5%
  and 97.5% quantiles of the random samples), `p_lower` and `p_upper`
  (one-sided permutation p-values `(1 + count) / (n_iter + 1)` for values as
  low or as high as observed);
- `simulated`: an `n_iter × 2` matrix of the simulated values (columns Δ+
  and Λ+);
- `S`, `pool_size`, `n_iter`, and `metric` (the metric recorded in `d`).

The expected Δ+ is the mean distance over all pairs of the master list,
which does not depend on the sample size; the expected Λ+ is the mean of the
simulated values. A Δ+ below the lower limit indicates that the sample is
less taxonomically distinct than expected from the pool. Pass `rng` (e.g.
`Random.Xoshiro(1)`) for reproducible results.
"""
function taxo_distinctness_test(
    sample,
    d;
    n_iter::Integer=999,
    rng::AbstractRNG=Random.default_rng(),
)
    d isa TaxonomicDistanceMatrix || throw(ArgumentError(
        "d must be a TaxonomicDistanceMatrix over the master list.",
    ))
    sample isa AbstractVector && all(taxon -> taxon isa AbstractString, sample) ||
        throw(ArgumentError("sample must be a vector of taxon names (strings)."))
    (n_iter isa Bool || n_iter < 1) && throw(ArgumentError(
        "n_iter must be a single positive integer.",
    ))
    sample = unique(String.(sample))
    full = d.values
    pool = d.taxa
    absent = setdiff(sample, pool)
    isempty(absent) || throw(ArgumentError(
        "Every sampled taxon must belong to the master list. " *
        "Not in d: $(join(repr.(absent), ", ")).",
    ))
    all(isfinite, _lower_triangle(full)) || throw(ArgumentError(
        "The master list distances must all be finite.",
    ))
    S = length(sample)
    pool_size = length(pool)
    (S < 2 || S > pool_size) && throw(ArgumentError(
        "The sample must contain between 2 and $(pool_size) taxa of the master list.",
    ))

    function stats_of(indices)
        omega = _lower_triangle(full[indices, indices])
        delta_plus = mean(omega)
        return (delta_plus, mean((omega .- delta_plus) .^ 2))
    end
    observed = stats_of([findfirst(==(taxon), pool) for taxon in sample])
    simulated = Matrix{Float64}(undef, n_iter, 2)
    for iteration in 1:n_iter
        values = stats_of(randperm(rng, pool_size)[1:S])
        simulated[iteration, 1] = values[1]
        simulated[iteration, 2] = values[2]
    end

    expected = (mean(_lower_triangle(full)), mean(simulated[:, 2]))
    tolerance = sqrt(eps(Float64))
    statistics = DataFrame(
        statistic=["Delta_plus", "Lambda_plus"],
        observed=[observed[1], observed[2]],
        expected=[expected[1], expected[2]],
        lower=[quantile(simulated[:, column], 0.025) for column in 1:2],
        upper=[quantile(simulated[:, column], 0.975) for column in 1:2],
        p_lower=[
            (1 + count(<=(observed[column] + tolerance), simulated[:, column])) / (n_iter + 1)
            for column in 1:2
        ],
        p_upper=[
            (1 + count(>=(observed[column] - tolerance), simulated[:, column])) / (n_iter + 1)
            for column in 1:2
        ],
    )
    return (
        statistics=statistics,
        simulated=simulated,
        S=S,
        pool_size=pool_size,
        n_iter=Int(n_iter),
        metric=d.metric,
    )
end

# Lower triangle of a square matrix in R's `x[lower.tri(x)]` order (column by
# column), so that sums are accumulated in the same order as in R.
function _lower_triangle(matrix::AbstractMatrix)
    n = size(matrix, 1)
    values = Vector{Float64}(undef, n * (n - 1) ÷ 2)
    position = 0
    for j in 1:n, i in (j + 1):n
        position += 1
        values[position] = matrix[i, j]
    end
    return values
end

_distinctness_name(name) = name isa Union{AbstractString,Symbol} ? String(name) : nothing

function _abundance_value(value)
    (value isa Real && !ismissing(value)) || throw(ArgumentError(
        "Abundances in comm must be non-negative and not missing.",
    ))
    abundance = Float64(value)
    (isnan(abundance) || abundance < 0) && throw(ArgumentError(
        "Abundances in comm must be non-negative and not missing.",
    ))
    return abundance
end

function _validated_taxon_names(names)
    all(name -> name !== nothing, names) && length(unique(names)) == length(names) ||
        throw(ArgumentError("The columns of comm must be named by unique taxa."))
    return String.(names)
end

"""
    _distinctness_samples(comm; taxa=nothing, samples=nothing)

Normalize the supported sample formats to `(sample_names, taxa, abundances)`
with an abundance matrix of one row per sample and one column per taxon.
"""
function _distinctness_samples(comm; taxa=nothing, samples=nothing)
    if comm isa AbstractVector{<:AbstractString} ||
       (comm isa AbstractVector && !isempty(comm) && all(item -> item isa AbstractString, comm))
        all(item -> !isempty(_trimws(item)), comm) || throw(ArgumentError(
            "comm cannot contain missing or empty taxon names.",
        ))
        unique_taxa = unique(String.(comm))
        return ["sample"], unique_taxa, ones(1, length(unique_taxa))
    end
    if comm isa AbstractVector && any(item -> item === nothing || ismissing(item), comm) &&
       all(item -> item === nothing || ismissing(item) || item isa AbstractString, comm)
        throw(ArgumentError("comm cannot contain missing or empty taxon names."))
    end
    if comm isa Union{AbstractDict,NamedTuple} || comm isa AbstractVector{<:Pair}
        entries = comm isa AbstractDict ? collect(pairs(comm)) :
            comm isa NamedTuple ? collect(pairs(comm)) : collect(comm)
        comm isa Dict && sort!(entries; by=entry -> string(first(entry)))
        entry_taxa = _validated_taxon_names([_distinctness_name(first(entry)) for entry in entries])
        abundances = reshape([_abundance_value(last(entry)) for entry in entries], 1, :)
        return ["sample"], entry_taxa, abundances
    end
    if comm isa AbstractVector && all(item -> item isa Real, comm)
        throw(ArgumentError(
            "A numeric comm vector must be named by taxon (use a Dict, a NamedTuple, " *
            "or `taxon => abundance` pairs).",
        ))
    end

    if comm isa AbstractDataFrame
        columns = names(comm)
        sample_names = if "sample" in columns && !(eltype(comm.sample) <: Union{Missing,Real})
            String.(string.(comm.sample))
        else
            nothing
        end
        sample_names === nothing || (columns = filter(!=("sample"), columns))
        taxon_names = _validated_taxon_names(columns)
        abundances = Matrix{Float64}(undef, nrow(comm), length(taxon_names))
        for (column, name) in enumerate(taxon_names), row in 1:nrow(comm)
            abundances[row, column] = _abundance_value(comm[row, name])
        end
        samples === nothing || (sample_names = String.(string.(samples)))
        return _sample_names(sample_names, nrow(comm)), taxon_names, abundances
    end
    if comm isa AbstractMatrix
        taxa === nothing && throw(ArgumentError(
            "The columns of comm must be named by unique taxa: pass `taxa` with the " *
            "taxon name of every column of a matrix.",
        ))
        length(taxa) == size(comm, 2) || throw(ArgumentError(
            "The columns of comm must be named by unique taxa: `taxa` must have one " *
            "name per column.",
        ))
        taxon_names = _validated_taxon_names([_distinctness_name(name) for name in taxa])
        abundances = [_abundance_value(value) for value in comm]
        sample_names = samples === nothing ? nothing : String.(string.(samples))
        return _sample_names(sample_names, size(comm, 1)), taxon_names, Matrix{Float64}(abundances)
    end
    throw(ArgumentError(
        "comm must be a vector of taxa, a dictionary (or named tuple, or pairs) of " *
        "abundances by taxon, or a numeric matrix or DataFrame.",
    ))
end

function _sample_names(sample_names, n_samples)
    if sample_names === nothing
        return n_samples == 1 ? ["sample"] : ["sample$(index)" for index in 1:n_samples]
    end
    length(sample_names) == n_samples || throw(ArgumentError(
        "There must be one sample name per row of comm.",
    ))
    return sample_names
end

# Full distance matrix restricted to and ordered by `taxa`.
function _distinctness_matrix(d, taxa)
    d isa TaxonomicDistanceMatrix || throw(ArgumentError(
        "d must be a TaxonomicDistanceMatrix.",
    ))
    absent = setdiff(taxa, d.taxa)
    isempty(absent) || throw(ArgumentError(
        "d has no distances for some taxa. Missing: $(join(repr.(absent), ", ")).",
    ))
    indices = [findfirst(==(taxon), d.taxa) for taxon in taxa]
    return d.values[indices, indices]
end

function _distinctness_indices(abundances::AbstractVector, full::AbstractMatrix)
    present = abundances .> 0
    S = count(present)
    empty = (
        S=S, Delta=NaN, Delta_star=NaN, Delta_plus=NaN, Lambda_plus=NaN, SDelta_plus=NaN,
    )
    S < 2 && return empty
    w = full[present, present]
    a = abundances[present]
    omega = _lower_triangle(w)
    all(isfinite, omega) || return empty

    xx = _lower_triangle(a * a')
    n = sum(a)
    weighted = sum(omega .* xx)
    delta_plus = mean(omega)
    return (
        S=S,
        Delta=weighted / (n * (n - 1) / 2),
        Delta_star=weighted / sum(xx),
        Delta_plus=delta_plus,
        Lambda_plus=mean((omega .- delta_plus) .^ 2),
        SDelta_plus=S * delta_plus,
    )
end
