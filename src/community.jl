# Community analyses: taxonomic dissimilarity between samples (Clarke,
# Somerfield & Chapman 2006) and standardized effect sizes of MPD and MNTD
# (Webb et al. 2002). Mirrors R's R/community.R.

"""
    taxo_dissimilarity(comm, d=nothing; index="gamma_plus", metric="ranked",
                       verbose=false, progress=true, taxa=nothing, samples=nothing)

Compute the taxonomic dissimilarity Γ+ (`index="gamma_plus"`, default) or Θ+
(`index="theta_plus"`) of Clarke, Somerfield and Chapman (2006) between every
pair of samples. Unlike Bray-Curtis or Jaccard dissimilarities, they credit
"near misses": two samples without any taxon in common are still similar when
their taxa are closely related.

`comm` holds the samples in any format accepted by `taxo_distinctness` (for
example a `DataFrame` or a matrix with `taxa`, one row per sample); values
greater than zero mean presence, abundances are not used. `d` is a
`TaxonomicDistanceMatrix` whose labels include every taxon of `comm`; when it
is `nothing`, distances are computed with `distance_matrix(taxa; metric)`.

Returns a `TaxonomicDistanceMatrix` between samples (labelled by sample name)
whose `method` is the index and whose `metric` is the metric of `d`
(`nothing` when unknown). Pairs involving a sample without taxa, or whose
nearest-relative distances are missing or infinite, are `NaN`.

With `ω_ij` the distance between taxon `i` of sample 1 and taxon `j` of
sample 2, and `s_1`, `s_2` the numbers of taxa in each sample:

- `Γ+ = [Σ_i min_j ω_ij + Σ_j min_i ω_ij] / (s_1 + s_2)`
- `Θ+ = [Σ_i min_j ω_ij / s_1 + Σ_j min_i ω_ij / s_2] / 2`

Taxa shared by both samples contribute zero. On a flat classification (every
pair of distinct taxa at the same distance) the indices reduce to the
presence/absence Bray-Curtis (Sorensen) and Kulczynski dissimilarities, times
that distance. Clarke et al. scale distances so that the largest possible
path length is 100; the indices are returned here in the units of `d`.

Reference: Clarke, Somerfield & Chapman (2006) J. Exp. Mar. Biol. Ecol. 330,
55-80.
"""
function taxo_dissimilarity(
    comm,
    d=nothing;
    index="gamma_plus",
    metric="ranked",
    verbose::Bool=false,
    progress::Bool=true,
    taxa=nothing,
    samples=nothing,
)
    (index isa AbstractString && index in ("gamma_plus", "theta_plus")) || throw(ArgumentError(
        "index must be \"gamma_plus\" or \"theta_plus\".",
    ))
    sample_names, taxon_names, abundances = _community_samples(comm; taxa=taxa, samples=samples)
    full, distance_metric = _community_distances(taxon_names, d, metric, verbose, progress)
    present = abundances .> 0
    n = length(sample_names)
    values = zeros(n, n)
    for l in 1:(n - 1), k in (l + 1):n
        value = _dissimilarity_pair(present[k, :], present[l, :], full, index)
        values[k, l] = value
        values[l, k] = value
    end
    return TaxonomicDistanceMatrix(values, sample_names; metric=distance_metric, method=index)
end

function _dissimilarity_pair(p::AbstractVector{Bool}, q::AbstractVector{Bool}, full, index)
    (any(p) && any(q)) || return NaN
    w = full[p, q]
    to_q = vec(minimum(w; dims=2))
    to_p = vec(minimum(w; dims=1))
    (all(isfinite, to_q) && all(isfinite, to_p)) || return NaN
    index == "gamma_plus" && return (sum(to_q) + sum(to_p)) / (length(to_q) + length(to_p))
    return (mean(to_q) + mean(to_p)) / 2
end

"""
    TaxodistSES

Result of `taxo_ses`, with fields `data` (a `DataFrame` with one row per
sample), `index`, `null_model`, `abundance_weighted`, `pool_size`, and
`metric`. Columns of `data` are also available as properties
(`result.NRI`).
"""
struct TaxodistSES
    data::DataFrame
    index::String
    null_model::String
    abundance_weighted::Bool
    pool_size::Int
    metric::Union{Nothing,String}
end

function Base.getproperty(result::TaxodistSES, name::Symbol)
    name in fieldnames(TaxodistSES) && return getfield(result, name)
    return getproperty(getfield(result, :data), name)
end

function Base.show(io::IO, ::MIME"text/plain", result::TaxodistSES)
    println(
        io,
        result.index == "mpd" ? "Net relatedness index (MPD)" : "Nearest taxon index (MNTD)",
    )
    println(
        io,
        "Null model: $(result.null_model), pool of $(result.pool_size) taxa",
        result.abundance_weighted ? ", abundance-weighted" : "",
        result.metric === nothing ? "" : "; metric: $(result.metric)",
    )
    show(io, MIME"text/plain"(), result.data; summary=false, eltypes=false, show_row_number=false)
end

"""
    taxo_ses(comm, d=nothing; index="mpd", null_model="pool", pool=nothing, n_iter=999,
             abundance_weighted=false, swaps=1000, metric="reciprocal", verbose=false,
             progress=true, taxa=nothing, samples=nothing, rng=Random.default_rng())

Compare the mean pairwise distance (MPD, `index="mpd"`) or mean
nearest-taxon distance (MNTD, `index="mntd"`) among the taxa of each sample
with its distribution under a null model, giving the standardized effect
sizes known as the net relatedness index (NRI) and the nearest taxon index
(NTI) of Webb et al. (2002), computed on taxonomic distances.

`comm` holds the samples in any format accepted by `taxo_distinctness` (a
vector of taxa, abundances by taxon, or a matrix or `DataFrame` with one row
per sample). `d` is a `TaxonomicDistanceMatrix` whose labels include every
taxon of `comm` and of `pool`; when it is `nothing`, it is computed with
`distance_matrix(taxa; metric)` (default `"reciprocal"`).

Null models (`null_model`):

- `"pool"` (default): each sample is replaced by the same number of taxa
  drawn at random, without replacement, from `pool`; with
  `abundance_weighted=true` the observed abundances are assigned to the drawn
  taxa. This is the null model of Webb et al. (2002). `pool` defaults to the
  labels of `d` (or, when `d` is `nothing`, the taxa of `comm`).
- `"independent_swap"`: the presence matrix of all samples is randomized with
  `swaps` independent swaps (Gotelli 2000), keeping the number of taxa of
  every sample and the number of samples in which every taxon occurs.
  Requires at least two samples.

With `abundance_weighted=true`, MPD weights each pair of taxa by the product
of their abundances (including each taxon with itself, which makes it Rao's
quadratic entropy) and MNTD weights each taxon by its abundance, as in the
picante package. Pass `rng` (e.g. `Random.Xoshiro(1)`) for reproducible
results.

Returns a `TaxodistSES` whose `data` has one row per sample and columns
`sample`, `S` (number of taxa), `observed`, `null_mean`, `null_sd`, `z`
(`(observed - null_mean) / null_sd`), `NRI` or `NTI` (`-z`), `p_lower` and
`p_upper` (one-sided permutation p-values for values as low or as high as
observed, `(1 + count) / (1 + n_iter)`), and `n_iter`. Samples with fewer than
two taxa have `NaN` statistics. Positive NRI and NTI indicate taxonomic
clustering; negative values, overdispersion.

References: Webb, Ackerly, McPeek & Donoghue (2002) Annu. Rev. Ecol. Syst.
33, 475-505; Gotelli (2000) Ecology 81, 2606-2621.
"""
function taxo_ses(
    comm,
    d=nothing;
    index="mpd",
    null_model="pool",
    pool=nothing,
    n_iter=999,
    abundance_weighted=false,
    swaps=1000,
    metric="reciprocal",
    verbose::Bool=false,
    progress::Bool=true,
    taxa=nothing,
    samples=nothing,
    rng::AbstractRNG=Random.default_rng(),
)
    (index isa AbstractString && index in ("mpd", "mntd")) || throw(ArgumentError(
        "index must be \"mpd\" or \"mntd\".",
    ))
    (null_model isa AbstractString && null_model in ("pool", "independent_swap")) || throw(ArgumentError(
        "null_model must be \"pool\" or \"independent_swap\".",
    ))
    _check_count(n_iter, "n_iter")
    _check_count(swaps, "swaps")
    abundance_weighted isa Bool || throw(ArgumentError(
        "abundance_weighted must be true or false.",
    ))
    n_iter = Int(n_iter)
    sample_names, taxon_names, abundances = _community_samples(comm; taxa=taxa, samples=samples)
    pool === nothing && (pool = d isa TaxonomicDistanceMatrix ? d.taxa : taxon_names)
    (pool isa AbstractVector && length(pool) >= 2 && all(taxon -> taxon isa AbstractString, pool)) ||
        throw(ArgumentError("pool must be a character vector of at least two taxa."))
    pool = unique(String.(pool))
    absent = setdiff(taxon_names, pool)
    if !isempty(absent) && null_model == "pool"
        throw(ArgumentError(
            "Every taxon of comm must belong to the pool. " *
            "Not in the pool: $(join(repr.(absent), ", ")).",
        ))
    end
    all_taxa = union(taxon_names, pool)
    full, distance_metric = _community_distances(all_taxa, d, metric, verbose, progress)
    n_samples = length(sample_names)
    null_model == "independent_swap" && n_samples < 2 && throw(ArgumentError(
        "The \"independent_swap\" null model needs at least two samples.",
    ))

    position = Dict(taxon => i for (i, taxon) in enumerate(all_taxa))
    comm_indices = [position[taxon] for taxon in taxon_names]
    pool_indices = [position[taxon] for taxon in pool]
    statistic(x, indices) =
        _community_statistic(x, full[indices, indices], index, abundance_weighted)

    observed = [statistic(abundances[k, :], comm_indices) for k in 1:n_samples]
    richness = [count(>(0), abundances[k, :]) for k in 1:n_samples]
    simulated = fill(NaN, n_samples, n_iter)
    for iteration in 1:n_iter
        if null_model == "pool"
            for k in 1:n_samples
                present = abundances[k, :] .> 0
                count(present) < 2 && continue
                drawn = randperm(rng, length(pool))[1:count(present)]
                simulated[k, iteration] = statistic(abundances[k, present], pool_indices[drawn])
            end
        else
            random = _independent_swap(abundances, Int(swaps), rng)
            for k in 1:n_samples
                simulated[k, iteration] = statistic(random[k, :], comm_indices)
            end
        end
    end

    null_mean = [mean(simulated[k, :]) for k in 1:n_samples]
    null_sd = [std(simulated[k, :]) for k in 1:n_samples]
    z = [isfinite(value) ? value : NaN for value in (observed .- null_mean) ./ null_sd]
    tolerance = sqrt(eps(Float64))
    function p_value(k, compare)
        row = simulated[k, :]
        (isnan(observed[k]) || any(isnan, row)) && return NaN
        return (1 + count(compare, row)) / (n_iter + 1)
    end
    p_lower = [p_value(k, <=(observed[k] + tolerance)) for k in 1:n_samples]
    p_upper = [p_value(k, >=(observed[k] - tolerance)) for k in 1:n_samples]

    data = DataFrame(
        sample=sample_names,
        S=richness,
        observed=observed,
        null_mean=null_mean,
        null_sd=null_sd,
        z=z,
    )
    data[!, index == "mpd" ? "NRI" : "NTI"] = -z
    data.p_lower = p_lower
    data.p_upper = p_upper
    data.n_iter = fill(n_iter, n_samples)
    return TaxodistSES(data, index, null_model, abundance_weighted, length(pool), distance_metric)
end

# -- Internal helpers ---------------------------------------------------------

function _check_count(value, name::AbstractString)
    (value isa Real && !(value isa Bool) && isfinite(value) && value >= 1 && isinteger(value)) ||
        throw(ArgumentError("$(name) must be a single positive integer."))
    return nothing
end

function _community_samples(comm; taxa=nothing, samples=nothing)
    parsed = _distinctness_samples(comm; taxa=taxa, samples=samples)
    isempty(parsed[2]) && throw(ArgumentError("comm contains no taxa."))
    return parsed
end

# Full distance matrix over `taxa`, computing it when `d` is `nothing`, and
# the metric recorded in `d`.
function _community_distances(taxa, d, metric, verbose::Bool, progress::Bool)
    if d === nothing
        d = distance_matrix(taxa; verbose=verbose, progress=progress, metric=_match_metric(metric))
    end
    return _distinctness_matrix(d, taxa), d.metric
end

"""
    _community_statistic(x, full, index, weighted)

MPD or MNTD of one sample (abundance vector `x` aligned with the rows of
`full`), following picante's definitions. `NaN` with fewer than two taxa or
non-finite distances.
"""
function _community_statistic(x::AbstractVector, full::AbstractMatrix, index::AbstractString, weighted::Bool)
    present = x .> 0
    count(present) < 2 && return NaN
    w = full[present, present]
    a = x[present]
    value = if index == "mpd"
        weighted ? sum(w .* (a * a')) / sum(a * a') : mean(_lower_triangle(w))
    else
        w[diagind(w)] .= Inf
        nearest = vec(minimum(w; dims=2))
        weighted ? sum(nearest .* a) / sum(a) : mean(nearest)
    end
    return isfinite(value) ? value : NaN
end

"""
    _independent_swap(samples, swaps, rng=Random.default_rng())

Independent swap randomization (Gotelli 2000): repeatedly picks two samples
and two taxa forming a checkerboard (each sample has one of the two taxa) and
swaps them, keeping every sample's richness and every taxon's number of
occurrences. Abundances move with the presences. Stops after `swaps`
successful swaps or `100 * swaps` attempts.
"""
function _independent_swap(samples::AbstractMatrix, swaps::Integer, rng::AbstractRNG=Random.default_rng())
    m = Matrix{Float64}(samples)
    n_row, n_col = size(m)
    done = 0
    attempts = 0
    while done < swaps && attempts < 100 * swaps
        attempts += 1
        r1, r2 = _two_distinct(rng, n_row)
        k1, k2 = _two_distinct(rng, n_col)
        a = m[r1, k1] > 0
        b = m[r1, k2] > 0
        c = m[r2, k1] > 0
        e = m[r2, k2] > 0
        if a && e && !b && !c
            m[r1, k2] = m[r1, k1]
            m[r1, k1] = 0
            m[r2, k1] = m[r2, k2]
            m[r2, k2] = 0
            done += 1
        elseif !a && !e && b && c
            m[r1, k1] = m[r1, k2]
            m[r1, k2] = 0
            m[r2, k2] = m[r2, k1]
            m[r2, k1] = 0
            done += 1
        end
    end
    return m
end

# Two distinct indices drawn at random from 1:n.
function _two_distinct(rng::AbstractRNG, n::Integer)
    first_index = rand(rng, 1:n)
    second_index = rand(rng, 1:(n - 1))
    second_index >= first_index && (second_index += 1)
    return first_index, second_index
end
