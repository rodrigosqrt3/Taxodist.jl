# Taxonomic dissimilarity (Clarke, Somerfield & Chapman 2006) and
# standardized effect sizes of MPD and MNTD (Webb et al. 2002).

using Taxodist: JSON3, mean
using Taxodist.Random: Xoshiro
using Taxodist.LinearAlgebra: I
using Taxodist.DataFrames: DataFrame, nrow

const COMMUNITY_LINEAGES = [
    "A" => ["Biota", "Animalia", "Chordata", "Mammalia", "Carnivora", "A"],
    "B" => ["Biota", "Animalia", "Chordata", "Mammalia", "Carnivora", "B"],
    "C" => ["Biota", "Animalia", "Chordata", "Mammalia", "Primates", "C"],
    "D" => ["Biota", "Animalia", "Arthropoda", "Insecta", "Diptera", "D"],
    "E" => ["Biota", "Animalia", "Arthropoda", "Insecta", "Coleoptera", "E"],
    "F" => ["Biota", "Plantae", "Tracheophyta", "Magnoliopsida", "Fagales", "F"],
]
community_d(metric="path") =
    distance_matrix(taxo_from_lineages(COMMUNITY_LINEAGES); metric=metric, progress=false)

# Abundance matrix with one row per sample, from rows of `taxon => value`.
function community_frame(rows::Pair...)
    taxa = unique(vcat([first.(last(row)) for row in rows]...))
    frame = DataFrame(sample=[first(row) for row in rows])
    for taxon in taxa
        frame[!, taxon] = [Float64(get(Dict(last(row)), taxon, 0)) for row in rows]
    end
    return frame
end

dissimilarity_values(matrix) = Taxodist._lower_triangle(matrix.values)

# Presence/absence references (vegan::vegdist(binary = TRUE)).
function sorensen(x, y)
    a, b, j = count(x), count(y), count(x .& y)
    return (a + b - 2j) / (a + b)
end
kulczynski(x, y) = 1 - (count(x .& y) / count(x) + count(x .& y) / count(y)) / 2

@testset "Gamma+ and Theta+ reduce to Sorensen and Kulczynski on a flat classification" begin
    rng = Xoshiro(11)
    taxa = ["sp$(i)" for i in 1:12]
    flat = TaxonomicDistanceMatrix(ones(12, 12) - I, taxa)
    presence = rand(rng, 6, 12) .< 0.4
    for row in 1:6
        any(presence[row, :]) || (presence[row, row] = true)
    end
    # Abundances are ignored.
    comm = presence .* rand(rng, 1:9, 6, 12)
    sample_labels = ["s$(i)" for i in 1:6]
    gamma = taxo_dissimilarity(comm, flat; taxa=taxa, samples=sample_labels)
    theta = taxo_dissimilarity(comm, flat; index="theta_plus", taxa=taxa, samples=sample_labels)
    expected_gamma = [sorensen(presence[k, :], presence[l, :]) for l in 1:5 for k in (l + 1):6]
    expected_theta = [kulczynski(presence[k, :], presence[l, :]) for l in 1:5 for k in (l + 1):6]
    @test dissimilarity_values(gamma) ≈ expected_gamma
    @test dissimilarity_values(theta) ≈ expected_theta
    @test gamma.taxa == sample_labels
    @test gamma.method == "gamma_plus"
    @test theta.method == "theta_plus"
    @test gamma.metric === nothing
    @test gamma == gamma'
end

@testset "Taxonomic dissimilarity follows its definition" begin
    d = community_d()
    comm = community_frame(
        "s1" => ["A" => 1, "C" => 1],
        "s2" => ["B" => 1],
        "s3" => ["D" => 1, "E" => 1],
    )
    gamma = taxo_dissimilarity(comm, d)
    theta = taxo_dissimilarity(comm, d; index="theta_plus")
    @test gamma["s1", "s2"] ≈ (d["A", "B"] + d["C", "B"] + d["B", "A"]) / 3
    @test theta["s1", "s2"] ≈ ((d["A", "B"] + d["C", "B"]) / 2 + d["B", "A"]) / 2
    # Equal sample sizes give Gamma+ = Theta+.
    @test gamma["s1", "s3"] ≈ theta["s1", "s3"]
    @test gamma.metric == "path"
    @test gamma["s1", "s1"] == 0
    shown = sprint(show, MIME"text/plain"(), gamma)
    @test occursin("method: gamma_plus", shown) && occursin("metric: path", shown)
    # Identical samples have zero dissimilarity.
    same = community_frame("a" => ["A" => 1, "B" => 1], "b" => ["A" => 1, "B" => 1])
    @test dissimilarity_values(taxo_dissimilarity(same, d)) == [0.0]
end

@testset "Taxonomic dissimilarity handles empty samples and missing distances" begin
    d = community_d()
    comm = DataFrame(sample=["s1", "s2", "s3"], A=[1, 0, 0], F=[0, 0, 1])
    values = dissimilarity_values(taxo_dissimilarity(comm, d))
    # Pairs (s2, s1) and (s3, s2) involve the empty sample s2.
    @test isnan(values[1]) && isnan(values[3])
    @test values[2] == d["A", "F"]
    disconnected = distance_matrix(
        taxo_from_lineages(["X" => ["Biota", "X"], "Y" => ["Other", "Y"]]); progress=false,
    )
    comm2 = DataFrame(sample=["a", "b"], X=[1, 0], Y=[0, 1])
    @test isnan(dissimilarity_values(taxo_dissimilarity(comm2, disconnected))[1])
    single = taxo_dissimilarity(DataFrame(sample=["a"], A=[1]), d)
    @test size(single) == (1, 1)
    @test_throws r"gamma_plus" taxo_dissimilarity(comm, d; index="jaccard")
    @test_throws r"gamma_plus" taxo_dissimilarity(comm, d; index=:gamma_plus)
    @test_throws r"no taxa" taxo_dissimilarity(String[], d)
    @test_throws r"TaxonomicDistanceMatrix" taxo_dissimilarity(comm, [0 1; 1 0])
end

# Stores lineages in the cache so that distance_matrix(names) works offline.
function cache_community_lineages()
    for (index, (name, lineage)) in enumerate(COMMUNITY_LINEAGES)
        Taxodist._taxodist_cache["id_$(name)"] = string(100 + index)
        Taxodist._taxodist_cache["lin_$(100 + index)"] = lineage
    end
end

@testset "Taxonomic dissimilarity computes distances when d is nothing" begin
    clear_cache()
    try
        cache_community_lineages()
        comm = DataFrame(sample=["a", "b"], A=[1, 0], B=[0, 1])
        out = taxo_dissimilarity(comm; metric="path", progress=false)
        @test dissimilarity_values(out) == [2.0]
        @test out.metric == "path"
    finally
        clear_cache()
    end
end

# Reference implementations written from picante's definitions.
function reference_mpd(x, m, weighted)
    p = x .> 0
    w = m[p, p]
    weighted || return mean(w[i, j] for j in axes(w, 2) for i in axes(w, 1) if i > j)
    weights = x[p] * x[p]'
    return sum(w .* weights) / sum(weights)
end
function reference_mntd(x, m, weighted)
    p = x .> 0
    w = m[p, p]
    nearest = [minimum(w[i, j] for i in axes(w, 1) if i != j) for j in axes(w, 2)]
    weighted || return mean(nearest)
    return sum(nearest .* x[p]) / sum(x[p])
end

@testset "MPD and MNTD follow picante's definitions" begin
    d = community_d()
    m = d.values
    comm = [3 1 0 2 0 1; 0 5 2 0 1 0]
    stat = Taxodist._community_statistic
    for k in 1:2, weighted in (false, true)
        x = Float64.(comm[k, :])
        @test stat(x, m, "mpd", weighted) ≈ reference_mpd(x, m, weighted)
        @test stat(x, m, "mntd", weighted) ≈ reference_mntd(x, m, weighted)
    end
    @test isnan(stat([1.0, 0.0], m[1:2, 1:2], "mpd", false))
    inf = m[1:2, 1:2]
    inf[1, 2] = inf[2, 1] = Inf
    @test isnan(stat([1.0, 1.0], inf, "mpd", false))
end

@testset "The pool null model has the pool mean as expected MPD" begin
    d = community_d()
    res = taxo_ses(["A", "B", "C"], d; n_iter=4000, rng=Xoshiro(3))
    m = d.values
    @test isapprox(res.null_mean[1], mean(Taxodist._lower_triangle(m)); rtol=0.02)
    @test res.observed[1] ≈ reference_mpd([1.0, 1.0, 1.0], m[1:3, 1:3], false)
    @test res.NRI[1] ≈ -(res.observed[1] - res.null_mean[1]) / res.null_sd[1]
    # Three closely related mammals are clustered relative to the pool.
    @test res.NRI[1] > 0
    @test res isa TaxodistSES
    @test res.pool_size == 6
    @test res.metric == "path"
    @test names(res.data) == ["sample", "S", "observed", "null_mean", "null_sd", "z",
                              "NRI", "p_lower", "p_upper", "n_iter"]
    @test res.n_iter == [4000]
    shown = sprint(show, MIME"text/plain"(), res)
    @test occursin("Net relatedness index (MPD)", shown)
    @test occursin("Null model: pool, pool of 6 taxa; metric: path", shown)
    @test occursin("sample", shown)
    # The same generator gives the same result.
    again = taxo_ses(["A", "B", "C"], d; n_iter=4000, rng=Xoshiro(3))
    @test again.null_mean == res.null_mean
end

@testset "taxo_ses supports MNTD, abundances, and custom pools" begin
    d = community_d()
    comm = community_frame(
        "s1" => ["A" => 3, "B" => 1],
        "s2" => ["C" => 2, "D" => 5],
        "s3" => ["A" => 1],
    )
    res = taxo_ses(
        comm, d; index="mntd", abundance_weighted=true,
        pool=["A", "B", "C", "D", "E"], n_iter=99, rng=Xoshiro(5),
    )
    @test names(res.data) == ["sample", "S", "observed", "null_mean", "null_sd", "z",
                              "NTI", "p_lower", "p_upper", "n_iter"]
    @test res.S == [2, 2, 1]
    @test isnan(res.observed[3])
    @test isnan(res.p_lower[3]) && isnan(res.p_upper[3])
    @test all(0 .< res.p_lower[1:2] .<= 1)
    @test all(0 .< res.p_upper[1:2] .<= 1)
    @test res.pool_size == 5
    shown = sprint(show, MIME"text/plain"(), res)
    @test occursin("abundance-weighted", shown)
    @test occursin("Nearest taxon index (MNTD)", shown)
    # A pool with an infinite distance gives missing null distributions.
    disconnected = distance_matrix(
        taxo_from_lineages(["X" => ["Biota", "X"], "Y" => ["Biota", "Y"], "Z" => ["Other", "Z"]]);
        progress=false,
    )
    partial = taxo_ses(["X", "Y"], disconnected; n_iter=50, rng=Xoshiro(1))
    @test isfinite(partial.observed[1])
    @test isnan(partial.null_mean[1]) && isnan(partial.p_lower[1])
end

@testset "The independent swap keeps richness and frequencies" begin
    rng = Xoshiro(9)
    comm = (rand(rng, 8, 10) .< 0.4) .* rand(rng, 1:5, 8, 10)
    random = Taxodist._independent_swap(comm, 200, rng)
    @test vec(sum(random .> 0; dims=2)) == vec(sum(comm .> 0; dims=2))
    @test vec(sum(random .> 0; dims=1)) == vec(sum(comm .> 0; dims=1))
    @test (random .> 0) != (comm .> 0)
    # Abundances move with presences.
    @test sort(random[random .> 0]) == sort(Float64.(comm[comm .> 0]))
    # A matrix without checkerboards cannot change.
    full = ones(2, 2)
    @test Taxodist._independent_swap(full, 5) == full
    for _ in 1:50
        first_index, second_index = Taxodist._two_distinct(rng, 3)
        @test first_index != second_index && first_index in 1:3 && second_index in 1:3
    end
end

@testset "taxo_ses runs the independent swap null model" begin
    d = community_d()
    comm = DataFrame(
        sample=["s1", "s2", "s3"],
        A=[1, 0, 1], B=[1, 1, 0], C=[0, 1, 0], D=[0, 1, 1], E=[0, 0, 1], F=[1, 0, 0],
    )
    res = taxo_ses(comm, d; null_model="independent_swap", n_iter=49, swaps=50, rng=Xoshiro(2))
    @test nrow(res.data) == 3
    @test all(isfinite, res.null_mean)
    @test res.null_model == "independent_swap"
    @test_throws r"at least two samples" taxo_ses(["A", "B"], d; null_model="independent_swap")
    # Taxa outside the pool are allowed by the swap model.
    @test taxo_ses(comm, d; null_model="independent_swap", pool=["A", "B"], n_iter=3,
                   swaps=5, rng=Xoshiro(1)).pool_size == 2
end

@testset "taxo_ses validates its arguments" begin
    d = community_d()
    @test_throws r"mpd" taxo_ses(["A", "B"], d; index="pd")
    @test_throws r"pool" taxo_ses(["A", "B"], d; null_model="frequency")
    @test_throws r"positive integer" taxo_ses(["A", "B"], d; n_iter=0)
    @test_throws r"positive integer" taxo_ses(["A", "B"], d; swaps=1.5)
    @test_throws r"positive integer" taxo_ses(["A", "B"], d; n_iter=true)
    @test_throws r"true or false" taxo_ses(["A", "B"], d; abundance_weighted=missing)
    @test_throws r"at least two" taxo_ses(["A", "B"], d; pool=["A"])
    @test_throws r"at least two" taxo_ses(["A", "B"], d; pool="AB")
    @test_throws r"belong to the pool" taxo_ses(["A", "B"], d; pool=["C", "D"])
    @test_throws r"no taxa" taxo_ses(String[], d)
    zero_sd = taxo_ses(["A", "B"], d; pool=["A", "B"], n_iter=5, rng=Xoshiro(1))
    @test isnan(zero_sd.z[1]) && isnan(zero_sd.NRI[1])
    one_iteration = taxo_ses(["A", "B"], d; n_iter=1, rng=Xoshiro(1))
    @test isnan(one_iteration.null_sd[1])
    # A user-built matrix has no metric.
    plain = TaxonomicDistanceMatrix(d.values, d.taxa)
    res = taxo_ses(["A", "B", "C"], plain; n_iter=5, rng=Xoshiro(1))
    @test res.metric === nothing
    @test !occursin("metric", sprint(show, MIME"text/plain"(), res))
end

@testset "taxo_ses computes distances when d is nothing" begin
    clear_cache()
    try
        cache_community_lineages()
        res = taxo_ses(Dict("A" => 1, "B" => 1, "D" => 1); n_iter=9, progress=false, rng=Xoshiro(1))
        @test res.pool_size == 3
        @test res.metric == "reciprocal"
    finally
        clear_cache()
    end
end

@testset "Community indices reproduce the shared fixture" begin
    fixture = JSON3.read(read(joinpath(@__DIR__, "fixtures", "expected_community.json"), String))
    lineages = [String(name) => String.(lineage) for (name, lineage) in pairs(fixture[:lineages])]
    resolution = taxo_from_lineages(lineages)
    taxa = [String(name) for name in keys(fixture[:comm][1])]
    comm = [Float64(row[Symbol(taxon)]) for row in fixture[:comm], taxon in taxa]
    sample_names = String.(fixture[:sample_names])
    expected_values(values) = [value === nothing ? NaN : Float64(value) for value in values]
    same(got, expected) = length(got) == length(expected) &&
        all(isnan(g) ? isnan(e) : isapprox(g, e; rtol=1e-12) for (g, e) in zip(got, expected))
    for metric in ("reciprocal", "path")
        d = distance_matrix(resolution; metric=metric, progress=false)
        full = Taxodist._distinctness_matrix(d, taxa)
        exp = fixture[Symbol(metric)]
        for index in ("gamma_plus", "theta_plus")
            got = taxo_dissimilarity(comm, d; index=index, taxa=taxa, samples=sample_names)
            @test same(dissimilarity_values(got), expected_values(exp[Symbol(index)]))
        end
        for index in ("mpd", "mntd"), weighted in (false, true)
            key = Symbol(index * (weighted ? "_weighted" : ""))
            got = [Taxodist._community_statistic(comm[k, :], full, index, weighted) for k in axes(comm, 1)]
            @test same(got, expected_values(exp[key]))
        end
    end
end
