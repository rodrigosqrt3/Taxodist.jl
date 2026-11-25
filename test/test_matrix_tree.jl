using Taxodist.DataFrames: nrow

const DINOSAURS = Dict(
    "Tyrannosaurus" => ["Biota", "Dinosauria", "Theropoda", "Tyrannosaurus"],
    "Velociraptor" => ["Biota", "Dinosauria", "Theropoda", "Velociraptor"],
    "Triceratops" => ["Biota", "Dinosauria", "Ornithischia", "Triceratops"],
    "Homo" => ["Biota", "Mammalia", "Homo"],
)

function random_lineages(n; depth=8, branching=3)
    lineages = Dict{String,Vector{String}}()
    state = UInt64(12345)
    next_int(k) = (state = state * 6364136223846793005 + 1442695040888963407; Int((state >> 33) % k) + 1)
    for i in 1:n
        len = next_int(depth - 1) + 1
        path = [next_int(branching) for _ in 1:len]
        lineages["t$(i)"] = vcat(["Biota"], ["N" * join(path[1:d]) for d in 1:len])
    end
    return lineages
end

@testset "Prefix-identifier distances match pairwise distances" begin
    lineages = random_lineages(40)
    lineages["dup"] = copy(lineages["t1"])
    lineages["ancestor"] = lineages["t2"][1:2]
    lineages["foreign"] = ["Natura", "Mineralia"]
    resolution = taxo_from_lineages(lineages)
    matrix = distance_matrix(resolution; progress=false)
    for i in eachindex(matrix.taxa), j in eachindex(matrix.taxa)
        expected = Taxodist._compute_distance(
            lineages[matrix.taxa[i]], lineages[matrix.taxa[j]],
        ).distance
        @test matrix[i, j] == expected
    end
    ids, missing_mask = Taxodist._prefix_node_ids([["Biota"], nothing])
    @test missing_mask == [false, true]
    @test isempty(ids[2])
    @test Taxodist._metric_value(0, 0, 0, :reciprocal) == Inf
end

# Brute-force MRCA depths and distances, as in R's tests of .mrca_index().
function brute_distance(a, b, metric)
    (a === nothing || b === nothing) && return NaN
    return Taxodist._compute_distance(a, b; metric=metric).distance
end

@testset "Sorted-prefix MRCA index agrees with pairwise comparison" begin
    state = Ref(UInt64(20261124))
    next_int(k) = (state[] = state[] * 6364136223846793005 + 1442695040888963407;
                   Int((state[] >> 33) % k) + 1)
    random_lineage() = ["n$(d)_$(next_int(2) - 1)" for d in 1:next_int(12)]
    lineages = Union{Nothing,Vector{String}}[random_lineage() for _ in 1:40]
    append!(lineages, Union{Nothing,Vector{String}}[
        nothing,                                  # unresolved taxon
        String[],                                 # empty lineage
        ["n1_0", "n2_0"],                         # ancestor of others
        ["n1_0", "n2_0"],                         # duplicate
        ["Other", "root"],                        # disconnected hierarchy
        ["n1_0", "n2_0", "n3_1", "n4_0", "n5_1"],
    ])
    n = length(lineages)
    labels = ["t$(i)" for i in 1:n]
    for metric in ("reciprocal", "path", "wu_palmer")
        fast = Taxodist._distance_matrix_from_lineages(labels, lineages; metric=metric)
        for j in 1:n, i in 1:n
            expected = i == j ? 0.0 : brute_distance(lineages[i], lineages[j], metric)
            @test isequal(fast[i, j], expected)
        end
    end

    ids, _ = Taxodist._prefix_node_ids(lineages)
    index = Taxodist._mrca_index(ids)
    for j in (1, 20, n - 1, n)
        expected = [
            l === nothing || lineages[j] === nothing ? 0 :
            Taxodist._common_prefix_depth(l, lineages[j]) for l in lineages
        ]
        @test Taxodist._mrca_depths_from(index, j) == expected
    end

    # The ranked metric goes through the same index on reduced lineages.
    ranked = Dict(
        "a" => ["Biota", "Animalia", "Felidae", "A"],
        "b" => ["Biota", "Animalia", "Canidae", "B"],
        "c" => ["Biota", "Plantae", "C"],
    )
    rank_vectors = Dict(
        "a" => [missing, "Kingdom", "Family", "Species"],
        "b" => [missing, "Kingdom", "Family", "Species"],
        "c" => [missing, "Kingdom", "Species"],
    )
    resolution = taxo_from_lineages(ranked; ranks=rank_vectors)
    matrix = distance_matrix(resolution; progress=false, metric="ranked")
    for x in matrix.taxa, y in matrix.taxa
        expected = Taxodist._compute_distance(
            (nodes=ranked[x], ranks=rank_vectors[x]),
            (nodes=ranked[y], ranks=rank_vectors[y]); metric="ranked",
        ).distance
        @test matrix[x, y] == expected
    end
end

@testset "cross_distances uses the same index as distance_matrix" begin
    state = Ref(UInt64(7))
    next_int(k) = (state[] = state[] * 6364136223846793005 + 1442695040888963407;
                   Int((state[] >> 33) % k) + 1)
    lineages = [
        "t$(i)" => ["c$(d)$(next_int(2) - 1)" for d in 1:(next_int(7) + 1)]
        for i in 1:25
    ]
    x = taxo_from_lineages(lineages[1:10])
    y = taxo_from_lineages(lineages[11:25])
    for metric in ("reciprocal", "path", "wu_palmer")
        full = distance_matrix(taxo_from_lineages(lineages); progress=false, metric=metric)
        cross = cross_distances(x, y; progress=false, metric=metric)
        @test isequal(cross.values, full.values[1:10, 11:25])
    end

    # Unresolved taxa on either side give NaN; empty lineages behave as in
    # distance_matrix.
    with_missing = Union{Nothing,Vector{String}}[["a", "b"], nothing, String[]]
    labels = ["p", "q", "r"]
    ids, missing_mask = Taxodist._prefix_node_ids(with_missing)
    @test missing_mask == [false, true, false]
    square = Taxodist._distance_matrix_from_lineages(labels, with_missing; metric="wu_palmer")
    @test isnan(square[1, 2]) && isnan(square[2, 3])
    @test square[1, 3] == 1.0 && square[3, 1] == 1.0   # no shared root
    @test square[3, 3] == 0.0
    # Two empty lineages: 1 - 0/0, as in R.
    @test isnan(Taxodist._distance_matrix_from_lineages(
        ["e1", "e2"], [String[], String[]]; metric="wu_palmer",
    )[1, 2])
    cross = cross_distances(
        taxo_from_lineages(["a" => ["r", "a"]]), [missing, "  "];
        progress=false,
    )
    @test all(isnan, cross.values)
    @test all(isnan, cross_distances([missing], taxo_from_lineages(["a" => ["r", "a"]]);
                                     progress=false).values)
end

@testset "MRCA index handles degenerate inputs" begin
    empty_index = Taxodist._mrca_index(Vector{Int}[])
    @test isempty(empty_index.order) && isempty(empty_index.adj)
    single = Taxodist._mrca_index(first(Taxodist._prefix_node_ids([["a", "b"]])))
    @test Taxodist._mrca_depths_from(single, 1) == [2]
    flat = Taxodist._mrca_index(first(Taxodist._prefix_node_ids([String[], String[]])))
    @test Taxodist._mrca_depths_from(flat, 1) == [0, 0]
    @test size(Taxodist._distance_matrix_from_lineages(String[], Vector{String}[])) == (0, 0)
end

@testset "Cross distances" begin
    lineages = random_lineages(12)
    names_ = sort(collect(keys(lineages)))
    x = taxo_from_lineages(Dict(k => lineages[k] for k in names_[1:5]))
    y = taxo_from_lineages(Dict(k => lineages[k] for k in names_[6:end]))
    full = distance_matrix(taxo_from_lineages(lineages); progress=false)
    cross = cross_distances(x, y; progress=false)
    @test size(cross) == (5, 7)
    @test Base.IndexStyle(TaxonomicCrossDistances) == Base.IndexCartesian()
    for (i, row) in enumerate(cross.rows), (j, column) in enumerate(cross.columns)
        @test cross[i, j] == full[row, column]
        @test cross[row, column] == full[row, column]
    end
    @test Matrix(cross) == cross.values
    @test_throws KeyError cross["Missing", cross.columns[1]]
    @test_throws KeyError cross[cross.rows[1], "Missing"]
    @test occursin("TaxonomicCrossDistances", sprint(show, MIME"text/plain"(), cross))

    bundle = taxo_bundle(taxo_from_lineages(DINOSAURS); progress=false)
    @test cross_distances(bundle, bundle; progress=false)["Homo", "Tyrannosaurus"] == 1.0

    clear_cache()
    Taxodist._taxodist_cache["id_Alpha"] = "1"
    Taxodist._taxodist_cache["lin_1"] = ["Biota", "Animalia", "Alpha"]
    Taxodist._taxodist_cache["id_Gone"] = nothing
    by_name = cross_distances(["Alpha", "Gone"], ["Alpha"]; progress=false)
    @test by_name["Alpha", "Alpha"] == 0.0
    @test isnan(by_name["Gone", "Alpha"])
    @test size(cross_distances(String[], ["Alpha"]; progress=false)) == (0, 1)
    clear_cache()
end

@testset "Classification trees" begin
    tree = taxo_tree(taxo_from_lineages(DINOSAURS); progress=false)
    @test tree isa TaxodistTree
    @test Taxodist.n_tips(tree) == 4
    @test tree.node_label[1] == "Biota"
    @test Set(tree.node_label) == Set(["Biota", "Dinosauria", "Theropoda"])
    @test tree.edge[1, 1] == 5
    @test all(>=(0), tree.edge_length)
    @test occursin("TaxodistTree with 4 tips", sprint(show, MIME"text/plain"(), tree))

    lineages = merge(DINOSAURS, Dict(
        "Dinosauria" => ["Biota", "Dinosauria"],
        "Rex1" => DINOSAURS["Tyrannosaurus"],
        "Rex2" => DINOSAURS["Tyrannosaurus"],
        "Theropoda" => ["Biota", "Dinosauria", "Theropoda"],
        "Theropoda2" => ["Biota", "Dinosauria", "Theropoda"],
    ))
    resolution = taxo_from_lineages(lineages)
    full_tree = taxo_tree(resolution; progress=false)
    expected = distance_matrix(resolution; progress=false)
    observed = cophenetic(full_tree)
    for a in expected.taxa, b in expected.taxa
        @test observed[a, b] ≈ expected[a, b] atol = 1e-12
    end
    @test "" in full_tree.node_label
    text = newick(full_tree)
    @test endswith(text, ")Biota;")
    @test occursin(")Tyrannosaurus:", text)

    quoted = taxo_tree(taxo_from_lineages(Dict(
        "Homo sapiens" => ["Biota", "Homo", "Homo sapiens"],
        "O'Brien" => ["Biota", "Homo", "O'Brien"],
        "Pan" => ["Biota", "Pan"],
    )); progress=false)
    @test occursin("'Homo sapiens'", newick(quoted))
    @test occursin("'O''Brien'", newick(quoted))

    clear_cache()
    for (index, (name, lineage)) in enumerate(DINOSAURS)
        Taxodist._taxodist_cache["id_$(name)"] = string(index)
        Taxodist._taxodist_cache["lin_$(index)"] = lineage
    end
    Taxodist._taxodist_cache["id_Missing"] = nothing
    named = @test_logs (:warn, r"dropped") taxo_tree(
        vcat(collect(keys(DINOSAURS)), ["Missing"]); progress=false,
    )
    @test Taxodist.n_nodes(named) == 3
    clear_cache()

    bundle = taxo_bundle(taxo_from_lineages(DINOSAURS); progress=false)
    @test Taxodist.n_nodes(taxo_tree(bundle)) == 3

    @test_throws ArgumentError taxo_tree(taxo_from_lineages(Dict("A" => ["Biota", "A"])))
    @test_throws ArgumentError taxo_tree(taxo_from_lineages(Dict(
        "A" => ["Biota", "A"], "B" => ["Natura", "B"],
    )))
end

@testset "Heatmap keyword arguments" begin
    matrix = distance_matrix(taxo_from_lineages(DINOSAURS); progress=false)
    @test taxo_heatmap(matrix; display_plot=false, color=:inferno, title="Custom") === matrix
    resolution = taxo_from_lineages(DINOSAURS)
    @test size(taxo_heatmap(resolution; display_plot=false, progress=false)) == (4, 4)
end

@testset "Bundle comparison" begin
    make_bundle(lineages; ids=nothing) = taxo_bundle(
        taxo_from_lineages(lineages; ids=ids); progress=false,
    )
    lineages = Dict("Alpha" => ["Biota", "A", "Alpha"], "Beta" => ["Biota", "B", "Beta"])
    same = compare_bundles(make_bundle(lineages), make_bundle(lineages))
    @test all(==("unchanged"), same.taxa.change)
    @test nrow(same.distances) == 0
    @test occursin("No taxon changed", sprint(show, MIME"text/plain"(), same))

    old = make_bundle(Dict(
        "Same" => ["Biota", "Animalia", "Same"],
        "Moved" => ["Biota", "Animalia", "Moved"],
        "Renumbered" => ["Biota", "Animalia", "Renumbered"],
        "Dropped" => ["Biota", "Plantae", "Dropped"],
    ); ids=Dict("Same" => "1", "Moved" => "2", "Renumbered" => "3", "Dropped" => "4"))
    new = make_bundle(Dict(
        "Same" => ["Biota", "Animalia", "Same"],
        "Moved" => ["Biota", "Plantae", "Moved"],
        "Renumbered" => ["Biota", "Animalia", "Renumbered"],
        "Added" => ["Biota", "Fungi", "Added"],
    ); ids=Dict("Same" => "1", "Moved" => "2", "Renumbered" => "33", "Added" => "5"))
    comparison = compare_bundles(old, new)
    changes = Dict(zip(comparison.taxa.input, comparison.taxa.change))
    @test changes == Dict(
        "Same" => "unchanged", "Moved" => "lineage_changed",
        "Renumbered" => "id_changed", "Dropped" => "removed", "Added" => "added",
    )
    pairs_ = Set(zip(comparison.distances.taxon_a, comparison.distances.taxon_b))
    @test length(pairs_) == 2
    @test Set(Set.(collect.(pairs_))) == Set([Set(["Same", "Moved"]), Set(["Moved", "Renumbered"])])
    shown = sprint(show, MIME"text/plain"(), comparison)
    @test occursin("Moved: lineage_changed", shown)

    mktempdir() do directory
        path = joinpath(directory, "old.json")
        write_taxodist_bundle(old, path)
        @test compare_bundles(path, old).taxa.change == fill("unchanged", 4)
    end

    resolved = make_bundle(Dict("Alpha" => ["Biota", "Alpha"], "Beta" => ["Biota", "Beta"]))
    data = copy(resolved.resolution.data)
    row = findfirst(==("Beta"), data.input)
    data[row, :status] = "retrieval_error"
    data[row, :id] = nothing
    data[row, :resolved_name] = nothing
    data[row, :lineage] = nothing
    data[row, :lineage_depth] = nothing
    data[row, :candidates] = Taxodist._empty_candidates()
    data[row, :n_candidates] = 0
    failed = taxo_bundle(
        Taxodist.TaxodistResolution(data, "The Taxonomicon", nothing, nothing);
        progress=false,
    )
    status = compare_bundles(resolved, failed)
    @test status.taxa.change[row] == "status_changed"
    @test nrow(status.distances) == 1
    @test isnan(status.distances.distance_new[1])
    @test nrow(compare_bundles(failed, failed).distances) == 0

    @test_throws ArgumentError compare_bundles(42, resolved)

    disconnected = Dict("A" => ["Root A", "A"], "B" => ["Root B", "B"])
    @test nrow(compare_bundles(make_bundle(disconnected), make_bundle(disconnected)).distances) == 0

    # Metrics of the two bundles.
    @test same.metric_old == "reciprocal" && same.metric_new == "reciprocal"
    @test !same.metric_changed
    @test !occursin("different", sprint(show, MIME"text/plain"(), same))
    resolution = taxo_from_lineages(lineages)
    reciprocal = taxo_bundle(resolution; progress=false)
    path = taxo_bundle(resolution; progress=false, metric="path")
    across = compare_bundles(reciprocal, path)
    @test across.metric_old == "reciprocal"
    @test across.metric_new == "path"
    @test across.metric_changed
    shown = sprint(show, MIME"text/plain"(), across)
    @test occursin("! The bundles use different metrics: \"reciprocal\" and \"path\".", shown)
    @test occursin("No taxon changed", shown)
    @test occursin("Distances changed: 1", shown)

    ranks = Dict("Alpha" => [missing, "Phylum", "Genus"], "Beta" => [missing, "Phylum", "Genus"])
    ranked_resolution = taxo_from_lineages(lineages; ranks=ranks)
    ranked_default = taxo_bundle(ranked_resolution; progress=false, metric="ranked")
    levels = copy(Taxodist.RANK_LEVELS[])
    ranked_other = try
        Taxodist.RANK_LEVELS[] = ["phylum", "genus"]
        taxo_bundle(ranked_resolution; progress=false, metric="ranked")
    finally
        Taxodist.RANK_LEVELS[] = levels
    end
    levels_changed = compare_bundles(ranked_default, ranked_other)
    @test levels_changed.metric_old == levels_changed.metric_new == "ranked"
    @test levels_changed.metric_changed
    @test occursin("! The bundles use different rank levels for the \"ranked\" metric.",
                   sprint(show, MIME"text/plain"(), levels_changed))
    @test !compare_bundles(ranked_default, ranked_default).metric_changed

    # Using a bundle skips recomputation; writing still checks distances.
    tampered = make_bundle(Dict("Alpha" => ["Biota", "Alpha"], "Beta" => ["Biota", "Beta"]))
    tampered.matrix.values[1, 2] = tampered.matrix.values[2, 1] = 0.75
    @test distance_matrix(tampered)[1, 2] == 0.75
    @test occursin("TaxodistBundle", sprint(show, MIME"text/plain"(), tampered))
    mktempdir() do directory
        @test_throws ArgumentError write_taxodist_bundle(tampered, joinpath(directory, "x.json"))
    end
end
