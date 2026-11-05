# Ranks, alternative distance metrics, schema 1.1 bundles, and taxonomic
# distinctness. Reference values marked "from R" were computed with the R
# implementation (taxodist 1.0.0) and agree with vegan::taxondive().

using Taxodist: mean, quantile
using Taxodist.Random: Xoshiro
using Taxodist.DataFrames: DataFrame, select, Not

const METRIC_FIXTURES = joinpath(@__DIR__, "fixtures")
metric_fixture(name) = joinpath(METRIC_FIXTURES, name)

const RANKED_LINEAGES = [
    "A" => ["Biota", "Animalia", "Chordata", "Mammalia", "Carnivora", "Felidae", "Panthera", "A"],
    "B" => ["Biota", "Animalia", "Chordata", "Mammalia", "Carnivora", "Felidae", "Panthera", "B"],
    "C" => ["Biota", "Animalia", "Chordata", "Mammalia", "Carnivora", "Canidae", "Canis", "C"],
    "D" => ["Biota", "Animalia", "Chordata", "Mammalia", "Primates", "Hominidae", "Homo", "D"],
    "E" => ["Biota", "Animalia", "Arthropoda", "Insecta", "Diptera", "Muscidae", "Musca", "E"],
]
const STANDARD_RANKS = [missing, "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

ranked_example() = taxo_from_lineages(
    RANKED_LINEAGES; ranks=[first(entry) => STANDARD_RANKS for entry in RANKED_LINEAGES],
)

subset_resolution(resolution, rows) = TaxodistResolution(
    resolution.data[rows, :], resolution.source, resolution.source_url, resolution.retrieved_at,
)

# Matrices of the 5-taxon example computed by R's distance_matrix().
const R_METRIC_MATRICES = Dict(
    "reciprocal" => [
        0 1/7 0.2 0.25 0.5
        1/7 0 0.2 0.25 0.5
        0.2 0.2 0 0.25 0.5
        0.25 0.25 0.25 0 0.5
        0.5 0.5 0.5 0.5 0
    ],
    "path" => [
        0 2 6 8 12
        2 0 6 8 12
        6 6 0 8 12
        8 8 8 0 12
        12 12 12 12 0
    ],
    "wu_palmer" => [
        0 0.125 0.375 0.5 0.75
        0.125 0 0.375 0.5 0.75
        0.375 0.375 0 0.5 0.75
        0.5 0.5 0.5 0 0.75
        0.75 0.75 0.75 0.75 0
    ],
    "ranked" => [
        0 1 3 4 6
        1 0 3 4 6
        3 3 0 4 6
        4 4 4 0 6
        6 6 6 6 0
    ],
)

# -- Ranks --------------------------------------------------------------------

@testset "The structured parser reads the rank of every fixture node" begin
    expected = Taxodist.JSON3.read(read(metric_fixture("expected_ranks.json"), String))
    @test length(expected) > 5
    for (id, ranks) in pairs(expected)
        want = Union{Missing,String}[rank === nothing ? missing : String(rank) for rank in ranks]
        parsed = Taxodist._lineage_and_ranks_from_page(
            read(metric_fixture("tree_$(id).html"), String), String(id),
        )
        @test isequal(parsed.ranks, want)
        @test parsed.nodes == Taxodist._lineage_from_page(
            read(metric_fixture("tree_$(id).html"), String), String(id),
        )
    end

    passer = Taxodist._lineage_and_ranks_from_page(
        read(metric_fixture("tree_59107.html"), String), "59107",
    )
    rank_of = Dict(zip(passer.nodes, passer.ranks))
    @test [rank_of[node] for node in [
        "Animalia", "Chordata", "Reptilia", "Aves", "Passeriformes",
        "Passeridae", "Passer", "Neornithes", "Biota",
    ]] == [
        "Kingdom", "Phylum", "Class", "Subclass", "Order", "Family", "Genus",
        "Clade", "Clade",
    ]
end

@testset "Ranks are inferred from the form of unlabelled names" begin
    @test Taxodist._rank_from_name("Homo sapiens") == "Species"
    @test Taxodist._rank_from_name("Homo sapiens sapiens") == "Subspecies"
    @test Taxodist._rank_from_name("Microtus (Agricola)") == "Subgenus"
    @test ismissing(Taxodist._rank_from_name("Natura"))
    @test Taxodist._normalize_rank("[plesion] Family †", "X") == "Family"
    @test Taxodist._normalize_rank("Taxon: Genus †", "Yi") == "Genus"
    @test Taxodist._normalize_rank("Clade \"", "X") == "Clade"
    @test Taxodist._rank_before_name("Superfamily ‡?Lipotoidea de Muizon", "Lipotoidea") ==
          "Superfamily"
    @test Taxodist._rank_before_name("Homo sapiens Linnaeus", "Homo sapiens") == "Species"

    # The page subject of a page without its own row gets its rank too.
    missing_self = """
    <div id="ctl00_divSubject">Taxon: Genus <b>†Delta?</b></div>
    <div id="ctl00_cphTabContent_divHierarchy">
    <a>Clade <b>Alpha</b></a><br>  <a>Family <b>Beta</b></a>
    <br>  <a><b>Sibling</b></a></div>
    """
    parsed = Taxodist._lineage_and_ranks_from_page(missing_self, "1")
    @test parsed.nodes == ["Alpha", "Beta", "Delta"]
    @test parsed.ranks == ["Clade", "Family", "Genus"]
end

@testset "The text fallback parser also records ranks" begin
    fake_html = "<html><body><div id=\"ctl00_divSubject\"><b>Homo sapiens</b></div>" *
                "<div id=\"divPageContent\">Natura - nature\nClade Biota Wagner 2004\n" *
                "Kingdom Animalia Linnaeus, 1758\n[crown] Clade Theria\n" *
                "Genus Homo Linnaeus, 1758\nHomo sapiens Linnaeus, 1758\n</div></body></html>"
    parsed = Taxodist._lineage_and_ranks_from_page(fake_html, "1")
    @test parsed.nodes == ["Natura", "Biota", "Animalia", "Theria", "Homo", "Homo sapiens"]
    @test isequal(parsed.ranks, [missing, "Clade", "Kingdom", "Clade", "Genus", "Species"])
    @test Taxodist._label_rank("Anamorph genus Epidermophyton Sabour.") == "Anamorph genus"
    @test Taxodist._label_rank("‡?[plesion] Family †Lipotidae") == "Family"
end

@testset "get_lineage_by_id caches the ranks of both variants" begin
    original = Taxodist._http_get[]
    clear_cache()
    html = read(metric_fixture("tree_59107.html"), String)
    try
        Taxodist._http_get[] = (url, headers; kwargs...) ->
            (status=200, body=Vector{UInt8}(codeunits(html)))
        cleaned = get_lineage_by_id("59107")
        raw = get_lineage_by_id("59107"; clean=false)
        cache = Taxodist._taxodist_cache
        @test length(cache["ranks_59107"]) == length(cleaned)
        @test length(cache["rawranks_59107"]) == length(raw)
        @test ismissing(first(cache["rawranks_59107"]))
        @test first(cache["ranks_59107"]) == "Clade"

        by_id = lineage_ranks(59107)
        @test by_id.node == cleaned
        @test isequal(by_id.rank, cache["ranks_59107"])
        @test lineage_ranks("59107"; clean=false).node == raw

        # Lineages cached without ranks (an older cache file) are fetched again.
        delete!(cache, "ranks_59107")
        @test length(Taxodist._ranks_by_id("59107")) == length(cleaned)
        @test haskey(cache, "ranks_59107")

        # Ranks survive a cache file round trip, with null for missing.
        mktempdir() do directory
            file = joinpath(directory, "cache.json")
            save_cache(file)
            saved = cache["rawranks_59107"]
            clear_cache()
            load_cache(file)
            @test isequal(Taxodist._taxodist_cache["rawranks_59107"], saved)
            write(file, """{"meta_parser_version": 2, "ranks_1": [1]}""")
            @test_throws ArgumentError load_cache(file)
        end
    finally
        Taxodist._http_get[] = original
        clear_cache()
    end
end

@testset "lineage_ranks returns nodes aligned with get_lineage" begin
    clear_cache()
    try
        cache = Taxodist._taxodist_cache
        cache["id_Passer"] = "59107"
        cache["lin_59107"] = ["Biota", "Animalia", "Passeridae", "Passer", "Passer domesticus"]
        cache["ranks_59107"] = ["Clade", "Kingdom", "Family", "Genus", "Species"]

        out = lineage_ranks("Passer")
        @test out.node == get_lineage("Passer")
        @test out.rank == ["Clade", "Kingdom", "Family", "Genus"]
        @test out.depth == 1:4

        # A queried name absent from the page is appended with its implied rank.
        cache["id_Passer montanus"] = "59107"
        out = lineage_ranks("Passer montanus")
        @test last(out.node) == "Passer montanus"
        @test last(out.rank) == "Species"
        @test lineage_ranks(missing) === nothing
        @test lineage_ranks(nothing) === nothing
        @test lineage_ranks("  ") === nothing
        @test_throws ArgumentError lineage_ranks(["A", "B"])

        # Misaligned cached ranks are reported as missing.
        cache["ranks_59107"] = ["Clade"]
        @test all(ismissing, lineage_ranks("59107").rank)
        @test all(ismissing, Taxodist._cached_ranks("59107", cache["lin_59107"]))
    finally
        clear_cache()
    end
end

@testset "Finalizing keeps ranks aligned with the truncated lineage" begin
    lineage = ["Biota", "[unranked] x", "Microtus", "Microtus (Microtus)", "Microtus arvalis"]
    ranks = ["Clade", missing, "Genus", "Subgenus", "Species"]
    @test Taxodist._finalize_ranks(ranks, lineage, "Microtus") == ["Clade", "Genus"]
    @test Taxodist._finalize_ranks(ranks, lineage, "Microtus arvalis") ==
          ["Clade", "Genus", "Subgenus", "Species"]
    @test Taxodist._finalize_ranks(ranks, lineage, "Microtus agrestis") ==
          ["Clade", "Genus", "Subgenus", "Species"]
end

@testset "taxo_resolve records the ranks of the selected lineage" begin
    original = Taxodist._http_get[]
    original_parser = Taxodist._search_parser[]
    clear_cache()
    try
        Taxodist._http_get[] = (url, headers; kwargs...) -> (status=200, body=UInt8[])
        cache = Taxodist._taxodist_cache
        cache["lin_7"] = ["Biota", "Animalia", "Passer", "Passer domesticus"]
        cache["ranks_7"] = ["Clade", "Kingdom", "Genus", "Species"]
        Taxodist._search_parser[] = html -> [(id="7", name="Passer (Genus)")]
        resolved = taxo_resolve(["Passer", "7", "Nothing here"]; progress=false)
        @test resolved.ranks[1] == ["Clade", "Kingdom", "Genus"]
        @test resolved.ranks[2] == ["Clade", "Kingdom", "Genus", "Species"]
        Taxodist._search_parser[] = html -> NamedTuple{(:id, :name),Tuple{String,String}}[]
        unresolved = taxo_resolve(["Nothing"]; progress=false)
        @test unresolved.ranks == [nothing]
    finally
        Taxodist._http_get[] = original
        Taxodist._search_parser[] = original_parser
        clear_cache()
    end
end

@testset "taxo_from_lineages validates and stores ranks" begin
    res = ranked_example()
    @test "ranks" in names(res.data)
    @test res.ranks[1][2] == "Kingdom"
    @test ismissing(res.ranks[1][1])

    plain = taxo_from_lineages(["A" => ["Biota", "A"]])
    @test all(ismissing, plain.ranks[1]) && length(plain.ranks[1]) == 2

    positional = taxo_from_lineages(
        ["A" => ["Biota", "A"], "B" => ["Biota", "B"]];
        ranks=[[nothing, "Genus"], ["Clade", "Genus"]],
    )
    @test isequal(positional.ranks[1], [missing, "Genus"])
    keyed = taxo_from_lineages(
        Dict("A" => ["Biota", "A"], "B" => ["Biota", "B"]);
        ranks=Dict("B" => ["Clade", "Genus"], "A" => [missing, "Genus"]),
    )
    @test keyed.input == ["A", "B"]
    @test keyed.ranks[2] == ["Clade", "Genus"]

    @test_throws r"one rank" taxo_from_lineages(["A" => ["Biota", "A"]]; ranks=["A" => ["Kingdom"]])
    @test_throws r"every lineage name" taxo_from_lineages(
        ["A" => ["Biota", "A"]]; ranks=Dict("B" => [missing, "Genus"]),
    )
    @test_throws ArgumentError taxo_from_lineages(["A" => ["Biota", "A"]]; ranks="x")
    @test_throws r"one rank" taxo_from_lineages(["A" => ["Biota", "A"]]; ranks=[[1, 2]])
    @test_throws r"Positional ranks" taxo_from_lineages(
        Dict("A" => ["Biota", "A"]); ranks=[[missing, "Genus"]],
    )
end

# -- Metrics ------------------------------------------------------------------

@testset "Every metric gives the documented values and matches R" begin
    res = ranked_example()
    get(metric) = distance_matrix(res; metric=metric, progress=false)

    reciprocal = get("reciprocal")
    @test reciprocal["A", "B"] ≈ 1 / 7
    @test reciprocal["A", "E"] ≈ 1 / 2
    path = get("path")
    @test (path["A", "B"], path["A", "C"], path["A", "E"]) == (2, 6, 12)
    wu = get("wu_palmer")
    @test wu["A", "B"] ≈ 1 - 14 / 16
    @test wu["A", "E"] ≈ 1 - 4 / 16
    ranked = get("ranked")
    # Clarke & Warwick steps: genus 1, family 2, order 3, class 4, phylum 5,
    # kingdom 6 (Biota has no standard rank).
    @test (ranked["A", "B"], ranked["A", "C"], ranked["A", "D"], ranked["A", "E"]) == (1, 3, 4, 6)

    for metric in Taxodist.DISTANCE_METRICS
        matrix = get(metric)
        @test matrix.metric == metric
        @test matrix.taxa == ["A", "B", "C", "D", "E"]
        @test Matrix(matrix) == R_METRIC_MATRICES[metric]
        @test all(iszero, [matrix[i, i] for i in 1:5])
        @test Matrix(matrix) == permutedims(Matrix(matrix))
        @test get(Symbol(metric)).metric == metric
    end
    @test occursin("metric: path", sprint(show, MIME"text/plain"(), path))
    @test TaxonomicDistanceMatrix(zeros(1, 1), ["A"]).metric === nothing
end

@testset "Reciprocal is ultrametric and path satisfies the triangle inequality" begin
    res = ranked_example()
    rec = Matrix(distance_matrix(res; progress=false))
    pth = Matrix(distance_matrix(res; metric="path", progress=false))
    for i in 1:5, j in 1:5, l in 1:5
        @test rec[i, l] <= max(rec[i, j], rec[j, l]) + 1e-12
        @test pth[i, l] <= pth[i, j] + pth[j, l] + 1e-12
    end
end

@testset "Metrics treat ancestors, identical lineages and disconnected roots" begin
    ancestor = ["Biota", "Animalia", "Chordata"]
    descendant = ["Biota", "Animalia", "Chordata", "Mammalia", "Homo"]
    other = ["Root2", "X"]
    cd(a, b, metric) = Taxodist._compute_distance(a, b; metric=metric)
    @test cd(ancestor, descendant, "path").distance == 2
    @test cd(ancestor, descendant, "wu_palmer").distance ≈ 1 - 6 / 8
    @test cd(descendant, descendant, "wu_palmer").distance == 0
    @test cd(descendant, descendant, "path").distance == 0
    @test cd(ancestor, other, "path").distance == Inf
    @test cd(ancestor, other, "reciprocal").distance == Inf
    @test cd(ancestor, other, "wu_palmer").distance == 1
    @test cd(ancestor, descendant, "path").metric == "path"
    @test cd(ancestor, descendant, "reciprocal").metric == "reciprocal"
    @test Taxodist._distance_value(ancestor, other; metric="wu_palmer") == 1
    @test Taxodist._prefix_distance([1, 2], [3], :ranked) == Inf

    # Ranked distances on reduced lineages; the MRCA refers to the full ones.
    ranked_a = (nodes=["Biota", "Animalia", "Felidae", "A"],
                ranks=[missing, "Kingdom", "Family", "Species"])
    ranked_b = (nodes=["Biota", "Animalia", "Canidae", "B"],
                ranks=[missing, "Kingdom", "Family", "Species"])
    result = cd(ranked_a, ranked_b, "ranked")
    @test result.distance == 2
    @test result.mrca == "Animalia"
    @test result.depth_a == 4
    @test cd(ranked_a, ranked_a, "ranked").distance == 0
    @test_throws r"requires the rank" cd(ranked_a, ["Biota", "B"], "ranked")
    @test_throws r"Unknown distance metric" cd(ancestor, descendant, "nope")
end

@testset "The ranked metric uses RANK_LEVELS" begin
    res = ranked_example()
    original = Taxodist.RANK_LEVELS[]
    try
        Taxodist.RANK_LEVELS[] = ["Kingdom", " phylum ", "class", "order"]
        ranked = distance_matrix(res; metric="ranked", progress=false)
        # Reduced lineages: root, Kingdom, Phylum, Class, Order, terminal.
        @test (ranked["A", "B"], ranked["A", "C"], ranked["A", "D"]) == (1, 1, 2)
        @test taxo_bundle(res; metric="ranked", progress=false).metric["rank_levels"] ==
              ["kingdom", "phylum", "class", "order"]
        Taxodist.RANK_LEVELS[] = String[]
        @test_throws ArgumentError distance_matrix(res; metric="ranked", progress=false)
    finally
        Taxodist.RANK_LEVELS[] = original
    end
    @test Taxodist._rank_levels() == collect(Taxodist.DEFAULT_RANK_LEVELS)

    # The last word of the rank is compared.
    @test Taxodist._project_lineage(
        ["Root", "F", "G", "x", "T"], [missing, "Family", "Anamorph genus", "Clade", missing],
        ["family", "genus"],
    ) == ["Root", "F", "G", "T"]
end

@testset "The ranked metric keeps taxa of different kingdoms connected" begin
    res = taxo_from_lineages(
        [
            "Homo" => ["Biota", "Animalia", "Chordata", "Homo"],
            "Quercus" => ["Biota", "Plantae", "Tracheophyta", "Quercus"],
        ];
        ranks=Dict(
            "Homo" => ["Clade", "Kingdom", "Phylum", "Genus"],
            "Quercus" => ["Clade", "Kingdom", "Phylum", "Genus"],
        ),
    )
    @test distance_matrix(res; metric="ranked", progress=false)["Homo", "Quercus"] == 3
end

@testset "The ranked metric requires ranks" begin
    plain = taxo_from_lineages(["A" => ["Biota", "A"], "B" => ["Biota", "B"]])
    @test_throws r"requires the rank" distance_matrix(plain; metric="ranked", progress=false)
    @test_throws r"Unknown distance metric" distance_matrix(plain; metric="nope", progress=false)
    @test_throws ArgumentError distance_matrix(plain; metric=1, progress=false)

    # A resolution without a ranks column (e.g. a schema 1.0 bundle).
    no_column = TaxodistResolution(
        select(plain.data, Not(:ranks)), plain.source, plain.source_url, plain.retrieved_at,
    )
    @test distance_matrix(no_column; progress=false).metric == "reciprocal"
    @test_throws r"requires the rank" distance_matrix(no_column; metric="ranked", progress=false)
end

@testset "metric is honoured by the pairwise and analysis functions" begin
    res = ranked_example()
    first_two = subset_resolution(res, 1:2)
    rest = subset_resolution(res, 3:5)
    for metric in ("path", "ranked", "wu_palmer")
        cross = cross_distances(first_two, rest; metric=metric, progress=false)
        full = Matrix(distance_matrix(res; metric=metric, progress=false))
        @test Matrix(cross) == full[1:2, 3:5]
    end

    cluster = taxo_cluster(res; metric="path", progress=false)
    @test cluster.dist.metric == "path"
    @test cluster.hclust !== nothing
    ordination = taxo_ordinate(res; metric="wu_palmer", progress=false)
    @test ordination.dist.metric == "wu_palmer"
    heatmap = taxo_heatmap(res; metric="ranked", progress=false, display_plot=false)
    @test heatmap.metric == "ranked"
end

@testset "taxo_distance, closest_relative and focal_distances accept a metric" begin
    lineages = Dict(
        "Tyrannosaurus" => ["Biota", "Dinosauria", "Theropoda", "Tyrannosaurus"],
        "Velociraptor" => ["Biota", "Dinosauria", "Theropoda", "Velociraptor"],
        "Triceratops" => ["Biota", "Dinosauria", "Ornithischia", "Triceratops"],
    )
    clear_cache()
    try
        cache = Taxodist._taxodist_cache
        for (index, (name, lineage)) in enumerate(sort(collect(lineages)))
            cache["id_$(name)"] = string(index)
            cache["lin_$(index)"] = lineage
            cache["ranks_$(index)"] = ["Clade", "Clade", "Order", "Genus"]
        end
        @test taxo_distance("Tyrannosaurus", "Triceratops"; metric="path").distance == 4
        @test taxo_distance("Tyrannosaurus", "Triceratops"; metric="path").metric == "path"
        @test taxo_distance("Tyrannosaurus", "Velociraptor"; metric="ranked").distance == 1
        @test taxo_distance("Tyrannosaurus", "Triceratops"; metric="ranked").distance == 2
        @test taxo_distance("Tyrannosaurus", "Triceratops").metric == "reciprocal"
        @test_throws ArgumentError taxo_distance("Tyrannosaurus", "Triceratops"; metric="x")

        relatives = closest_relative(
            "Tyrannosaurus", ["Triceratops", "Velociraptor"]; metric="path",
        )
        @test relatives.taxon == ["Velociraptor", "Triceratops"]
        @test relatives.distance == [2, 4]

        focal = focal_distances(
            "Tyrannosaurus", ["Triceratops", "Velociraptor", "Tyrannosaurus"];
            metric="wu_palmer", progress=false,
        )
        @test focal.data.distance ≈ [0, 1 - 6 / 8, 1 - 4 / 8]
        ranked_focal = focal_distances(
            "Tyrannosaurus", ["Triceratops"]; metric="ranked", progress=false,
        )
        @test ranked_focal.data.distance == [2]
        @test ranked_focal.data.mrca == ["Dinosauria"]

        names = ["Tyrannosaurus", "Velociraptor", "Triceratops"]
        by_name = distance_matrix(names; metric="ranked", progress=false)
        @test by_name["Tyrannosaurus", "Triceratops"] == 2
        @test Matrix(cross_distances(names, names; metric="ranked", progress=false)) ==
              Matrix(by_name)
    finally
        clear_cache()
    end
end

# -- Bundles (schema 1.1) -----------------------------------------------------

@testset "Bundles record the metric and ranks and survive a JSON round trip" begin
    res = ranked_example()
    names_map = Dict(
        "reciprocal" => "inverse_mrca_depth", "path" => "path_length",
        "wu_palmer" => "wu_palmer", "ranked" => "ranked_path_length",
    )
    mktempdir() do directory
        for metric in Taxodist.DISTANCE_METRICS
            bundle = taxo_bundle(res; metric=metric, progress=false)
            @test bundle.schema_version == "1.1"
            @test bundle.matrix.metric == metric
            @test bundle.metric["name"] == names_map[metric]
            @test haskey(bundle.metric, "rank_levels") == (metric == "ranked")
            file = joinpath(directory, "$(metric).json")
            write_taxodist_bundle(bundle, file)
            back = read_taxodist_bundle(file)
            @test back.metric["name"] == names_map[metric]
            @test back.matrix.metric == metric
            @test Matrix(back.matrix) == Matrix(bundle.matrix)
            @test isequal(back.resolution.ranks[1], res.ranks[1])
        end
        ranked = taxo_bundle(res; metric="ranked", progress=false)
        @test ranked.metric["rank_levels"] == collect(Taxodist.DEFAULT_RANK_LEVELS)
        @test_throws r"Unknown distance metric" taxo_bundle(res; metric="nope", progress=false)
    end
end

@testset "Bundle JSON keeps missing ranks as null" begin
    mktempdir() do directory
        file = joinpath(directory, "bundle.json")
        write_taxodist_bundle(taxo_bundle(ranked_example(); progress=false), file)
        raw = Taxodist.JSON3.read(read(file, String))
        @test raw["taxa"][1]["lineage_ranks"][1] === nothing
        @test raw["taxa"][1]["lineage_ranks"][2] == "Kingdom"
        @test raw["schema_version"] == "1.1"
        @test Set(String.(keys(raw["metric"]))) ==
              Set(["name", "definition", "root_depth", "common_ancestor"])
    end
end

@testset "Bundle validation recomputes distances with the stored metric" begin
    res = ranked_example()
    mktempdir() do directory
        file = joinpath(directory, "bundle.json")
        write_taxodist_bundle(taxo_bundle(res; metric="path", progress=false), file)
        raw = read(file, String)
        write(file, replace(raw, "\"path_length\"" => "\"inverse_mrca_depth\""))
        @test_throws r"stored distances" read_taxodist_bundle(file)
        write(file, replace(raw, "\"path_length\"" => "\"unknown_metric\""))
        @test_throws r"unknown distance metric" read_taxodist_bundle(file)

        ranked = taxo_bundle(res; metric="ranked", progress=false)
        delete!(ranked.metric, "rank_levels")
        @test_throws r"rank_levels" validate_taxodist_bundle(ranked; check_distances=false)

        # Ranks must align with lineages; the ranks of a missing lineage are ignored.
        broken = taxo_bundle(res; progress=false)
        broken.resolution.data.ranks[1] = ["Kingdom"]
        @test_throws r"lineage ranks" validate_taxodist_bundle(broken; check_distances=false)
        broken.resolution.data.ranks[1] = nothing
        @test_throws r"lineage ranks" validate_taxodist_bundle(broken; check_distances=false)
    end
end

@testset "Schema 1.0 bundles are still read and must use the reciprocal metric" begin
    res = taxo_from_lineages(["A" => ["Biota", "X", "A"], "B" => ["Biota", "Y", "B"]])
    mktempdir() do directory
        file = joinpath(directory, "bundle.json")
        write_taxodist_bundle(taxo_bundle(res; progress=false), file)
        raw = Taxodist.JSON3.read(read(file, String), Dict{String,Any})
        raw["schema_version"] = "1.0"
        for record in raw["taxa"]
            delete!(record, "lineage_ranks")
        end
        write(file, Taxodist.JSON3.write(raw))
        old = read_taxodist_bundle(file)
        @test old.schema_version == "1.0"
        @test !("ranks" in names(old.resolution.data))
        @test old.matrix.metric == "reciprocal"

        raw["metric"]["name"] = "path_length"
        raw["matrix"]["values"] = [[0, 4], [4, 0]]
        write(file, Taxodist.JSON3.write(raw))
        @test_throws r"only the reciprocal" read_taxodist_bundle(file)

        raw["schema_version"] = "2.0"
        write(file, Taxodist.JSON3.write(raw))
        @test_throws r"Unsupported" read_taxodist_bundle(file)
    end
end

@testset "Bundles written by R are read and reproduced" begin
    bundle = read_taxodist_bundle(metric_fixture("bundle_ranked_r.json"))
    @test bundle.schema_version == "1.1"
    @test bundle.software["language"] == "R"
    @test bundle.metric["name"] == "ranked_path_length"
    @test bundle.matrix.metric == "ranked"
    @test bundle.matrix.taxa == ["A", "B", "C", "D", "E"]
    @test Matrix(bundle.matrix) == R_METRIC_MATRICES["ranked"]
    @test isequal(bundle.resolution.ranks[1], STANDARD_RANKS)
    # The same analysis in Julia gives the same matrix.
    julia_bundle = taxo_bundle(ranked_example(); metric="ranked", progress=false)
    @test Matrix(julia_bundle.matrix) == Matrix(bundle.matrix)
    @test distance_matrix(bundle.resolution; metric="ranked", progress=false) ==
          bundle.matrix
end

# -- Taxonomic distinctness ----------------------------------------------------

@testset "taxo_distinctness reproduces the hand-computed indices" begin
    d = TaxonomicDistanceMatrix([0 1 3; 1 0 3; 3 3 0], ["a", "b", "c"])
    out = taxo_distinctness(Dict("a" => 2, "b" => 1, "c" => 1), d)
    row = out.data[1, :]
    # pairs: ab = 1 (x 2*1), ac = 3 (x 2*1), bc = 3 (x 1*1); n = 4
    @test row.sample == "sample"
    @test row.S == 3
    @test row.Delta ≈ (1 * 2 + 3 * 2 + 3 * 1) / 6
    @test row.Delta_star ≈ (1 * 2 + 3 * 2 + 3 * 1) / 5
    @test row.Delta_plus ≈ 7 / 3
    @test row.Lambda_plus ≈ mean(([1, 3, 3] .- 7 / 3) .^ 2)
    @test row.SDelta_plus ≈ 7
    @test out.metric === nothing
    @test names(out.data) ==
          ["sample", "S", "Delta", "Delta_star", "Delta_plus", "Lambda_plus", "SDelta_plus"]

    presence = taxo_distinctness(["a", "b", "c", "a"], d).data
    @test presence.S == [3]
    @test presence.Delta ≈ presence.Delta_plus
    @test presence.Delta_star ≈ presence.Delta_plus
    @test taxo_distinctness((a=2, b=1, c=1), d).data.Delta == out.data.Delta
    @test taxo_distinctness(["a" => 2, "b" => 1, "c" => 1], d).data.Delta == out.data.Delta
end

@testset "taxo_distinctness matches R (and vegan::taxondive)" begin
    res = ranked_example()
    d = distance_matrix(res; metric="ranked", progress=false)
    comm = [5 3 1 0 2; 1 0 2 4 0; 1 1 1 1 1]
    out = taxo_distinctness(comm, d; taxa=["A", "B", "C", "D", "E"], samples=["s1", "s2", "s3"])
    data = out.data
    @test out.metric == "ranked"
    @test data.sample == ["s1", "s2", "s3"]
    @test data.S == [4, 3, 5]
    # From R's taxo_distinctness().
    @test data.Delta ≈ [2.6727272727272728, 2.5714285714285716, 4.2999999999999998] rtol = 1e-12
    @test data.Delta_star ≈ [3.5853658536585367, 3.8571428571428572, 4.2999999999999998] rtol = 1e-12
    @test data.Delta_plus ≈ [4.166666666666667, 3.6666666666666665, 4.2999999999999998] rtol = 1e-12
    @test data.Lambda_plus ≈ [3.8055555555555554, 0.22222222222222224, 2.6099999999999999] rtol = 1e-12
    @test data.SDelta_plus ≈ [16.666666666666668, 11, 21.5] rtol = 1e-12

    frame = DataFrame(
        sample=["s1", "s2", "s3"], A=[5, 1, 1], B=[3, 0, 1], C=[1, 2, 1], D=[0, 4, 1], E=[2, 0, 1],
    )
    @test taxo_distinctness(frame, d).data == data
    unnamed = taxo_distinctness(frame[:, 2:end], d).data
    @test unnamed.sample == ["sample1", "sample2", "sample3"]
    @test unnamed.Delta == data.Delta
    @test taxo_distinctness(comm .> 0, d; taxa=["A", "B", "C", "D", "E"]).data.Delta ≈
          data.Delta_plus
end

@testset "taxo_distinctness handles small samples, data frames and bad input" begin
    d = TaxonomicDistanceMatrix([0 2; 2 0], ["a", "b"])
    one = taxo_distinctness(Dict("a" => 3, "b" => 0), d).data
    @test one.S == [1]
    @test isnan(one.Delta_plus[1])

    frame = DataFrame(sample=["x", "y"], a=[1, 0], b=[1, 1])
    out = taxo_distinctness(frame, d).data
    @test out.sample == ["x", "y"]
    @test out.Delta_plus[1] == 2
    @test isnan(out.Delta_plus[2])

    infinite = TaxonomicDistanceMatrix([0 Inf; Inf 0], ["a", "b"])
    @test isnan(taxo_distinctness(["a", "b"], infinite).data.Delta_plus[1])
    unresolved = TaxonomicDistanceMatrix([0 NaN; NaN 0], ["a", "b"])
    @test isnan(taxo_distinctness(["a", "b"], unresolved).data.Delta_plus[1])

    @test_throws r"no distances" taxo_distinctness(["a", "z"], d)
    @test_throws r"named by taxon" taxo_distinctness([1, 2], d)
    @test_throws r"non-negative" taxo_distinctness(Dict("a" => -1, "b" => 1), d)
    @test_throws r"non-negative" taxo_distinctness(Dict("a" => missing, "b" => 1), d)
    @test_throws r"must be" taxo_distinctness([[1]], d)
    @test_throws r"missing or empty" taxo_distinctness(["a", " "], d)
    @test_throws r"named by unique taxa" taxo_distinctness([1 2], d)
    @test_throws r"named by unique taxa" taxo_distinctness([1 2], d; taxa=["a", "a"])
    @test_throws r"TaxonomicDistanceMatrix" taxo_distinctness(["a", "b"], [0 2; 2 0])
    @test_throws r"one sample name" taxo_distinctness([1 2], d; taxa=["a", "b"], samples=["x", "y"])
end

@testset "taxo_distinctness computes distances when d is not supplied" begin
    clear_cache()
    try
        cache = Taxodist._taxodist_cache
        for (index, (name, lineage)) in enumerate(RANKED_LINEAGES)
            cache["id_$(name)"] = string(index)
            cache["lin_$(index)"] = lineage
            cache["ranks_$(index)"] = STANDARD_RANKS
        end
        d = distance_matrix(ranked_example(); metric="ranked", progress=false)
        expected = taxo_distinctness(["A", "B", "E"], d)
        got = taxo_distinctness(["A", "B", "E"]; progress=false)
        @test got.data == expected.data
        @test got.metric == "ranked"
        @test taxo_distinctness(["A", "B", "E"]; metric="path", progress=false).metric == "path"
    finally
        clear_cache()
    end
end

@testset "taxo_distinctness_test compares a sample with random draws from the pool" begin
    pool = distance_matrix(ranked_example(); metric="ranked", progress=false)
    test = taxo_distinctness_test(["A", "B", "C"], pool; n_iter=499, rng=Xoshiro(42))
    stats = test.statistics
    @test stats.statistic == ["Delta_plus", "Lambda_plus"]
    @test names(stats) ==
          ["statistic", "observed", "expected", "lower", "upper", "p_lower", "p_upper"]
    @test stats.observed[1] ≈ mean([1, 3, 3])
    # The expected Delta+ is the mean distance of the whole pool.
    full = Matrix(pool)
    @test stats.expected[1] ≈ mean(full[i, j] for j in 1:5 for i in (j + 1):5)
    @test stats.expected[2] ≈ mean(test.simulated[:, 2])
    @test all(stats.lower .<= stats.upper)
    @test all(0 .< stats.p_lower .<= 1)
    @test all(0 .< stats.p_upper .<= 1)
    @test size(test.simulated) == (499, 2)
    @test (test.S, test.pool_size, test.n_iter, test.metric) == (3, 5, 499, "ranked")
    # The three closest taxa form the least distinct sample.
    @test stats.p_lower[1] <= 0.2
    # Quantiles and p-values follow R's definitions.
    @test stats.lower[1] == quantile(test.simulated[:, 1], 0.025)
    @test stats.p_upper[1] ==
          (1 + count(>=(stats.observed[1] - sqrt(eps())), test.simulated[:, 1])) / 500

    a = taxo_distinctness_test(["A", "E"], pool; n_iter=50, rng=Xoshiro(1))
    b = taxo_distinctness_test(["A", "E"], pool; n_iter=50, rng=Xoshiro(1))
    @test a.simulated == b.simulated
    @test a.statistics == b.statistics
end

@testset "taxo_distinctness_test validates its input" begin
    pool = distance_matrix(ranked_example(); metric="ranked", progress=false)
    @test_throws r"master list" taxo_distinctness_test(["A", "Z"], pool)
    @test_throws r"between 2" taxo_distinctness_test(["A"], pool)
    @test_throws r"TaxonomicDistanceMatrix" taxo_distinctness_test(["A", "B"], Matrix(pool))
    @test_throws r"n_iter" taxo_distinctness_test(["A", "B"], pool; n_iter=0)
    @test_throws r"taxon names" taxo_distinctness_test([1, 2], pool)
    infinite = TaxonomicDistanceMatrix([0 Inf; Inf 0], ["a", "b"])
    @test_throws r"finite" taxo_distinctness_test(["a", "b"], infinite)
end

@testset "internal guards reject invalid input" begin
    @test_throws r"Unknown distance metric" Taxodist._metric_value(1, 2, 2, :nope)
    @test Taxodist._common_prefix_depth(["a", "b", "c"], ["a", "b", "d"]) == 2
    @test_throws r"missing or empty" Taxodist._distinctness_samples(Any["A", missing])
    @test_throws r"missing or empty" Taxodist._distinctness_samples(Any[nothing, "A"])
end
