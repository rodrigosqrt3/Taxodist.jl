function seed_lineage(name, id, lineage)
    Taxodist._taxodist_cache["id_$(name)"] = String(id)
    Taxodist._taxodist_cache["lin_$(id)"] = String.(lineage)
end

@testset "Labeled distance matrices" begin
    clear_cache()
    seed_lineage("Alpha", "1", ["Biota", "Animalia", "CladeA", "Alpha"])
    seed_lineage("Beta", "2", ["Biota", "Animalia", "CladeA", "Beta"])
    seed_lineage("Gamma", "3", ["Biota", "Animalia", "CladeB", "Gamma"])

    matrix = distance_matrix(["Alpha", "Beta", "Gamma"]; progress=false)

    @test size(matrix) == (3, 3)
    @test Base.IndexStyle(TaxonomicDistanceMatrix) == Base.IndexCartesian()
    @test matrix.taxa == ["Alpha", "Beta", "Gamma"]
    @test [matrix[i, i] for i in 1:3] == zeros(3)
    @test matrix["Alpha", "Beta"] == 1 / 3
    @test matrix[1, 3] == 1 / 2
    @test Matrix(matrix) == transpose(Matrix(matrix))
    @test_throws KeyError matrix["Missing", "Beta"]
    @test_throws KeyError matrix["Alpha", "Missing"]
    @test size(distance_matrix(String[]; progress=false)) == (0, 0)
    duplicate = distance_matrix(["Alpha", "Alpha"]; progress=false)
    @test duplicate.taxa == ["Alpha", "Alpha"]
    @test duplicate[1, 2] == 0.0

    @test_throws ArgumentError TaxonomicDistanceMatrix(zeros(2, 3), ["A", "B"])
    @test_throws ArgumentError TaxonomicDistanceMatrix(zeros(2, 2), ["A"])
    @test TaxonomicDistanceMatrix(zeros(2, 2), ["A", "A"]).taxa == ["A", "A"]

    shown = sprint(show, MIME"text/plain"(), matrix)
    @test occursin("TaxonomicDistanceMatrix", shown)
    @test occursin("Alpha, Beta, Gamma", shown)
    empty_shown = sprint(
        show,
        MIME"text/plain"(),
        TaxonomicDistanceMatrix(zeros(0, 0), String[]),
    )
    @test !occursin("taxa:", empty_shown)

    verbose_output = mktemp() do path, io
        redirect_stdout(io) do
            distance_matrix(["Alpha", "Beta"]; verbose=true, progress=true)
        end
        seekstart(io)
        read(io, String)
    end
    @test occursin("Retrieving lineage 1/2: Alpha", verbose_output)

    clear_cache()
end

@testset "Unresolved taxa remain explicit" begin
    clear_cache()
    seed_lineage("Alpha", "1", ["Biota", "Animalia", "Alpha"])
    Taxodist._taxodist_cache["id_Missing"] = nothing

    matrix = distance_matrix(["Alpha", "Missing"]; progress=false)
    @test isnan(matrix["Alpha", "Missing"])
    @test matrix["Missing", "Missing"] == 0.0

    seed_lineage("Mineral", "3", ["Mineralia", "Mineral"])
    disconnected = distance_matrix(["Alpha", "Mineral"]; progress=false)
    @test isinf(disconnected["Alpha", "Mineral"])

    clear_cache()
end

@testset "Closest and focal distances" begin
    clear_cache()
    seed_lineage("Alpha", "1", ["Biota", "Animalia", "CladeA", "Alpha"])
    seed_lineage("Beta", "2", ["Biota", "Animalia", "CladeA", "Beta"])
    seed_lineage("Gamma", "3", ["Biota", "Animalia", "CladeB", "Gamma"])
    Taxodist._taxodist_cache["id_Missing"] = nothing

    closest = closest_relative("Alpha", ["Gamma", "Missing", "Beta"])
    @test closest.taxon == ["Beta", "Gamma", "Missing"]
    @test closest.distance[1:2] == [1 / 3, 1 / 2]
    @test isnan(closest.distance[3])
    @test size(closest_relative("Alpha", String[])) == (0, 2)
    @test closest_relative("Missing", ["Alpha"]) === nothing

    focal = focal_distances("Alpha", ["Gamma", "Alpha", "Missing", "Beta"])
    @test focal.focal == "Alpha"
    @test focal.data.taxon == ["Alpha", "Beta", "Gamma", "Missing"]
    @test focal.data.distance[1] == 0.0
    @test focal.data.mrca[1] == "Alpha"
    @test focal.data.mrca_depth[1] == 4
    @test isnan(focal.data.distance[4])
    @test focal_distances("Missing", ["Alpha"]) === nothing

    @test lineage_depth("Alpha") == 4
    @test lineage_depth("Missing") === nothing

    coverage = check_coverage(["Alpha", "Missing"])
    @test coverage.taxon == ["Alpha", "Missing"]
    @test coverage.covered == [true, false]

    @test taxo_distance("Alpha", "Beta").mrca == "CladeA"
    @test mrca("Alpha", "Beta") == "CladeA"
    @test taxo_distance("Missing", "Alpha") === nothing
    @test taxo_distance("Alpha", "Missing") === nothing
    @test mrca("Missing", "Alpha") === nothing

    clear_cache()
end

function example_analysis_matrix()
    return TaxonomicDistanceMatrix(
        [
            0.0 0.2 0.5
            0.2 0.0 0.3
            0.5 0.3 0.0
        ],
        ["A", "B", "C"],
    )
end

@testset "Hierarchical clustering" begin
    matrix = example_analysis_matrix()
    result = taxo_cluster(matrix)

    @test result.hclust !== nothing
    @test result.dist === matrix
    @test result.method == :average
    @test length(result.hclust.order) == 3
    @test sort(result.hclust.order) == [1, 2, 3]

    singleton = TaxonomicDistanceMatrix(zeros(1, 1), ["A"])
    @test taxo_cluster(singleton).hclust === nothing

    missing_matrix = TaxonomicDistanceMatrix([0.0 NaN; NaN 0.0], ["A", "B"])
    @test taxo_cluster(missing_matrix).hclust === nothing

    infinite_matrix = TaxonomicDistanceMatrix([0.0 Inf; Inf 0.0], ["A", "B"])
    @test taxo_cluster(infinite_matrix).hclust === nothing
    @test_throws ArgumentError taxo_cluster(matrix; method="unknown")
    @test taxo_cluster(matrix; method="SINGLE").method == :single
    @test_throws ArgumentError Taxodist._analysis_matrix(1)

    clear_cache()
    seed_lineage("A", "10", ["Biota", "A"])
    seed_lineage("B", "11", ["Biota", "B"])
    from_taxa = Taxodist._analysis_matrix(["A", "B"]; progress=false)
    @test from_taxa.taxa == ["A", "B"]
    clear_cache()
end

@testset "Principal coordinates analysis" begin
    matrix = example_analysis_matrix()
    result = taxo_ordinate(matrix; k=2)

    @test result.points !== nothing
    @test result.dist === matrix
    @test size(result.points) == (3, 3)
    @test names(result.points) == ["taxon", "PC1", "PC2"]
    @test result.points.taxon == ["A", "B", "C"]
    @test length(result.GOF) == 2
    @test length(result.eig) == 3

    two_taxa = TaxonomicDistanceMatrix([0.0 0.5; 0.5 0.0], ["A", "B"])
    reduced = taxo_ordinate(two_taxa; k=2)
    @test names(reduced.points) == ["taxon", "PC1"]

    @test_throws ArgumentError taxo_ordinate(matrix; k=0)
    @test_throws ArgumentError taxo_ordinate(matrix; k=1.5)
    @test_throws ArgumentError taxo_ordinate(matrix; k=true)
    @test_throws ArgumentError taxo_ordinate(matrix; k=Inf)

    singleton = TaxonomicDistanceMatrix(zeros(1, 1), ["A"])
    @test taxo_ordinate(singleton).points === nothing

    missing_matrix = TaxonomicDistanceMatrix([0.0 NaN; NaN 0.0], ["A", "B"])
    @test taxo_ordinate(missing_matrix).points === nothing

    infinite_matrix = TaxonomicDistanceMatrix([0.0 Inf; Inf 0.0], ["A", "B"])
    @test taxo_ordinate(infinite_matrix).points === nothing

    zero_matrix = TaxonomicDistanceMatrix(zeros(3, 3), ["A", "B", "C"])
    zero_result = taxo_ordinate(zero_matrix; k=2)
    @test names(zero_result.points) == ["taxon"]
    @test all(isnan, zero_result.GOF)
end

function offline_resolution_080()
    return taxo_from_lineages(
        Dict(
            "Alpha" => ["Biota", "Animalia", "Alpha"],
            "Beta" => ["Biota", "Animalia", "Beta"],
        );
        source="Curated study",
    )
end

@testset "Portable taxodist bundles" begin
    resolution = offline_resolution_080()
    bundle = taxo_bundle(resolution; progress=false)

    @test bundle isa TaxodistBundle
    @test bundle.schema_version == "1.0"
    @test bundle.source["name"] == "Curated study"
    @test bundle.source["url"] === nothing
    @test bundle.software["version"] == "0.8.0"
    @test bundle.software["language"] == "Julia"
    @test bundle.metric["name"] == "inverse_mrca_depth"
    @test bundle.resolution === resolution
    @test bundle.matrix isa TaxonomicDistanceMatrix
    @test distance_matrix(bundle) === bundle.matrix
    @test Taxodist._analysis_matrix(bundle) === bundle.matrix
    @test validate_taxodist_bundle(bundle) === bundle
    @test occursin("TaxodistBundle", sprint(show, MIME"text/plain"(), bundle))

    mktempdir() do directory
        file = joinpath(directory, "bundle.json")
        @test write_taxodist_bundle(bundle, file) == abspath(file)
        restored = read_taxodist_bundle(file)
        @test restored isa TaxodistBundle
        @test restored.source["name"] == "Curated study"
        @test restored.resolution.input == resolution.input
        @test restored.resolution.lineage == resolution.lineage
        @test restored.matrix.taxa == bundle.matrix.taxa
        @test restored.matrix.values == bundle.matrix.values

        compact = joinpath(directory, "compact.json")
        write_taxodist_bundle(bundle, compact; pretty=false)
        @test isfile(compact)

        empty_file = joinpath(directory, "empty.json")
        empty_bundle = taxo_bundle(
            taxo_from_lineages(Dict{String,Vector{String}}());
            progress=false,
        )
        write_taxodist_bundle(empty_bundle, empty_file)
        restored_empty = read_taxodist_bundle(empty_file)
        @test size(restored_empty.matrix) == (0, 0)
        @test size(restored_empty.resolution.data, 1) == 0
    end
end

@testset "Bundle missing values and portable infinities" begin
    rows = [
        (
            input="Alpha",
            resolved_name="Alpha",
            id="1",
            status="resolved",
            n_candidates=1,
            lineage_depth=2,
            lineage=["Biota", "Alpha"],
            candidates=Taxodist.DataFrame(id=["1"], name=["Alpha"]),
        ),
        (
            input="Missing",
            resolved_name=nothing,
            id=nothing,
            status="unresolved",
            n_candidates=0,
            lineage_depth=nothing,
            lineage=nothing,
            candidates=Taxodist.DataFrame(id=String[], name=String[]),
        ),
    ]
    resolution = Taxodist._make_resolution(
        rows;
        source="The Taxonomicon",
        source_url="http://taxonomicon.taxonomy.nl",
    )
    bundle = taxo_bundle(resolution; progress=false)
    @test isnan(bundle.matrix["Alpha", "Missing"])

    disconnected = taxo_bundle(taxo_from_lineages(Dict(
        "Animal" => ["Biota", "Animal"],
        "Mineral" => ["Natura", "Mineralia", "Mineral"],
    )); progress=false)
    @test isinf(disconnected.matrix["Animal", "Mineral"])

    mktempdir() do directory
        missing_file = joinpath(directory, "missing.json")
        write_taxodist_bundle(bundle, missing_file)
        raw_missing = Taxodist.JSON3.read(read(missing_file, String))
        @test raw_missing.matrix.values[1][2] === nothing
        restored_missing = read_taxodist_bundle(missing_file)
        @test restored_missing.resolution.id[2] === nothing
        @test restored_missing.resolution.lineage[2] === nothing

        infinite_file = joinpath(directory, "infinite.json")
        write_taxodist_bundle(disconnected, infinite_file)
        raw_infinite = Taxodist.JSON3.read(read(infinite_file, String))
        @test raw_infinite.matrix.values[1][2] == "Infinity"
        @test isinf(read_taxodist_bundle(infinite_file).matrix[1, 2])
    end
end

@testset "Bundle validation rejects inconsistent objects" begin
    bad(change, pattern) = begin
        bundle = taxo_bundle(offline_resolution_080(); progress=false)
        change(bundle)
        error = try
            validate_taxodist_bundle(bundle)
            nothing
        catch caught
            caught
        end
        @test error isa ArgumentError
        error isa Exception && @test occursin(pattern, sprint(showerror, error))
    end

    @test_throws ArgumentError validate_taxodist_bundle(Dict())
    bad(bundle -> delete!(bundle, "metric"), "required fields")
    bad(bundle -> (bundle.schema_version = "2.0"), "Unsupported")
    bad(bundle -> (bundle.resolution = Taxodist.DataFrame()), "wrong type")
    bad(bundle -> (bundle.resolution.status[1] = "mystery"), "unknown")
    bad(bundle -> (bundle.resolution.n_candidates[1] = 2), "count mismatch")
    bad(bundle -> (bundle.resolution.id[1] = nothing), "incomplete resolved")
    bad(bundle -> (bundle.resolution.lineage_depth[1] = 99), "depth mismatch")
    bad(bundle -> (bundle.matrix.taxa[1] = "Wrong"), "matrix labels")
    bad(bundle -> (bundle.matrix.values[1, 2] = 0.75), "stored distances")
end

@testset "Bundle JSON input validation and cross-language reading" begin
    mktempdir() do directory
        @test_throws ArgumentError read_taxodist_bundle(joinpath(directory, "missing.json"))

        invalid = joinpath(directory, "invalid.json")
        write(invalid, "not json")
        @test_throws ArgumentError read_taxodist_bundle(invalid)

        foreign = joinpath(directory, "r-bundle.json")
        write(foreign, """
        {
          "format": "taxodist_bundle",
          "schema_version": "1.0",
          "created_at": "2026-09-23 12:00:00 UTC",
          "source": {"name": "Curated study", "url": null, "retrieved_at": null},
          "software": {"name": "taxodist", "version": "0.8.0", "language": "R"},
          "metric": {
            "name": "inverse_mrca_depth",
            "definition": "0 for identical lineages; otherwise 1 / depth(MRCA)",
            "root_depth": 1,
            "common_ancestor": "continuous common lineage prefix"
          },
          "taxa": [{
            "input": "Alpha",
            "resolved_name": "Alpha",
            "id": "custom:Alpha",
            "status": "resolved",
            "n_candidates": 1,
            "lineage_depth": 2,
            "lineage": ["Biota", "Alpha"],
            "candidates": [{"id": "custom:Alpha", "name": "Alpha"}]
          }],
          "matrix": {"labels": ["Alpha"], "values": [[0]]}
        }
        """)
        restored = read_taxodist_bundle(foreign)
        @test restored.software["language"] == "R"
        @test restored.source["url"] === nothing
        @test restored.source["retrieved_at"] === nothing
        @test restored.resolution.retrieved_at === nothing
        @test restored.matrix[1, 1] == 0.0

        wrong_format = joinpath(directory, "wrong.json")
        write(wrong_format, "{\"format\":\"wrong\",\"schema_version\":\"1.0\"}")
        @test_throws ArgumentError read_taxodist_bundle(wrong_format)
    end
end