# Regression tests for Julia-specific bugs fixed in 0.9.1.

@testset "taxo_from_lineages preserves order and rejects ambiguous ids" begin
    ordered = taxo_from_lineages(
        ["Beta" => ["Biota", "Beta"], "Alpha" => ["Biota", "Alpha"]];
        ids=["b", "a"],
    )
    @test ordered.input == ["Beta", "Alpha"]
    @test ordered.id == ["b", "a"]
    @test ordered.resolved_name == ["Beta", "Alpha"]

    # A Dict has no defined order: labels are sorted and positional ids are
    # rejected, while ids keyed by label still work.
    unordered = Dict("Zeta" => ["Biota", "Zeta"], "Alpha" => ["Biota", "Alpha"], "Mu" => ["Biota", "Mu"])
    @test taxo_from_lineages(unordered).input == ["Alpha", "Mu", "Zeta"]
    error = try
        taxo_from_lineages(unordered; ids=["1", "2", "3"])
        nothing
    catch caught
        caught
    end
    @test error isa ArgumentError
    @test occursin("Positional ids", sprint(showerror, error))
    keyed = taxo_from_lineages(unordered; ids=Dict("Zeta" => "z", "Alpha" => "a", "Mu" => "m"))
    @test Dict(keyed.input .=> keyed.id) == Dict("Alpha" => "a", "Mu" => "m", "Zeta" => "z")

    @test_throws ArgumentError taxo_from_lineages(["A" => ["Root"], "A" => ["Root"]])
    @test_throws ArgumentError taxo_from_lineages([1 => ["Root"]])
    @test_throws ArgumentError taxo_from_lineages(Dict(1 => ["Root"]))
    @test_throws ArgumentError taxo_from_lineages(["A" => "Root"])
    @test size(taxo_from_lineages(Pair{String,Vector{String}}[]), 1) == 0
end

@testset "Numeric taxa are treated as IDs" begin
    clear_cache()
    Taxodist._taxodist_cache["lin_50841"] = ["Biota", "Dinosauria", "Tyrannosaurus"]
    Taxodist._taxodist_cache["lin_1"] = ["Biota", "Dinosauria", "Velociraptor"]
    Taxodist._taxodist_cache["lin_2"] = ["Mineralia", "Rock"]
    Taxodist._taxodist_cache["id_Gone"] = nothing

    @test get_lineage_by_id(1.0) == ["Biota", "Dinosauria", "Velociraptor"]
    @test get_lineage_by_id(50841.0) == get_lineage_by_id(50841)
    @test get_lineage_by_id(1.5) === nothing
    @test get_lineage_by_id(missing) === nothing
    @test get_lineage(50841.0) == ["Biota", "Dinosauria", "Tyrannosaurus"]

    result = taxo_distance(50841, 1.0)
    @test result.distance == 1 / 2
    @test result.taxon_a == "50841"
    @test result.taxon_b == "1"
    @test mrca(50841, 1) == "Dinosauria"
    @test shared_clades(50841.0, 1) == ["Biota", "Dinosauria"]
    @test compare_lineages(50841, 1).mrca_depth == 2
    @test taxo_path(50841, 1).taxon_a == "50841"
    @test filter_clade([50841, 2], "Dinosauria") == ["50841"]

    closest = closest_relative(50841, [2, "Gone", 1.0])
    @test closest.taxon == ["1", "2", "Gone"]
    @test closest.distance[1:2] == [1 / 2, Inf]
    @test isnan(closest.distance[3])

    # Unresolved taxa (NaN) are sorted after disconnected ones (Inf).
    focal = focal_distances(50841, ["Gone", 2, 50841.0]; progress=false)
    @test focal.focal == "50841"
    @test focal.data.taxon == ["50841", "2", "Gone"]
    @test isinf(focal.data.distance[2])
    @test isnan(focal.data.distance[3])

    progress_output = mktemp() do path, io
        redirect_stdout(io) do
            focal_distances(50841, [1]; verbose=true, progress=true)
        end
        seekstart(io)
        read(io, String)
    end
    @test occursin("Computing focal distance 1/1: 1", progress_output)

    matrix = distance_matrix([50841, 1.0, missing]; progress=false)
    @test matrix.taxa == ["50841", "1", "missing"]
    @test matrix[1, 2] == 1 / 2
    @test isnan(matrix[1, 3]) && isnan(matrix[3, 2])
    @test matrix[3, 3] == 0.0
    clear_cache()
end

@testset "taxo_resolve accepts a single name" begin
    clear_cache()
    Taxodist._taxodist_cache["lin_7"] = ["Biota", "Seven"]
    single = taxo_resolve("7"; progress=false)
    @test single isa TaxodistResolution
    @test single.status == ["resolved"]
    @test single.input == ["7"]
    clear_cache()
end

@testset "Clustering methods" begin
    matrix = TaxonomicDistanceMatrix([0.0 0.2 0.5; 0.2 0.0 0.3; 0.5 0.3 0.0], ["A", "B", "C"])
    for method in ("single", "complete", "average", "ward", "ward_presquared")
        @test taxo_cluster(matrix; method=method).hclust !== nothing
    end
    @test taxo_cluster(matrix; method="ward.D2").method == :ward
    @test taxo_cluster(matrix; method="ward.D").method == :ward_presquared
    @test taxo_cluster(matrix; method=:complete).method == :complete
    error = try
        taxo_cluster(matrix; method="median")
        nothing
    catch caught
        caught
    end
    @test error isa ArgumentError
    message = sprint(showerror, error)
    @test all(name -> occursin(name, message), ["single", "complete", "average", "ward"])
end

@testset "Bundle JSON numbers are validated" begin
    record(; n_candidates=1, lineage_depth=2, candidates="""[{"id": "1", "name": "A"}]""") = """
    {"format": "taxodist_bundle", "schema_version": "1.0", "created_at": null,
     "source": {"name": "x", "url": null, "retrieved_at": null},
     "software": {}, "metric": {"name": "inverse_mrca_depth"},
     "taxa": [{"input": "A", "resolved_name": "A", "id": "1", "status": "resolved",
               "n_candidates": $(n_candidates), "lineage_depth": $(lineage_depth),
               "lineage": ["Biota", "A"], "candidates": $(candidates)}],
     "matrix": {"labels": ["A"], "values": [[0]]}}
    """
    mktempdir() do directory
        file = joinpath(directory, "bundle.json")
        write(file, record())
        @test read_taxodist_bundle(file) isa TaxodistBundle
        for bad in (
            record(n_candidates="1.5"),
            record(n_candidates="true"),
            record(lineage_depth="2.5"),
            record(candidates="{}"),
        )
            write(file, bad)
            @test_throws ArgumentError read_taxodist_bundle(file)
        end
    end
    @test_throws ArgumentError Taxodist._matrix_from_bundle_json(Dict(
        "labels" => ["A", "B"], "values" => Any[Any[0, true], Any[true, 0]],
    ))
    @test_throws ArgumentError Taxodist._resolution_from_bundle_json(
        Dict("a" => 1), Dict("name" => "x", "url" => nothing, "retrieved_at" => nothing),
    )
end

@testset "Cache access is thread safe" begin
    @test Taxodist._cache_lock isa ReentrantLock
    original = Taxodist._http_get[]
    clear_cache()
    try
        Taxodist._http_get[] = (url, headers; kwargs...) -> begin
            id = match(r"id=([0-9]+)", url).captures[1]
            body = "<a href='TaxonTree.aspx?id=1'>Biota</a><a href='TaxonTree.aspx?id=$(id)'>T$(id)</a>"
            (status=200, body=Vector{UInt8}(codeunits(body)))
        end
        tasks = [Threads.@spawn get_lineage_by_id(string(id)) for id in 100:163]
        lineages = fetch.(tasks)
        @test all(lineages[index] == ["Biota", "T$(99 + index)"] for index in eachindex(lineages))
        @test count(key -> startswith(key, "lin_"), keys(Taxodist._taxodist_cache)) == 64
    finally
        Taxodist._http_get[] = original
        clear_cache()
    end
end

@testset "Plots is an optional extension" begin
    @test Base.get_extension(Taxodist, :TaxodistPlotsExt) !== nothing
    @test Taxodist._plots() isa Module
end
