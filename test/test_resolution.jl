html_response(body; status=200) = (status=status, body=Vector{UInt8}(codeunits(body)))

@testset "Lineage finalization" begin
    lineage = ["Biota", "Dinosauria", "[unranked] x", "Tyrannosaurus", "Tyrannosaurus rex"]
    expected = ["Biota", "Dinosauria", "Tyrannosaurus"]
    @test Taxodist._finalize_lineage(lineage, "Tyrannosaurus") == expected
    @test Taxodist._finalize_lineage(lineage, "tyrannosaurus") == expected
    @test Taxodist._finalize_lineage(["Biota", "Blattidae"], "Periplaneta") ==
          ["Biota", "Blattidae", "Periplaneta"]

    binomial = ["Biota", "Homo", "Homo sapiens", "Homo sapiens sapiens"]
    @test Taxodist._finalize_lineage(binomial, "homo sapiens") ==
          ["Biota", "Homo", "Homo sapiens"]
    @test Taxodist._finalize_lineage(["Biota", "Homo"], "Homo sapiens") ==
          ["Biota", "Homo", "Homo sapiens"]
end

@testset "Name normalization in get_lineage" begin
    clear_cache()
    Taxodist._taxodist_cache["id_Tyrannosaurus"] = "1"
    Taxodist._taxodist_cache["id_tyrannosaurus"] = "1"
    Taxodist._taxodist_cache["lin_1"] = ["Biota", "Dinosauria", "Tyrannosaurus"]
    expected = ["Biota", "Dinosauria", "Tyrannosaurus"]
    @test get_lineage("Tyrannosaurus") == expected
    @test get_lineage("  Tyrannosaurus ") == expected
    @test get_lineage("tyrannosaurus") == expected
    @test taxo_distance("Tyrannosaurus", "tyrannosaurus").distance == 0.0

    @test_throws ArgumentError get_lineage(["A", "B"])
    @test get_lineage(nothing) === nothing
    @test get_lineage(missing) === nothing
    @test get_lineage("   ") === nothing

    Taxodist._taxodist_cache["lin_50841"] = ["Biota", "Tyrannosaurus"]
    @test get_lineage(50841) == ["Biota", "Tyrannosaurus"]
    @test get_lineage_by_id(50841) == ["Biota", "Tyrannosaurus"]
    @test get_lineage_by_id(nothing) === nothing
    @test_throws ArgumentError get_lineage_by_id(["1", "2"])
    clear_cache()
end

@testset "Separate caches for cleaned and raw lineages" begin
    original = Taxodist._http_get[]
    clear_cache()
    calls = Ref(0)
    try
        page = """
        <html><body><div id="ctl00_divSubject"><b>Foo</b></div>
        <div id="divPageContent"><div>Natura</div><div>actualia</div>
        <div>Biota</div><div>Animalia</div><div>Foo</div></div></body></html>
        """
        Taxodist._http_get[] = (url, headers; kwargs...) -> begin
            calls[] += 1
            html_response(page)
        end
        raw = get_lineage_by_id("999"; clean=false)
        cleaned = get_lineage_by_id("999"; clean=true)
        @test raw == ["Natura", "actualia", "Biota", "Animalia", "Foo"]
        @test cleaned == ["Biota", "Animalia", "Foo"]
        @test calls[] == 1
    finally
        Taxodist._http_get[] = original
        clear_cache()
    end
end

@testset "Retries with backoff" begin
    original = Taxodist._http_get[]
    original_sleep = Taxodist._sleep[]
    pauses = Float64[]
    try
        Taxodist._sleep[] = seconds -> push!(pauses, seconds)
        responses = Any[ErrorException("offline"), html_response(""; status=503), html_response("ok")]
        Taxodist._http_get[] = (url, headers; kwargs...) -> begin
            next = popfirst!(responses)
            next isa Exception && throw(next)
            next
        end
        @test Taxodist._request_html("http://example.test") == "ok"
        @test pauses == [1, 2]

        Taxodist._http_get[] = (url, headers; kwargs...) -> html_response(""; status=404)
        empty!(pauses)
        @test Taxodist._request_html("http://example.test") === nothing
        @test isempty(pauses)
    finally
        Taxodist._http_get[] = original
        Taxodist._sleep[] = original_sleep
    end
end

@testset "ID lookup statuses and coverage" begin
    original = Taxodist._http_get[]
    clear_cache()
    search = """
    <table><tr><td>Alpha</td><td><a class="Valid" href="TaxonTree.aspx?id=101">tree</a></td></tr></table>
    """
    try
        Taxodist._http_get[] = (url, headers; kwargs...) -> html_response(""; status=503)
        status = @test_logs (:warn, r"Cannot reach") Taxodist._lookup_taxonomicon_id("Alpha")
        @test status == (id=nothing, status="retrieval_error")

        Taxodist._http_get[] = (url, headers; kwargs...) ->
            occursin("TaxonList", url) ? html_response(search) : html_response("<html></html>")
        @test Taxodist._lookup_taxonomicon_id("Alpha").status == "retrieval_error"

        Taxodist._http_get[] = (url, headers; kwargs...) -> occursin("TaxonList", url) ?
            html_response(search) :
            html_response("<a href='TaxonTree.aspx?id=1'>Natura</a><a href='TaxonTree.aspx?id=101'>Alpha</a>")
        @test Taxodist._lookup_taxonomicon_id("Alpha").status == "not_found"

        clear_cache()
        Taxodist._taxodist_cache["id_Found"] = "1"
        Taxodist._taxodist_cache["id_Absent"] = nothing
        Taxodist._http_get[] = (url, headers; kwargs...) -> html_response(""; status=503)
        coverage = @test_logs (:warn, r"Cannot reach") check_coverage(["Found", "Absent", "Offline"])
        @test coverage.covered[1] === true
        @test coverage.covered[2] === false
        @test ismissing(coverage.covered[3])
    finally
        Taxodist._http_get[] = original
        clear_cache()
    end
end

@testset "taxo_resolve consistency, normalization, and within" begin
    original = Taxodist._http_get[]
    clear_cache()
    searches = String[]
    try
        function resolve_http(url, headers; kwargs...)
            if occursin("TaxonList", url)
                push!(searches, url)
                body = if occursin("Periplaneta", url)
                    "<table><tr><td>Periplaneta</td><td><a class='Valid' href='TaxonTree.aspx?id=2'>t</a></td></tr></table>"
                elseif occursin("Blatta", url)
                    "<table><tr><td>Blatta</td><td><a class='Valid' href='TaxonTree.aspx?id=3'>t</a></td></tr></table>"
                else
                    """
                    <table>
                      <tr><td>Nereis worm</td><td><a class='Valid' href='TaxonTree.aspx?id=20'>t</a></td></tr>
                      <tr><td>Nereis plant</td><td><a class='Valid' href='TaxonTree.aspx?id=21'>t</a></td></tr>
                    </table>
                    """
                end
                return html_response(body)
            end
            body = if occursin("id=2&", url)
                "<a href='TaxonTree.aspx?id=1'>Biota</a><a href='TaxonTree.aspx?id=8'>Insecta</a><a href='TaxonTree.aspx?id=9'>Blattidae</a>"
            elseif occursin("id=3&", url)
                "<a href='TaxonTree.aspx?id=1'>Biota</a><a href='TaxonTree.aspx?id=8'>Insecta</a><a href='TaxonTree.aspx?id=9'>Blattidae</a><a href='TaxonTree.aspx?id=3'>Blatta</a>"
            elseif occursin("id=20&", url)
                "<a href='TaxonTree.aspx?id=1'>Biota</a><a href='TaxonTree.aspx?id=5'>Animalia</a><a href='TaxonTree.aspx?id=6'>Annelida</a><a href='TaxonTree.aspx?id=20'>Nereis</a>"
            else
                "<a href='TaxonTree.aspx?id=1'>Biota</a><a href='TaxonTree.aspx?id=7'>Plantae</a><a href='TaxonTree.aspx?id=21'>Nereis</a>"
            end
            return html_response(body)
        end
        Taxodist._http_get[] = resolve_http

        resolved = taxo_resolve(["Periplaneta", "Blatta", " Periplaneta "]; progress=false)
        @test resolved.lineage[1] == ["Biota", "Insecta", "Blattidae", "Periplaneta"]
        @test resolved.resolved_name[1] == "Periplaneta"
        @test resolved.input == ["Periplaneta", "Blatta", " Periplaneta "]
        @test count(url -> occursin("Periplaneta", url), searches) == 1
        by_resolution = distance_matrix(resolved; progress=false)
        by_name = distance_matrix(["Periplaneta", "Blatta"]; progress=false)
        @test by_resolution[1, 2] == by_name[1, 2] == 1 / 3

        animal = taxo_resolve(["Nereis"]; within="annelida", progress=false)
        @test animal.status == ["resolved"]
        @test animal.id == ["20"]
        fungus = taxo_resolve(["Nereis"]; within="Fungi", progress=false)
        @test fungus.status == ["unresolved"]
        @test fungus.n_candidates == [0]
        @test_throws ArgumentError taxo_resolve(["Nereis"]; within=" ", progress=false)
    finally
        Taxodist._http_get[] = original
        clear_cache()
    end
end

@testset "Persistent cache and cache_info names" begin
    clear_cache()
    Taxodist._taxodist_cache["lin_50841"] = ["Biota", "Tyrannosaurus"]
    Taxodist._taxodist_cache["lin_12345"] = ["Biota", "Velociraptor"]
    Taxodist._taxodist_cache["lin_999"] = String[]
    info = cache_info()
    @test info.ids == ["12345", "50841", "999"]
    @test info.taxa == ["Velociraptor", "Tyrannosaurus", nothing]

    mktempdir() do directory
        withenv("TAXODIST_CACHE_DIR" => joinpath(directory, "nested")) do
            path = save_cache()
            @test isfile(path)
            clear_cache()
            load_cache()
            @test Taxodist._taxodist_cache["lin_50841"] == ["Biota", "Tyrannosaurus"]
        end
        withenv(
            "TAXODIST_CACHE_DIR" => nothing,
            "LOCALAPPDATA" => directory,
            "XDG_CACHE_HOME" => directory,
        ) do
            @test Taxodist._platform() in (:windows, :apple, :unix)
            @test Taxodist._default_cache_file(:windows) ==
                  joinpath(directory, "taxodist", "Cache", "lineage-cache.json")
            @test Taxodist._default_cache_file(:apple) ==
                  joinpath(homedir(), "Library", "Caches", "taxodist", "lineage-cache.json")
            @test Taxodist._default_cache_file(:unix) ==
                  joinpath(directory, "taxodist", "lineage-cache.json")
        end
        withenv("TAXODIST_CACHE_DIR" => nothing, "XDG_CACHE_HOME" => nothing) do
            @test Taxodist._default_cache_file(:unix) ==
                  joinpath(homedir(), ".cache", "taxodist", "lineage-cache.json")
        end
    end
    @test_throws ArgumentError Taxodist._validated_cache(Dict("rawlin_1" => "Biota"))
    clear_cache()
end
