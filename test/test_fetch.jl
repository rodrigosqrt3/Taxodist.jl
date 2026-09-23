@testset "Lineage label cleaning" begin
    examples = [
        "Subphylum † Craniata auct." => "Craniata",
        "Family † Dromaeosauridae Gauthier, 1986" => "Dromaeosauridae",
        "Clade Romeriida Gauthier," => "Romeriida",
        "Homo sapiens" => "Homo sapiens",
        "[crown] Clade Dinosauria" => "Dinosauria",
    ]

    for (raw, expected) in examples
        @test Taxodist._clean_lineage_label(raw) == expected
    end

    @test Taxodist._normalise_space("  Alpha\n  Beta  ") == "Alpha Beta"
    @test Taxodist._clean_lineage_label("Species Alpha A.B.") == "Alpha"
    @test Taxodist._clean_lineage_label("Species Alpha von Author") == "Alpha"
    @test Taxodist._clean_lineage_label("Species Alpha (Author)") == "Alpha"
    @test Taxodist._clean_lineage_label("Species Alpha (1985)") == "Alpha"
    @test Taxodist._clean_lineage_label("Species Alpha [note]") == "Alpha"
    @test Taxodist._clean_lineage_label("Species Alpha A.") == "Alpha"
    @test Taxodist._clean_lineage_label("Species Alpha (informal)") == "Alpha"
    @test Taxodist._clean_lineage_label("Species Alpha \"note\"") == "Alpha"
    @test Taxodist._clean_lineage_label("\"invalid") == ""
end

@testset "Search HTML parsing" begin
    html = """
    <html><body><table>
      <tr>
        <td>Carnotaurus - animal - dinosaur</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=12345&src=0">tree</a></td>
      </tr>
      <tr>
        <td>Carnotaurus - duplicate</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=12345&src=0">tree</a></td>
      </tr>
      <tr>
        <td>Venus - astronomical object</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=999&src=0">tree</a></td>
      </tr>
      <tr>
        <td>Invalid taxon</td>
        <td><a class="Invalid" href="TaxonTree.aspx?id=555&src=0">tree</a></td>
      </tr>
    </table></body></html>
    """

    results = Taxodist._parse_search_results(html)

    @test length(results) == 1
    @test results[1].id == "12345"
    @test occursin("Carnotaurus", results[1].name)

    prefixed = Taxodist._parse_search_results("""
    <table><tr><td>N | T | P | R | B | L Alpha</td>
    <td><a class="Valid" href="TaxonTree.aspx?id=7">tree</a></td></tr></table>
    """)
    @test occursin("Alpha", prefixed[1].name)
    @test !startswith(prefixed[1].name, "N | T")
    @test isempty(Taxodist._parse_search_results("<table><tr><td>none</td></tr></table>"))
    @test isempty(Taxodist._parse_search_results("""
    <table>
      <tr><td>wrong page</td><td><a class="Valid" href="Other.aspx?id=1">tree</a></td></tr>
      <tr><td>missing id</td><td><a class="Valid" href="TaxonTree.aspx">tree</a></td></tr>
    </table>
    """))
end

@testset "Lineage HTML parsing" begin
    html = """
    <html><body>
      <a href="TaxonTree.aspx?id=0&src=0">Natura</a>
      <a href="TaxonTree.aspx?id=1&src=0">Biota</a>
      <a href="TaxonTree.aspx?id=2&src=0">Kingdom Animalia</a>
      <a href="TaxonTree.aspx?id=3&src=0">Clade Dinosauria</a>
      <a href="TaxonTree.aspx?id=99&src=0">Carnotaurus</a>
      <a href="TaxonTree.aspx?id=100&src=0">SomeChild</a>
    </body></html>
    """

    clean = Taxodist._parse_lineage_html(html, "99"; clean=true)
    unclean = Taxodist._parse_lineage_html(html, "99"; clean=false)

    @test clean == ["Biota", "Animalia", "Dinosauria", "Carnotaurus"]
    @test unclean == ["Natura", "Biota", "Animalia", "Dinosauria", "Carnotaurus"]
    @test !("SomeChild" in clean)

    noisy = """
    <a href="Other.aspx?id=1">Ignored</a>
    <a href="TaxonTree.aspx">Ignored</a>
    <a href="TaxonTree.aspx?id=1">Biota</a>
    <a href="TaxonTree.aspx?id=2">Grade</a>
    <a href="TaxonTree.aspx?id=3">Population sample</a>
    <a href="TaxonTree.aspx?id=4">Clade Alpha</a>
    <a href="TaxonTree.aspx?id=5">Clade Alpha</a>
    <a href="TaxonTree.aspx?id=6">Target</a>
    """
    @test Taxodist._parse_lineage_html(noisy, "6") == ["Biota", "Alpha", "Target"]

    # Real Taxonomicon pages expose hierarchy information in the complete
    # content text; some intermediate nodes are plain text rather than links.
    structured = """
    <html><body>
      <div id="ctl00_divSubject"><b>Carnotaurus</b></div>
      <div id="divPageContent">
        <div>Natura</div>
        <div>Biota</div>
        <div>Kingdom Animalia</div>
        <div>Clade Dinosauria</div>
        <div>Clade Neotheropoda</div>
        <div><a href="TaxonTree.aspx?id=99">Genus Carnotaurus</a></div>
        <div>Child taxon that must be excluded</div>
      </div>
    </body></html>
    """
    @test Taxodist._parse_lineage_html(structured, "99") == [
        "Biota", "Animalia", "Dinosauria", "Neotheropoda", "Carnotaurus",
    ]

    @test Taxodist._contains_distinct_word("Genus Gallus Brisson", "Gallus")
    @test !Taxodist._contains_distinct_word("Gallusian clade", "Gallus")

    # NullNode is the third concrete HTMLNode type in Gumbo. It carries no
    # text, so the generic structured-text fallback must ignore it.
    io = IOBuffer()
    @test isnothing(Taxodist._write_structured_text(io, Taxodist.Gumbo.NullNode()))
    @test isempty(String(take!(io)))
end

@testset "HTTP request handling" begin
    original = Taxodist._http_get[]
    try
        function mock_ok(url, headers; kwargs...)
            @test headers == ["User-Agent" => "taxodist Julia package/0.8.0"]
            @test kwargs[:status_exception] === false
            @test kwargs[:retry] === false
            @test kwargs[:request_timeout] == 30
            return (status=200, body=Vector{UInt8}(codeunits("<html>ok</html>")))
        end
        Taxodist._http_get[] = mock_ok
        @test Taxodist._request_html("http://example.test") == "<html>ok</html>"

        Taxodist._http_get[] = (url, headers; kwargs...) ->
            (status=503, body=UInt8[])
        @test Taxodist._request_html("http://example.test"; verbose=true) === nothing

        function mock_error(url, headers; kwargs...)
            error("offline")
        end
        Taxodist._http_get[] = mock_error
        @test Taxodist._request_html("http://example.test"; verbose=true) === nothing
    finally
        Taxodist._http_get[] = original
    end
end

@testset "Retrieval API with deterministic HTTP" begin
    original = Taxodist._http_get[]
    clear_cache()
    try
        search_html = """
        <table>
          <tr><td>Other</td><td><a class="Valid" href="TaxonTree.aspx?id=10">tree</a></td></tr>
          <tr><td>Alpha</td><td><a class="Valid" href="TaxonTree.aspx?id=20">tree</a></td></tr>
        </table>
        """
        lineage_10 = """
        <a href="TaxonTree.aspx?id=1">Biota</a>
        <a href="TaxonTree.aspx?id=10">Other</a>
        """
        lineage_20 = """
        <a href="TaxonTree.aspx?id=1">Biota</a>
        <a href="TaxonTree.aspx?id=2">Animalia</a>
        <a href="TaxonTree.aspx?id=20">Alpha</a>
        """
        function mock_taxonomicon(url, headers; kwargs...)
            body = occursin("TaxonList", url) ? search_html :
                occursin("id=10", url) ? lineage_10 : lineage_20
            return (status=200, body=Vector{UInt8}(codeunits(body)))
        end
        Taxodist._http_get[] = mock_taxonomicon

        results = taxo_search("Alpha beta"; verbose=true)
        @test names(results) == ["id", "name"]
        @test results.id == ["10", "20"]
        @test get_taxonomicon_id("Alpha"; verbose=true) == "20"
        @test get_taxonomicon_id("Alpha"; verbose=true) == "20"
        @test get_lineage_by_id("20"; verbose=true) == ["Biota", "Animalia", "Alpha"]
        copied = get_lineage_by_id("20")
        push!(copied, "mutation")
        @test get_lineage_by_id("20")[end] == "Alpha"

        clear_cache()
        @test_logs (:warn, r"Multiple valid biological entries") get_taxonomicon_id("Ambiguous") == "10"

        Taxodist._http_get[] = (url, headers; kwargs...) ->
            (status=200, body=Vector{UInt8}(codeunits("<html></html>")))
        clear_cache()
        @test taxo_search("Missing"; verbose=true) === nothing
        @test get_taxonomicon_id("Missing") === nothing
        @test get_lineage_by_id("30") === nothing

        nonbiological_search = """
        <table><tr><td>Rock</td><td><a class="Valid" href="TaxonTree.aspx?id=40">tree</a></td></tr></table>
        """
        nonbiological_lineage = """
        <a href="TaxonTree.aspx?id=1">Mineralia</a>
        <a href="TaxonTree.aspx?id=40">Rock</a>
        """
        Taxodist._http_get[] = (url, headers; kwargs...) -> begin
            body = occursin("TaxonList", url) ? nonbiological_search : nonbiological_lineage
            (status=200, body=Vector{UInt8}(codeunits(body)))
        end
        clear_cache()
        @test get_taxonomicon_id("Rock") === nothing

        empty_lineage_search = """
        <table><tr><td>Empty</td><td><a class="Valid" href="TaxonTree.aspx?id=50">tree</a></td></tr></table>
        """
        Taxodist._http_get[] = (url, headers; kwargs...) -> begin
            body = occursin("TaxonList", url) ? empty_lineage_search : "<html></html>"
            (status=200, body=Vector{UInt8}(codeunits(body)))
        end
        clear_cache()
        @test get_taxonomicon_id("Empty") === nothing
        Taxodist._taxodist_cache["id_Empty"] = "50"
        @test get_lineage("Empty") === nothing

        Taxodist._http_get[] = (url, headers; kwargs...) ->
            (status=500, body=UInt8[])
        @test taxo_search("Offline") === nothing
    finally
        Taxodist._http_get[] = original
        clear_cache()
    end
end

@testset "Lineage retrieval and name filtering" begin
    clear_cache()
    @test get_lineage_by_id("") === nothing
    @test get_lineage_by_id("abc") === nothing

    Taxodist._taxodist_cache["lin_1"] = ["Natura", "Biota", "Animalia"]
    @test get_lineage("1") == ["Natura", "Biota", "Animalia"]

    Taxodist._taxodist_cache["id_Alpha"] = "2"
    Taxodist._taxodist_cache["lin_2"] = [
        "Biota", "Some clade", "[informal]", "Alpha", "Child",
    ]
    @test get_lineage("Alpha") == ["Biota", "Alpha"]

    Taxodist._taxodist_cache["id_Homo sapiens"] = "3"
    Taxodist._taxodist_cache["lin_3"] = ["Biota", "Animalia", "Homo", "Homo sapiens", "Child taxon"]
    @test get_lineage("Homo sapiens") == ["Biota", "Animalia", "Homo", "Homo sapiens"]

    Taxodist._taxodist_cache["id_Added"] = "4"
    Taxodist._taxodist_cache["lin_4"] = ["Biota", "Animalia"]
    @test get_lineage("Added") == ["Biota", "Animalia", "Added"]

    Taxodist._taxodist_cache["id_Missing"] = nothing
    @test get_lineage("Missing") === nothing
    clear_cache()
end

@testset "Cache management" begin
    Taxodist.clear_cache()
    empty_info = Taxodist.cache_info()
    @test empty_info.n_lineages == 0
    @test empty_info.n_ids == 0

    Taxodist._taxodist_cache["id_Carnotaurus"] = "12345"
    Taxodist._taxodist_cache["lin_12345"] = ["Biota", "Animalia", "Carnotaurus"]

    populated = Taxodist.cache_info()
    @test populated.n_lineages == 1
    @test populated.n_ids == 1

    mktempdir() do directory
        file = joinpath(directory, "cache.json")
        @test isnothing(Taxodist.save_cache(file))
        Taxodist.clear_cache()
        @test isnothing(Taxodist.load_cache(file))
        @test Taxodist._taxodist_cache["id_Carnotaurus"] == "12345"
        @test Taxodist._taxodist_cache["lin_12345"] == [
            "Biota", "Animalia", "Carnotaurus",
        ]
    end

    Taxodist.clear_cache()
end

@testset "Cache validation is transactional" begin
    Taxodist.clear_cache()
    Taxodist._taxodist_cache["id_existing"] = "1"

    mktempdir() do directory
        file = joinpath(directory, "invalid.json")
        write(file, "{\"id_A\": 123}")
        @test_throws ArgumentError Taxodist.load_cache(file)
        @test Taxodist._taxodist_cache == Dict{String,Any}("id_existing" => "1")

        @test_throws ArgumentError Taxodist.load_cache(joinpath(directory, "missing.json"))
        write(file, "[]")
        @test_throws ArgumentError Taxodist.load_cache(file)
    end

    @test_throws ArgumentError Taxodist._validated_cache(Dict("" => "value"))
    @test_throws ArgumentError Taxodist._validated_cache(Dict("id_A" => ""))
    @test_throws ArgumentError Taxodist._validated_cache(Dict("lin_1" => "Biota"))
    @test_throws ArgumentError Taxodist._validated_cache(Dict("lin_1" => Any["Biota", 1]))
    @test_throws ArgumentError Taxodist._validated_cache(
        Dict("resolved_lineage_1_A" => Any["Biota", 1]),
    )
    valid_resolved = Taxodist._validated_cache(
        Dict("resolved_lineage_1_A" => AbstractString["Biota", "Animalia", "A"]),
    )
    @test valid_resolved["resolved_lineage_1_A"] == ["Biota", "Animalia", "A"]
    @test valid_resolved["resolved_lineage_1_A"] isa Vector{String}
    @test Taxodist._validated_cache(Dict("other" => 1))["other"] == 1

    Taxodist.clear_cache()
end

@testset "Final resolved lineage cache" begin
    clear_cache()
    Taxodist._taxodist_cache["id_Alpha"] = "20"
    Taxodist._taxodist_cache["lin_20"] = ["Biota", "Animalia", "Alpha"]

    first = get_lineage("Alpha")
    Taxodist._taxodist_cache["lin_20"] = ["Biota", "Changed", "Alpha"]
    second = get_lineage("Alpha")

    @test first == ["Biota", "Animalia", "Alpha"]
    @test second == first
    @test haskey(Taxodist._taxodist_cache, "resolved_lineage_1_Alpha")
    clear_cache()
end

@testset "Auditable batch resolution" begin
    original = Taxodist._http_get[]
    clear_cache()
    search_calls = Ref(0)
    lineage_calls = Ref(0)
    try
        function resolution_http(url, headers; kwargs...)
            if occursin("TaxonList", url)
                search_calls[] += 1
                body = if occursin("Alpha", url)
                    """
                    <table><tr><td>Alpha</td><td><a class="Valid" href="TaxonTree.aspx?id=1">tree</a></td></tr></table>
                    """
                elseif occursin("Nereis", url)
                    """
                    <table>
                      <tr><td>Nereis one</td><td><a class="Valid" href="TaxonTree.aspx?id=2">tree</a></td></tr>
                      <tr><td>Nereis two</td><td><a class="Valid" href="TaxonTree.aspx?id=3">tree</a></td></tr>
                    </table>
                    """
                elseif occursin("Offline", url)
                    return (status=503, body=UInt8[])
                else
                    "<html></html>"
                end
                return (status=200, body=Vector{UInt8}(codeunits(body)))
            end

            lineage_calls[] += 1
            body = if occursin("id=1", url)
                "<a href='TaxonTree.aspx?id=10'>Biota</a><a href='TaxonTree.aspx?id=1'>Alpha</a>"
            elseif occursin("id=2", url)
                "<a href='TaxonTree.aspx?id=10'>Biota</a><a href='TaxonTree.aspx?id=11'>Animalia</a><a href='TaxonTree.aspx?id=2'>Nereis</a>"
            elseif occursin("id=3", url)
                "<a href='TaxonTree.aspx?id=10'>Biota</a><a href='TaxonTree.aspx?id=11'>Animalia</a><a href='TaxonTree.aspx?id=3'>Nereis</a>"
            elseif occursin("id=99", url)
                "<a href='TaxonTree.aspx?id=10'>Biota</a><a href='TaxonTree.aspx?id=99'>Direct</a>"
            else
                "<html></html>"
            end
            return (status=200, body=Vector{UInt8}(codeunits(body)))
        end
        Taxodist._http_get[] = resolution_http

        result = @test_logs (:warn, r"Ambiguous taxon names") taxo_resolve(
            ["Alpha", "Nereis", "Missing", "Offline", "Alpha"];
            progress=false,
        )
        @test result isa TaxodistResolution
        @test result.status == [
            "resolved", "ambiguous", "unresolved", "retrieval_error", "resolved",
        ]
        @test result.id == ["1", "2", nothing, nothing, "1"]
        @test result.n_candidates == [1, 2, 0, 0, 1]
        @test size(result.candidates[2], 1) == 2
        @test search_calls[] == 4
        @test lineage_calls[] == 3
        @test result.source == "The Taxonomicon"
        @test summary_counts(result) == (
            resolved=2,
            ambiguous=1,
            unresolved=1,
            retrieval_error=1,
        )

        first = taxo_resolve(["Nereis"]; ambiguity="first", progress=false)
        @test first.status == ["ambiguous"]
        @test first.id == ["2"]
        @test_throws ArgumentError taxo_resolve(
            ["Nereis"];
            ambiguity="error",
            progress=false,
        )

        direct = taxo_resolve(["99"]; progress=false)
        failed = taxo_resolve(["404"]; progress=false)
        @test direct.resolved_name == ["Direct"]
        @test direct.id == ["99"]
        @test failed.status == ["retrieval_error"]
    finally
        Taxodist._http_get[] = original
        clear_cache()
    end

    @test_throws ArgumentError taxo_resolve([1, 2]; progress=false)
    @test_throws ArgumentError taxo_resolve(["Alpha", ""]; progress=false)
    @test_throws ArgumentError taxo_resolve(["Alpha"]; ambiguity="guess", progress=false)
end

@testset "Offline resolutions" begin
    resolution = taxo_from_lineages(
        Dict(
            "Alpha" => ["Biota", "Animalia", "Alpha"],
            "Beta" => ["Biota", "Animalia", "Beta"],
        );
        source="Curated study",
    )
    @test resolution isa TaxodistResolution
    @test Set(resolution.id) == Set(["custom:Alpha", "custom:Beta"])
    @test resolution.lineage_depth == [3, 3]
    @test resolution.source == "Curated study"
    @test resolution.source_url === nothing
    shown = sprint(show, MIME"text/plain"(), resolution)
    @test occursin("TaxodistResolution", shown)

    matrix = distance_matrix(resolution; progress=false)
    @test matrix["Alpha", "Beta"] == 1 / 2

    named = taxo_from_lineages(
        Dict("Alpha" => ["Root", "Alpha"], "Beta" => ["Root", "Beta"]);
        ids=Dict("Beta" => "B", "Alpha" => "A"),
    )
    identifier_by_input = Dict(named.input .=> named.id)
    @test identifier_by_input == Dict("Alpha" => "A", "Beta" => "B")

    @test size(taxo_from_lineages(Dict{String,Vector{String}}()).data, 1) == 0
    @test_throws ArgumentError taxo_from_lineages(Dict("" => ["Root"]))
    @test_throws ArgumentError taxo_from_lineages(Dict("Alpha" => String[]))
    @test_throws ArgumentError taxo_from_lineages(Dict("Alpha" => ["Root", ""]))
    @test_throws ArgumentError taxo_from_lineages(
        Dict("Alpha" => ["Root"]);
        source="",
    )
    @test_throws ArgumentError taxo_from_lineages(
        Dict("Alpha" => ["Root"]);
        ids=String[],
    )
    @test_throws ArgumentError taxo_from_lineages(
        Dict("Alpha" => ["Root"], "Beta" => ["Root"]);
        ids=["same", "same"],
    )
end