# Real TaxonTree and TaxonList pages, trimmed to the parts the parser reads,
# are shared with the R and Python implementations in test/fixtures.
# expected_lineages.json holds the raw lineage (before trimming at Biota)
# expected from each tree page. expected_text_fallback.json holds what R's
# fallback text parser returns for the same pages with the hierarchy block
# disabled, so the Julia text parser can be checked against R exactly.

const FIXTURES = joinpath(@__DIR__, "fixtures")
fixture_html(name) = read(joinpath(FIXTURES, name), String)
fixture_lineage(id) = Taxodist._lineage_from_page(fixture_html("tree_$(id).html"), id)
read_expected(name) = Dict(
    String(key) => String.(collect(value))
    for (key, value) in pairs(Taxodist.JSON3.read(fixture_html(name)))
)

@testset "Structured parser reproduces every fixture lineage" begin
    expected = read_expected("expected_lineages.json")
    @test length(expected) > 5
    for (id, lineage) in expected
        @test fixture_lineage(id) == lineage
    end
end

@testset "Structured parser keeps clean names despite authors and markers" begin
    passer = fixture_lineage("59107")
    @test all(name -> name in passer, ["Passeres", "Passerida", "Terrestrornithes"])
    @test !any(name -> occursin(r"[ ,(]", name), filter(!=("Passer"), passer))
    @test last(passer) == "Passer"
    @test !any(name -> startswith(name, "Passer "), passer)

    microtus = fixture_lineage("62669")
    @test "Supramyomorpha" in microtus
    @test last(microtus) == "Microtus"

    lipotes = fixture_lineage("68670")
    @test lipotes[(end - 2):end] == ["Lipotoidea", "Lipotidae", "Lipotes"]

    anamorph = fixture_lineage("204388")
    @test last(anamorph) == "Epidermophyton"
end

@testset "Structured parser appends the subject when its own row is missing" begin
    yi = fixture_lineage("7248674")
    @test yi[(end - 1):end] == ["Scansoriopterygidae", "Yi"]
    @test !("[plesion]" in yi)
end

@testset "Structured parser handles species pages" begin
    homo = fixture_lineage("66295")
    @test homo[(end - 1):end] == ["Homo", "Homo sapiens"]
end

@testset "Hierarchy rows: nesting, skipped classes, and missing bold" begin
    html = """
    <div id="ctl00_divSubject"><b>Gamma</b></div>
    <div id="ctl00_cphTabContent_divHierarchy"><p>
    <a href="x"><b>Alpha</b><span><b>Nested</b></span></a>
    <br>  <span class="Small"><b>Skipped</b></span><a><i>no bold</i></a>
    <a class="Ordinal"><b>Ordinal</b></a><a><b>Beta</b></a>
    <br>    <span class="Valid">Genus <b> Gamma   x </b></span>
    <br>      <a><b>Child</b></a>
    </p></div>
    """
    @test Taxodist._lineage_from_page(html, "1") == ["Alpha", "Beta", "Gamma x"]

    # Without the page's own row, keep the strictly increasing run of
    # indentation and append the subject (markers removed).
    missing_self = """
    <div id="ctl00_divSubject">Genus <b>†Delta?</b></div>
    <div id="ctl00_cphTabContent_divHierarchy">
    <a><b>Alpha</b></a><br>  <a><b>Beta</b></a>
    <br>  <a><b>Sibling</b></a></div>
    """
    @test Taxodist._lineage_from_page(missing_self, "1") == ["Alpha", "Beta", "Delta"]

    # An empty hierarchy block falls back to text parsing.
    empty_block = """
    <div id="ctl00_divSubject"><b>Beta</b></div>
    <div id="divPageContent"><div id="ctl00_cphTabContent_divHierarchy"></div>
    Natura
    Clade Beta Author
    </div>
    """
    @test Taxodist._lineage_from_page(empty_block, "1") == ["Natura", "Beta"]
end

@testset "Text fallback parser matches R on real pages" begin
    expected = read_expected("expected_text_fallback.json")
    for (id, lineage) in expected
        html = replace(
            fixture_html("tree_$(id).html"),
            "id=\"ctl00_cphTabContent_divHierarchy\"" => "id=\"hierarchy_removed\"";
            count=1,
        )
        @test Taxodist._lineage_from_page(html, id) == lineage
    end
end

@testset "get_lineage_by_id uses the structured parser and caches both variants" begin
    original = Taxodist._http_get[]
    clear_cache()
    html = fixture_html("tree_59107.html")
    try
        Taxodist._http_get[] = (url, headers; kwargs...) ->
            (status=200, body=Vector{UInt8}(codeunits(html)))
        cleaned = get_lineage_by_id("59107")
        raw = get_lineage_by_id("59107"; clean=false)
        @test first(cleaned) == "Biota"
        @test first(raw) == "Natura"
        @test raw[(end - length(cleaned) + 1):end] == cleaned
    finally
        Taxodist._http_get[] = original
        clear_cache()
    end
end

@testset "_finalize_lineage keeps name-ancestors of species and subspecies" begin
    finalize = Taxodist._finalize_lineage
    @test finalize(["Biota", "Homo", "Homo sapiens"], "Homo sapiens sapiens") ==
          ["Biota", "Homo", "Homo sapiens", "Homo sapiens sapiens"]
    @test finalize(
        ["Biota", "Microtus", "Microtus (Microtus)", "Microtus arvalis"], "Microtus arvalis",
    ) == ["Biota", "Microtus", "Microtus (Microtus)", "Microtus arvalis"]
    # A subgenus is a child of a queried genus and is dropped.
    @test finalize(["Biota", "Microtus", "Microtus (Microtus)"], "Microtus") ==
          ["Biota", "Microtus"]
    # Species of another genus are not name-ancestors.
    @test finalize(["Biota", "Homo", "Pan troglodytes"], "Homo sapiens") ==
          ["Biota", "Homo", "Homo sapiens"]
    # Truncation happens at the last match.
    @test finalize(["Biota", "Alpha", "Beta", "Alpha", "Child"], "Alpha") ==
          ["Biota", "Alpha", "Beta", "Alpha"]
end

@testset "Fallback text parser extracts names regardless of author format" begin
    extract = Taxodist._extract_scientific_name
    raw = [
        "Ichthyosauria de Blainville, 1835", "Foraminifera d'Orbigny, 1826",
        "Xenus Özdikmen, 2009", "Xenus O'Brien, 1990",
        "Xidae van der Hoeven, 1850", "Homo sapiens Linnaeus, 1758",
        "Homo sapiens sapiens Linnaeus, 1758",
        "Alismatales R.Br. ex Bercht. & J.Presl", "Firmicutes corrig. Gibbons",
        "Scombriformes - mackerels", "Romeriida Gauthier,",
        "Accipitres (C. Linnaeus, 1758)", "neomura Cavalier-Smith, 1987",
        "Pan-Aves Gauthier & de Queiroz, 2001",
    ]
    @test extract(raw) == [
        "Ichthyosauria", "Foraminifera", "Xenus", "Xenus", "Xidae",
        "Homo sapiens", "Homo sapiens sapiens", "Alismatales", "Firmicutes",
        "Scombriformes", "Romeriida", "Accipitres", "neomura", "Pan-Aves",
    ]
    # Quoted and bracketed lines are left for the caller's filters.
    @test extract(["\"Setaphyta\"", "[unranked]"]) == ["\"Setaphyta\"", "[unranked]"]
    @test ismissing(extract(missing))
    @test extract("") == ""
end

@testset "Fallback text parser strips double daggers and qualified markers" begin
    fake_html = "<html><body><div id=\"ctl00_divSubject\"><b>Lipotes</b></div>" *
                "<div id=\"divPageContent\">Natura - nature\nClade Biota Wagner 2004\n" *
                "Superfamily ‡?Lipotoidea de Muizon, 1988\n" *
                "[plesion] Family †Lipotidae Zhou, 1978\n" *
                "Genus ‡?Lipotesᵀ Miller, 1918\n</div></body></html>"
    @test Taxodist._lineage_from_page(fake_html, "1") ==
          ["Natura", "Biota", "Lipotoidea", "Lipotidae", "Lipotes"]

    # Rank prefixes require whitespace after the rank, and the qualified
    # prefix accepts only the ranks of R's list.
    @test Taxodist._clean_lineage_label("Orderella magna Author") == "Orderella magna"
    @test Taxodist._clean_lineage_label("[stem] Subspecies Homo sapiens idaltu White") ==
          "Homo sapiens idaltu"
    @test Taxodist._clean_lineage_label("[crown] Legion Alpha") == "Alpha"
    @test Taxodist._clean_lineage_label("Anamorph genus Epidermophyton Sabour.") ==
          "Epidermophyton"
    @test Taxodist._clean_lineage_label("[crown]Clade Alpha") == "[crown]Clade Alpha"

    # The cutoff is case-sensitive: a lower-case occurrence does not end the
    # lineage early.
    case_html = "<div id=\"ctl00_divSubject\"><b>Alpha</b></div><div id=\"divPageContent\">" *
                "Natura\nClade alpha-group\nClade Beta\nGenus Alpha\nGenus Child\n</div>"
    @test Taxodist._lineage_from_page(case_html, "1") ==
          ["Natura", "alpha-group", "Beta", "Alpha"]

    # Without #divPageContent, links are read up to the last link to the page.
    links_html = """
    <div id="ctl00_divSubject"><b>Target</b></div>
    <a href="TaxonTree.aspx?id=1">Biota</a>
    <a href="TaxonTree.aspx?id=6">Target</a>
    <a href="TaxonTree.aspx?id=2">Kingdom Animalia</a>
    <a href="TaxonTree.aspx?id=6&src=0">Target</a>
    <a href="TaxonTree.aspx?id=7">Child</a>
    """
    @test Taxodist._lineage_from_page(links_html, "6") == ["Biota", "Target", "Animalia"]
end

@testset "Regular-expression escaping and whole-word matching" begin
    for text in ["Aus (Bus) cus", "a.b", "x|y", "a[b]", "a{2}", "a+b", "a^b\$", "a\\b", "a*?"]
        @test occursin(Regex("^" * Taxodist._regex_escape(text) * "\$"), text)
    end
    @test Taxodist._contains_word("Microtus (Agricola)", "Microtus (Agricola)")
    @test Taxodist._contains_word("Aus (Bus) cus", "Aus (Bus)")
    @test !Taxodist._contains_word("Homopus", "Homo")
    @test Taxodist._contains_word("Émys", "Émys")
    @test !Taxodist._contains_word("xémys", "mys")
    @test !Taxodist._contains_word("Genus gallus", "Gallus")
    @test Taxodist._contains_word("Genus gallus", "Gallus"; ignore_case=true)
    @test Taxodist._contains_word("Genus Gallus Brisson", "Gallus")
    @test !Taxodist._contains_word("Gallusian clade", "Gallus")
end

@testset "taxo_resolve reports non-biological candidates as unresolved" begin
    original = Taxodist._http_get[]
    original_parser = Taxodist._search_parser[]
    clear_cache()
    try
        Taxodist._http_get[] = (url, headers; kwargs...) -> (status=200, body=UInt8[])
        Taxodist._taxodist_cache["lin_1"] = ["Natura", "Mundus", "Venus"]
        Taxodist._search_parser[] = html -> [(id="1", name="Venus - planet")]
        venus = taxo_resolve("Venus"; progress=false)
        @test venus.status == ["unresolved"]
        # Unresolved records carry no candidates, as bundle validation requires.
        @test venus.n_candidates == [0]

        Taxodist._search_parser[] = html -> [(id="abc", name="Broken")]
        broken = taxo_resolve("Broken"; progress=false)
        @test broken.status == ["retrieval_error"]
    finally
        Taxodist._http_get[] = original
        Taxodist._search_parser[] = original_parser
        clear_cache()
    end
end

@testset "Invalid or rejected homonyms are ranked after valid ones" begin
    descriptions = [
        "Amanita (Genus) Dill., nom. inval. see Amanita Pers.",
        "Amanita (Genus) Pers. (1797), nom. cons. (fungus)",
        "Amanita (Genus) Foo, nom. rej.",
    ]
    @test Taxodist._candidate_order(descriptions) == [2, 1, 3]
    @test Taxodist._candidate_order(["a", "b"]) == [1, 2]
    @test Taxodist._candidate_order(["x NOM. NUD. y", "x nom.illeg. y", "valid"]) == [3, 1, 2]

    original = Taxodist._http_get[]
    original_parser = Taxodist._search_parser[]
    clear_cache()
    try
        Taxodist._http_get[] = (url, headers; kwargs...) -> (status=200, body=UInt8[])
        Taxodist._taxodist_cache["lin_1"] = ["Biota", "Plantae", "Amanita"]
        Taxodist._taxodist_cache["lin_2"] = ["Biota", "Fungi", "Amanita"]
        Taxodist._search_parser[] = html -> [
            (id="1", name="Amanita (Genus) Dill., nom. inval."),
            (id="2", name="Amanita (Genus) Pers., nom. cons."),
        ]
        resolved = taxo_resolve(["Amanita"]; ambiguity="first", progress=false)
        @test resolved.id == ["2"]
        @test resolved.candidates[1].id == ["2", "1"]
        @test (@test_logs (:warn, r"Multiple valid") get_taxonomicon_id("Amanita")) == "2"
    finally
        Taxodist._http_get[] = original
        Taxodist._search_parser[] = original_parser
        clear_cache()
    end
end

@testset "Cache files from an older parser keep only taxon IDs" begin
    clear_cache()
    mktempdir() do directory
        file = joinpath(directory, "cache.json")
        write(file, """{"id_Passer": "59107", "lin_59107": ["Biota", "Passeres (C."],
                        "resolved_lineage_1_Passer": ["Biota", "Passer"]}""")
        @test_logs (:warn, r"older version") load_cache(file)
        @test Taxodist._taxodist_cache["id_Passer"] == "59107"
        @test !haskey(Taxodist._taxodist_cache, "lin_59107")
        @test !haskey(Taxodist._taxodist_cache, "resolved_lineage_1_Passer")

        clear_cache()
        write(file, "{}")
        @test_logs min_level=Base.CoreLogging.Warn load_cache(file)

        write(file, """{"id_A": "1", "meta_parser_version": 1}""")
        @test_logs min_level=Base.CoreLogging.Warn load_cache(file)
        @test !haskey(Taxodist._taxodist_cache, "meta_parser_version")

        clear_cache()
        Taxodist._taxodist_cache["lin_1"] = ["Biota", "A"]
        save_cache(file)
        @test Taxodist.JSON3.read(read(file, String))["meta_parser_version"] == 2
        clear_cache()
        @test_logs min_level=Base.CoreLogging.Warn load_cache(file)
        @test Taxodist._taxodist_cache["lin_1"] == ["Biota", "A"]
        @test !haskey(Taxodist._taxodist_cache, "meta_parser_version")
    end
    clear_cache()
end

@testset "Requests are throttled by MIN_INTERVAL" begin
    original = Taxodist.MIN_INTERVAL[]
    try
        Taxodist.MIN_INTERVAL[] = 0.2
        Taxodist._throttle()
        start = time()
        Taxodist._throttle()
        @test time() - start >= 0.15

        Taxodist.MIN_INTERVAL[] = 0.0
        start = time()
        Taxodist._throttle()
        @test time() - start < 0.1

        # Every HTTP attempt is throttled.
        Taxodist.MIN_INTERVAL[] = 0.2
        original_get = Taxodist._http_get[]
        try
            Taxodist._http_get[] = (url, headers; kwargs...) -> (status=200, body=UInt8[])
            Taxodist._request_html("http://example.test")
            start = time()
            Taxodist._request_html("http://example.test")
            @test time() - start >= 0.15
        finally
            Taxodist._http_get[] = original_get
        end
    finally
        Taxodist.MIN_INTERVAL[] = original
    end
end

@testset "Search fixtures are parsed like R's .parse_search_rows" begin
    parse(name) = Taxodist._parse_search_results(fixture_html("search_$(name).html"))
    homo = parse("Homo_sapiens")
    @test [row.id for row in homo] == ["66295"]
    @test homo[1].name == "Homo sapiens (Species) Linnaeus, 1758 - human (mammal)"

    nereis = parse("Nereis")
    @test [row.id for row in nereis] == ["1109123", "15133"]
    @test [row.name for row in nereis] == [
        "Nereis (Genus) Hübner, 1806, nom. illeg. see (Genus) Mechanitis Fabricius, 1807 (butterfly)",
        "Nereis (Genus) Linnaeus, 1748, nom. inval. see (Genus) Nereis Linnaeus, 1758 (animal)",
    ]

    tyrannosaurus = parse("Tyrannosaurus")
    @test [row.id for row in tyrannosaurus] == ["50841"]
    @test tyrannosaurus[1].name == "†Tyrannosaurus (Genus) Osborn, 1905 (dinosaur)"

    # The class must be exactly "Valid", and only the first such TaxonTree
    # link of a row is considered.
    rows = Taxodist._parse_search_results("""
    <table>
      <tr><td>Multi</td><td><a class="Valid Extra" href="TaxonTree.aspx?id=1">t</a></td></tr>
      <tr><td>NoId</td><td><a class="Valid" href="TaxonTree.aspx">t</a>
          <a class="Valid" href="TaxonTree.aspx?id=2">t</a></td></tr>
      <tr><td>Astronomy course</td><td><a class="Valid" href="TaxonTree.aspx?id=3">t</a></td></tr>
      <tr><td>Minor Planet </td><td><a class="Valid" href="TaxonTree.aspx?id=4">t</a></td></tr>
    </table>
    """)
    @test [row.id for row in rows] == ["3"]
end
