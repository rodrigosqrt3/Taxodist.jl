# Open Tree Taxonomy and NCBI Taxonomy sources. Requests are answered from
# real responses stored in fixtures/sources, shared with R and Python (see
# tools/record_source_fixtures.R in the R package); a request without a
# stored response behaves as a failed request (a network error).

using Taxodist: JSON3, HTTP
using Taxodist.Random: Xoshiro
using Taxodist.DataFrames: nrow

const SOURCE_FIXTURES = joinpath(@__DIR__, "fixtures", "sources")

source_response(text::AbstractString; status=200) = (status=status, body=Vector{UInt8}(codeunits(text)))

function source_fixture_response(file)
    path = joinpath(SOURCE_FIXTURES, file)
    isfile(path) || throw(Base.IOError("no stored response for $(file)", 0))
    return source_response(read(path, String))
end

function url_parameter(url, name)
    found = match(Regex("[?&]$(name)=([^&]*)"), url)
    return HTTP.URIs.unescapeuri(found.captures[1])
end

function fixture_http(method, url, headers, body; kwargs...)
    endswith(url, "/taxonomy/about") && return source_fixture_response("ott_about.json")
    endswith(url, "/tnrs/match_names") && return source_fixture_response("ott_match_names.json")
    if endswith(url, "/taxonomy/taxon_info")
        return source_fixture_response("ott_taxon_info_$(JSON3.read(body)[:ott_id]).json")
    end
    if occursin("esearch.fcgi", url)
        term = url_parameter(url, "term")
        name = lowercase(match(r"^\"(.*)\"\[", term).captures[1])
        field = occursin("[Scientific Name]", term) ? "scientific" : "all"
        return source_fixture_response(
            "ncbi_esearch_$(replace(name, r"[^a-z0-9]+" => "_"))_$(field).json",
        )
    end
    if occursin("efetch.fcgi", url)
        taxa = map(split(url_parameter(url, "id"), ",")) do id
            path = joinpath(SOURCE_FIXTURES, "ncbi_taxon_$(id).xml")
            isfile(path) || return ""
            return match(r"(?s)<Taxon>.*</Taxon>", read(path, String)).match
        end
        return source_response("<TaxaSet>" * join(taxa) * "</TaxaSet>")
    end
    throw(Base.IOError("no stored response for $(url)", 0))
end

failing_http(method, url, headers, body; kwargs...) = throw(Base.IOError("offline", 0))

# Answers source requests with `request` while running `f`, starting from an
# empty cache.
function with_sources(f, request=fixture_http)
    original = Taxodist._http_request[]
    original_sleep = Taxodist._throttle_sleep[]
    clear_cache()
    try
        Taxodist._http_request[] = request
        # NCBI requests are always spaced (0.4 s); skip the pauses.
        Taxodist._throttle_sleep[] = seconds -> nothing
        return f()
    finally
        Taxodist._http_request[] = original
        Taxodist._throttle_sleep[] = original_sleep
        clear_cache()
    end
end

function capture_stdout(f)
    path, io = mktemp()
    try
        redirect_stdout(f, io)
        close(io)
        return read(path, String)
    finally
        isopen(io) && close(io)
        rm(path; force=true)
    end
end

const SOURCE_NAMES = [
    "Homo", "Pan", "Mus", "Drosophila", "Nereis", "Quercus",
    "Fakeosaurus", "human", "Homo sapiens",
]
const SOURCE_IDS = Dict("ott" => ["ott:770309", "417957"], "ncbi" => ["ncbi:9606", "9598"])

json_or_nothing(value) = value === nothing ? nothing : value

@testset "Resolutions reproduce the recorded expectations (shared fixture)" begin
    expected = JSON3.read(read(joinpath(SOURCE_FIXTURES, "expected_resolutions.json"), String))
    @test sort(String.(collect(keys(expected)))) ==
          ["ncbi", "ncbi_within_metazoa", "ott", "ott_within_metazoa"]
    for (raw_key, exp) in pairs(expected)
        key = String(raw_key)
        source = replace(key, r"_within_metazoa$" => "")
        within = occursin("within", key) ? "Metazoa" : nothing
        with_sources() do
            res = taxo_resolve(
                vcat(SOURCE_NAMES, SOURCE_IDS[source]);
                source=source, progress=false, ambiguity="first", within=within,
            )
            @test res.source == exp[:source]
            @test res.source_url == exp[:source_url]
            @test length(res) == length(exp[:records])
            for (i, rec) in enumerate(exp[:records])
                @testset "$(key) $(rec[:input])" begin
                    @test res.input[i] == rec[:input]
                    @test res.status[i] == rec[:status]
                    @test res.id[i] == json_or_nothing(rec[:id])
                    @test res.resolved_name[i] == json_or_nothing(rec[:resolved_name])
                    @test res.n_candidates[i] == rec[:n_candidates]
                    @test res.lineage_depth[i] == json_or_nothing(rec[:lineage_depth])
                    if rec[:lineage] === nothing
                        @test res.lineage[i] === nothing
                        @test res.ranks[i] === nothing
                    else
                        @test res.lineage[i] == String.(rec[:lineage])
                        @test isequal(
                            res.ranks[i],
                            [rank === nothing ? missing : rank for rank in rec[:ranks]],
                        )
                    end
                    @test res.candidates[i].id == String[candidate[:id] for candidate in rec[:candidates]]
                    @test res.candidates[i].name == String[candidate[:name] for candidate in rec[:candidates]]
                end
            end
        end
    end
end

@testset "External sources resolve the expected taxa" begin
    with_sources() do
        ott = taxo_resolve(["Homo", "Nereis", "Fakeosaurus"]; source="ott", progress=false)
        @test isequal(ott.id, ["ott:770309", "ott:305738", nothing])
        @test ott.status == ["resolved", "resolved", "unresolved"]
        @test ott.lineage[1][(end - 2):end] == ["Hominidae", "Homininae", "Homo"]
        @test ott.lineage[1][1] == "life"
        @test last(ott.ranks[1]) == "Genus"
        @test startswith(ott.source, "Open Tree Taxonomy 3.7")

        ncbi = @test_logs (:warn, r"Drosophila \(2 candidates\)") taxo_resolve(
            ["Mus", "Drosophila", "human"]; source=:ncbi, progress=false,
        )
        # The nominotypical subgenus Mus (Mus) is not counted as a homonym.
        @test ncbi.id[1] == "ncbi:10088"
        @test ncbi.status[1] == "resolved"
        # Drosophila is a true homonym (a fly and a fungus); the fly has the lower ID.
        @test ncbi.id[2] == "ncbi:7215"
        @test ncbi.status[2] == "ambiguous"
        # A common name is matched through "All Names".
        @test ncbi.resolved_name[3] == "Homo sapiens"
        @test ncbi.lineage[1][1] == "cellular organisms"
        @test ncbi.source == "NCBI Taxonomy"
        @test ncbi.source_url == "https://www.ncbi.nlm.nih.gov/taxonomy"

        # A single name, and identifiers with or without the prefix.
        @test taxo_resolve("Pan"; source="ott", progress=false).id == ["ott:417957"]
        ids = taxo_resolve(["9598", "ncbi:9606"]; source="ncbi", progress=false)
        @test ids.id == ["ncbi:9598", "ncbi:9606"]
        @test ids.resolved_name == ["Pan troglodytes", "Homo sapiens"]
        # An identifier too large to be an OTT id cannot be retrieved.
        @test taxo_resolve(["123456789012345678901234567890"]; source="ott", progress=false).status ==
              ["retrieval_error"]
        output = capture_stdout(() -> taxo_resolve(["Homo"]; source="ott"))
        @test occursin("Resolving taxa [1/1]: Homo", output)
    end
end

@testset "Resolutions from other sources feed every matrix function" begin
    with_sources() do
        res = taxo_resolve(["Homo", "Pan", "Mus", "Quercus"]; source="ott", progress=false)
        d = distance_matrix(res; progress=false)
        @test d isa TaxonomicDistanceMatrix
        @test d["Homo", "Pan"] < d["Homo", "Mus"] < d["Homo", "Quercus"]
        @test taxo_tree(res; progress=false) isa TaxodistTree

        bundle = taxo_bundle(res; progress=false)
        @test startswith(bundle.source["name"], "Open Tree Taxonomy")
        mktempdir() do directory
            file = joinpath(directory, "bundle.json")
            write_taxodist_bundle(bundle, file)
            back = read_taxodist_bundle(file)
            @test back.source["name"] == bundle.source["name"]
            @test back.resolution.id == res.id
            @test back.matrix.values == d.values
        end

        bundle2 = taxo_bundle(["Homo", "Pan"]; source="ncbi", progress=false)
        @test bundle2.source["name"] == "NCBI Taxonomy"
        @test bundle2.resolution.id == ["ncbi:9605", "ncbi:9596"]
    end
end

@testset "Failed requests are retrieval errors" begin
    with_sources(failing_http) do
        ott = taxo_resolve(["Homo", "770309"]; source="ott", progress=false)
        @test ott.status == ["retrieval_error", "retrieval_error"]
        @test ott.source == "Open Tree Taxonomy"
        @test ott.source_url == "https://tree.opentreeoflife.org/about/taxonomy-version"
        ncbi = taxo_resolve("Homo"; source="ncbi", progress=false)
        @test ncbi.status == ["retrieval_error"]
    end
end

@testset "A candidate whose lineage cannot be retrieved is a retrieval error" begin
    function no_lineages(method, url, headers, body; kwargs...)
        occursin("taxon_info", url) && return source_response(""; status=500)
        return fixture_http(method, url, headers, body)
    end
    with_sources(no_lineages) do
        res = taxo_resolve("Homo"; source="ott", progress=false)
        @test res.status == ["retrieval_error"]
        @test res.candidates[1].id == ["ott:770309"]
        @test res.n_candidates == [1]
        inside = taxo_resolve("Homo"; source="ott", progress=false, within="Metazoa")
        @test inside.status == ["retrieval_error"]
        @test inside.candidates[1].id == ["ott:770309"]
    end
end

@testset "within restricts candidates of external sources" begin
    with_sources() do
        plants = taxo_resolve(["Homo", "Quercus"]; source="ott", progress=false, within="Chloroplastida")
        @test plants.status == ["unresolved", "resolved"]
        @test plants.n_candidates == [0, 1]
    end
end

@testset "taxo_sources describes every source" begin
    sources = taxo_sources()
    @test sources.id == ["taxonomicon", "ott", "ncbi"]
    @test all(!isempty, sources.citation)
    @test names(sources) == ["id", "name", "url", "id_prefix", "ranks", "licence", "citation"]
    @test sources.id_prefix == ["", "ott:", "ncbi:"]
end

@testset "Sources are validated" begin
    @test_throws r"Unknown taxonomic source \"gbif\"" taxo_resolve("Homo"; source="gbif")
    @test_throws r"single string" taxo_resolve("Homo"; source=["ott", "ncbi"])
    @test_throws r"single string" taxo_resolve("Homo"; source=missing)
    @test_throws r"Unknown taxonomic source" taxo_bundle(["Homo"]; source="gbif")
end

@testset "Source ranks follow the Taxonomicon convention" begin
    ranks = [
        Taxodist._normalize_source_rank(rank) for rank in
        ["genus", "no rank", "clade", "cellular root", "", nothing, missing, "Superfamily",
         "no rank - terminal", " acellular root "]
    ]
    @test isequal(ranks, ["Genus", missing, "Clade", missing, missing, missing, missing,
                          "Superfamily", missing, missing])
end

@testset "Nested candidates of the same name are not homonyms" begin
    outer = Taxodist._outermost_candidates
    @test outer(Any[["a", "Mus"], ["a", "Mus", "Mus (Mus)"], ["b", "Mus"], nothing]) ==
          [true, false, true, true]
    @test outer(Any[["a", "x"]]) == [true]
end

@testset "The Taxonomicon ignores nominotypical subgenera" begin
    original = Taxodist._http_get[]
    clear_cache()
    search = """
    <table>
      <tr><td>Mus (Mus)</td><td><a class='Valid' href='TaxonTree.aspx?id=2'>t</a></td></tr>
      <tr><td>Mus</td><td><a class='Valid' href='TaxonTree.aspx?id=1'>t</a></td></tr>
    </table>
    """
    try
        Taxodist._http_get[] = (url, headers; kwargs...) -> source_response(search)
        Taxodist._taxodist_cache["lin_1"] = ["Biota", "Rodentia", "Mus"]
        Taxodist._taxodist_cache["lin_2"] = ["Biota", "Rodentia", "Mus", "Mus (Mus)"]
        res = @test_logs taxo_resolve("Mus"; progress=false)
        @test res.status == ["resolved"]
        @test res.id == ["1"]
        @test res.n_candidates == [1]
        @test res.source == "The Taxonomicon"
        @test (@test_logs get_taxonomicon_id("Mus")) == "1"
    finally
        Taxodist._http_get[] = original
        clear_cache()
    end
end

@testset "The NCBI request rate depends on the API key" begin
    original_key = Taxodist.NCBI_API_KEY[]
    original_interval = Taxodist.MIN_INTERVAL[]
    try
        Taxodist.MIN_INTERVAL[] = 0.0
        Taxodist.NCBI_API_KEY[] = ""
        @test Taxodist._source_interval("ncbi") == 0.4
        Taxodist.NCBI_API_KEY[] = "abc"
        @test Taxodist._source_interval("ncbi") == 0.12
        @test occursin("api_key=abc", Taxodist._ncbi_url("esearch.fcgi", ["db" => "taxonomy"]))
        @test Taxodist._source_interval("ott") == 0
        Taxodist.MIN_INTERVAL[] = 0.5
        @test Taxodist._source_interval("ncbi") == 0.5
        Taxodist.MIN_INTERVAL[] = NaN
        @test Taxodist._source_interval("ott") == 0
        Taxodist.NCBI_API_KEY[] = missing
        @test Taxodist._ncbi_api_key() == ""
        Taxodist.NCBI_API_KEY[] = nothing
        withenv("NCBI_API_KEY" => "fromenv") do
            @test Taxodist._ncbi_api_key() == "fromenv"
        end
        withenv("NCBI_API_KEY" => nothing) do
            @test Taxodist._ncbi_api_key() == ""
            url = Taxodist._ncbi_url("esearch.fcgi", ["term" => "\"Homo sapiens\"[All Names]"])
            @test url == "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" *
                         "term=%22Homo%20sapiens%22%5BAll%20Names%5D&tool=taxodist"
        end
    finally
        Taxodist.NCBI_API_KEY[] = original_key
        Taxodist.MIN_INTERVAL[] = original_interval
    end
end

@testset "Requests to each source are throttled independently" begin
    original_sleep = Taxodist._throttle_sleep[]
    waits = Float64[]
    try
        Taxodist._throttle_sleep[] = seconds -> push!(waits, seconds)
        empty!(Taxodist._last_requests)
        Taxodist._throttle("ncbi", 10.0)
        Taxodist._throttle("ott", 10.0)
        @test isempty(waits)
        Taxodist._throttle("ncbi", 10.0)
        @test length(waits) == 1 && 9 < waits[1] <= 10
        Taxodist._throttle("ott", 0.0)
        Taxodist._throttle("ott", NaN)
        @test length(waits) == 1
    finally
        Taxodist._throttle_sleep[] = original_sleep
        empty!(Taxodist._last_requests)
    end
end

@testset "External requests are built, retried, and parsed" begin
    original = Taxodist._http_request[]
    original_sleep = Taxodist._sleep[]
    pauses = Float64[]
    captured = Any[]
    try
        Taxodist._sleep[] = seconds -> push!(pauses, seconds)
        Taxodist._http_request[] = (method, url, headers, body; kwargs...) -> begin
            push!(captured, (method=method, url=url, headers=headers, body=body))
            source_response("{\"a\":1}")
        end
        res = Taxodist._source_request("ott", "POST", "http://example.org"; body=(a=1,))
        @test res == (status=200, text="{\"a\":1}")
        @test captured[1].method == "POST"
        @test captured[1].body == "{\"a\":1}"
        @test ("User-Agent" => Taxodist._user_agent) in captured[1].headers
        @test ("Content-Type" => "application/json") in captured[1].headers
        Taxodist._source_request("ott", "GET", "http://example.org")
        @test captured[2].body == ""
        @test length(captured[2].headers) == 1

        responses = Any[Base.IOError("offline", 0), source_response(""; status=429), source_response("ok")]
        Taxodist._http_request[] = (method, url, headers, body; kwargs...) -> begin
            next = popfirst!(responses)
            next isa Exception && throw(next)
            next
        end
        @test Taxodist._source_request("ncbi", "GET", "http://example.org") == (status=200, text="ok")
        @test pauses == [1, 2]

        Taxodist._http_request[] = (args...; kwargs...) -> source_response("gone"; status=404)
        empty!(pauses)
        @test Taxodist._source_request("ott", "GET", "http://example.org") == (status=404, text="gone")
        @test isempty(pauses)
        Taxodist._http_request[] = (args...; kwargs...) -> source_response("down"; status=503)
        @test Taxodist._source_request("ott", "GET", "http://example.org").status == 503
        empty!(pauses)
        Taxodist._http_request[] = (args...; kwargs...) -> source_response("slow down"; status=429)
        @test Taxodist._source_request("ncbi", "GET", "http://example.org").status == 429
        @test pauses == [1, 2]
        Taxodist._http_request[] = (args...; kwargs...) -> throw(Base.IOError("offline", 0))
        @test Taxodist._source_request("ott", "GET", "http://example.org") == (status=nothing, text=nothing)
        Taxodist._http_request[] = (args...; kwargs...) -> error("bug")
        @test_throws ErrorException Taxodist._source_request("ott", "GET", "http://example.org")
    finally
        Taxodist._http_request[] = original
        Taxodist._sleep[] = original_sleep
    end

    @test Taxodist._parse_json(nothing) === nothing
    @test Taxodist._parse_json("{not json") === nothing
    @test Taxodist._parse_json("{\"a\": 1}")[:a] == 1
    @test Taxodist._json_get([1, 2], :a) === nothing
end

@testset "NCBI merged identifiers and malformed responses are handled" begin
    clear_cache()
    try
        Taxodist._ncbi_parse_taxa("not xml")
        Taxodist._ncbi_parse_taxa("<Other><Taxon><TaxId>1</TaxId></Taxon></Other>")
        @test Taxodist._source_cache_get("src_ncbi_lin_1") === nothing
        xml = """
        <?xml version="1.0" ?>
        <!DOCTYPE TaxaSet PUBLIC "-//NLM//DTD Taxon, 14th January 2002//EN" "https://www.ncbi.nlm.nih.gov/entrez/query/DTD/taxon.dtd">
        <!-- recorded -->
        <TaxaSet><Taxon><TaxId>2</TaxId><ScientificName>B</ScientificName>
        <Rank>genus</Rank><AkaTaxIds><TaxId>3</TaxId></AkaTaxIds>
        <LineageEx><Taxon><TaxId>1</TaxId><ScientificName>A</ScientificName>
        <Rank>no rank</Rank></Taxon></LineageEx></Taxon></TaxaSet>
        """
        Taxodist._ncbi_parse_taxa(xml)
        merged = Taxodist._source_cache_get("src_ncbi_lin_3")
        @test merged.nodes == ["A", "B"]
        @test isequal(merged.ranks, [missing, "Genus"])
        @test isequal(Taxodist._source_cache_get("src_ncbi_lin_2"), merged)
    finally
        clear_cache()
    end
end

@testset "Minimal XML reader" begin
    parse_xml = Taxodist._parse_xml
    root = parse_xml("""<?xml version="1.0"?><a><!-- c --><?pi x?><b>x &amp; y &lt;&gt; &quot;&apos; &#65;&#x42; é</b><c/><d><![CDATA[<raw> é]]></d><e> <f>1</f>2</e></a>""")
    @test root.name == "a"
    @test Taxodist._xml_child_text(root, "b") == "x & y <> \"' AB é"
    @test Taxodist._xml_child_text(root, "c") == ""
    @test Taxodist._xml_child_text(root, "d") == "<raw> é"
    @test Taxodist._xml_child_text(root, "e") == " 12"
    @test Taxodist._xml_child_text(root, "missing") == ""
    @test isempty(Taxodist._xml_children(root, "missing"))
    @test parse_xml("<é></é>").name == "é"
    @test parse_xml("<a ></a >").name == "a"
    for bad in ["", "text", "<a>", "<a></b>", "<a>&bogus;</a>", "<a></a>trailing", "<a><!-- x</a>",
                "<></>", "<a><![CDATA[x</a>"]
        @test_throws ArgumentError parse_xml(bad)
    end
end

@testset "Source records survive save_cache and load_cache" begin
    with_sources() do
        taxo_resolve(["Homo", "Mus"]; source="ott", progress=false)
        taxo_resolve(["Homo", "Mus"]; source="ncbi", progress=false)
        snapshot = Taxodist._cache_snapshot()
        mktempdir() do directory
            file = joinpath(directory, "cache.json")
            redirect_stdout(devnull) do
                save_cache(file)
                clear_cache()
                load_cache(file)
            end
            for (key, value) in snapshot
                startswith(key, "src_") || continue
                @test isequal(Taxodist._taxodist_cache[key], value)
            end
            # Restored records are used without any request.
            Taxodist._http_request[] = failing_http
            res = taxo_resolve(["Homo", "Mus"]; source="ott", progress=false)
            @test res.status == ["resolved", "resolved"]
            @test startswith(res.source, "Open Tree Taxonomy 3.7")
            ncbi = taxo_resolve(["Homo", "Mus"]; source="ncbi", progress=false)
            @test ncbi.id == ["ncbi:9605", "ncbi:10088"]
        end
    end
    value = Taxodist._source_cache_value
    invalid = r"malformed source entry"
    @test_throws invalid value("src_ott_lin_1", JSON3.read("""{"nodes": ["a"], "ranks": []}"""))
    @test_throws invalid value("src_ncbi_lin_1", JSON3.read("""[1]"""))
    @test_throws invalid value("src_ncbi_search_x", JSON3.read("""[1]"""))
    @test_throws invalid value("src_ott_match_x", JSON3.read("""[{"id": "ott:1", "name": "x"}]"""))
    @test_throws invalid value("src_ott_about", JSON3.read("""{"name": "x"}"""))
    @test value("src_other", 1) == 1
    @test isempty(value("src_ncbi_search_x", JSON3.read("[]")))
end

@testset "compare_sources compares distances across sources" begin
    with_sources() do
        cmp = compare_sources(
            ["Homo", "Pan", "Mus", "Quercus", "Fakeosaurus", "Homo"];
            sources=["ott", "ncbi"], permutations=99, rng=Xoshiro(42), progress=false,
        )
        @test cmp isa TaxodistSourceComparison
        @test cmp.shared == ["Homo", "Pan", "Mus", "Quercus"]
        @test nrow(cmp.taxa) == 6
        @test all(column -> column in names(cmp.taxa), ["status_ott", "id_ncbi", "depth_ncbi"])
        @test ismissing(cmp.taxa.id_ott[5]) && ismissing(cmp.taxa.depth_ncbi[5])
        @test cmp.taxa.status_ncbi == ["resolved", "resolved", "resolved", "resolved", "unresolved", "resolved"]
        @test cmp.agreement.n_taxa == [4]
        @test cmp.agreement.spearman[1] > 0.9
        @test 0 < cmp.agreement.p_value[1] <= 1
        @test size(cmp.matrices["ott"]) == (4, 4)
        @test cmp.matrices["ncbi"].taxa == cmp.shared
        @test cmp.metric == "reciprocal" && cmp.permutations == 99
        @test cmp.resolutions["ncbi"].source == "NCBI Taxonomy"
        shown = sprint(show, MIME"text/plain"(), cmp)
        @test occursin("Taxa resolved in every source: 4 of 6", shown)
        @test occursin("Sources: ott, ncbi", shown)
        @test occursin("Median lineage depth: ott 33, ncbi 30", shown)
        @test occursin("spearman", shown)

        no_test = compare_sources(
            ["Homo", "Pan", "Mus"]; sources=["ott", "ncbi"], permutations=0,
            within=(ott="Metazoa",), progress=false,
        )
        @test ismissing(no_test.agreement.p_value[1])
        shared_within = compare_sources(
            ["Homo", "Pan", "Mus"]; sources=["ott", "ncbi"], permutations=0.5,
            within="Metazoa", progress=false,
        )
        @test shared_within.shared == ["Homo", "Pan", "Mus"]
        @test shared_within.permutations == 0
        by_dict = compare_sources(
            ["Homo", "Pan", "Quercus"]; sources=["ott", "ncbi"], permutations=0,
            within=Dict("ncbi" => "Metazoa"), progress=false, metric="path",
        )
        @test by_dict.shared == ["Homo", "Pan"]
        @test ismissing(by_dict.agreement.spearman[1])
        by_pairs = compare_sources(
            ["Homo", "Pan", "Quercus"]; sources=["ott", "ncbi"], permutations=0,
            within=["ott" => "Metazoa"], progress=false,
        )
        @test by_pairs.taxa.status_ott[3] == "unresolved"
        @test by_pairs.taxa.status_ncbi[3] == "resolved"
        single = compare_sources(
            ["Homo", "Pan", "Quercus"]; sources=["ott", "ncbi"], permutations=0,
            within=["Metazoa"], progress=false,
        )
        @test single.shared == ["Homo", "Pan"]
        none = compare_sources(
            ["Fakeosaurus", "Pan", "Quercus"]; sources=["ott", "ncbi"], permutations=0,
            within="Fungi", progress=false,
        )
        @test isempty(none.shared)
        @test occursin("Median lineage depth: ott NA, ncbi NA", sprint(show, MIME"text/plain"(), none))
    end
end

@testset "compare_sources validates its arguments" begin
    @test_throws r"at least two" compare_sources(["a", "b", "c"]; sources=["ott"])
    @test_throws r"at least two" compare_sources(["a", "b", "c"]; sources=["ott", "ott"])
    @test_throws r"at least two" compare_sources(["a", "b", "c"]; sources="ott")
    @test_throws r"Unknown" compare_sources(["a", "b", "c"]; sources=["ott", "gbif"])
    @test_throws r"at least three" compare_sources(["a", "b"]; sources=["ott", "ncbi"])
    @test_throws r"non-negative" compare_sources(["a", "b", "c"]; sources=["ott", "ncbi"], permutations=-1)
    @test_throws r"non-negative" compare_sources(["a", "b", "c"]; sources=["ott", "ncbi"], permutations=true)
    @test_throws r"non-negative" compare_sources(["a", "b", "c"]; sources=["ott", "ncbi"], permutations=Inf)
    @test_throws r"named by source" compare_sources(
        ["Homo", "Pan", "Mus"]; sources=["ott", "ncbi"], within=["a", "b"], progress=false,
    )
    @test_throws r"named by source" compare_sources(
        ["Homo", "Pan", "Mus"]; sources=["ott", "ncbi"], within=1, progress=false,
    )
end

@testset "The Mantel statistic handles degenerate inputs" begin
    mantel = Taxodist._mantel
    labels = ["x", "y", "z"]
    a = TaxonomicDistanceMatrix([0 1 2; 1 0 3; 2 3 0], labels)
    b = TaxonomicDistanceMatrix([0 Inf 2; Inf 0 3; 2 3 0], labels)
    @test isequal(mantel(a, b, 10, Xoshiro(1)), (r=missing, p=missing))
    flat = TaxonomicDistanceMatrix([0 1 1; 1 0 1; 1 1 0], labels)
    @test ismissing(mantel(a, flat, 10, Xoshiro(1)).r)
    same = mantel(a, a, 9, Xoshiro(1))
    @test same.r ≈ 1
    @test 0 < same.p <= 1
    @test Taxodist._average_ranks([3.0, 1.0, 3.0, 2.0]) == [3.5, 1.0, 3.5, 2.0]
end

@testset "Verbose resolution reports each request" begin
    with_sources() do
        output = capture_stdout(() -> taxo_resolve(["Homo", "Pan"]; source="ott", progress=false, verbose=true))
        @test occursin("Matching 2 names in the Open Tree Taxonomy", output)
        output = capture_stdout(() -> taxo_resolve("Mus"; source="ott", progress=false, verbose=true))
        @test occursin("Matching 1 name in the Open Tree Taxonomy", output)
        @test occursin("Retrieving lineage of ott:1068778", output)
        output = capture_stdout(() -> taxo_resolve("Quercus"; source="ncbi", progress=false, verbose=true))
        @test occursin("Searching NCBI Taxonomy for Quercus", output)
        @test occursin("Retrieving 1 lineage from NCBI Taxonomy", output)
        output = capture_stdout(() -> taxo_resolve("Homo"; source="ncbi", progress=false, verbose=true))
        @test occursin("Retrieving 1 lineage from NCBI Taxonomy", output)
        output = capture_stdout(() -> Taxodist._ncbi_fetch(["ncbi:9596", "ncbi:9598"]; verbose=true))
        @test occursin("Retrieving 2 lineages from NCBI Taxonomy", output)
    end
end

@testset "Missing metadata and failed lineage batches are tolerated" begin
    function partial(method, url, headers, body; kwargs...)
        endswith(url, "/taxonomy/about") && return source_response("{\"version\": \"9.9\", \"source\": \"ott9.9\"}")
        occursin("efetch.fcgi", url) && return source_response(""; status=503)
        return fixture_http(method, url, headers, body)
    end
    with_sources(partial) do
        about = Taxodist._ott_about()
        @test about.name == "Open Tree Taxonomy 9.9 (ott9.9)"
        @test about.url == "https://tree.opentreeoflife.org/about/taxonomy-version"
        # The name was found but its lineage could not be retrieved.
        res = taxo_resolve("Quercus"; source="ncbi", progress=false)
        @test res.status == ["retrieval_error"]
        @test res.n_candidates == [0]
    end
    function odd(method, url, headers, body; kwargs...)
        endswith(url, "/taxonomy/about") && return source_response("{}")
        endswith(url, "/tnrs/match_names") && return source_response("[]")
        occursin("esearch.fcgi", url) && return source_response("{\"esearchresult\": {}}")
        return fixture_http(method, url, headers, body)
    end
    with_sources(odd) do
        @test Taxodist._ott_about().name == "Open Tree Taxonomy  ()"
        res = taxo_resolve(["Homo"]; source="ott", progress=false)
        @test res.status == ["unresolved"]
        @test taxo_resolve(["Homo"]; source="ncbi", progress=false).status == ["unresolved"]
    end
    # A fungus homonym whose lineage is missing is dropped from the candidates.
    function no_fungus(method, url, headers, body; kwargs...)
        if occursin("efetch.fcgi", url)
            url = replace(url, "%2C2081351" => "", "2081351%2C" => "")
        end
        return fixture_http(method, url, headers, body)
    end
    with_sources(no_fungus) do
        res = taxo_resolve(["Drosophila"]; source="ncbi", progress=false)
        @test res.status == ["resolved"]
        @test res.id == ["ncbi:7215"]
    end
    # A synonym match comes after a direct match.
    function synonyms(method, url, headers, body; kwargs...)
        if endswith(url, "/tnrs/match_names")
            return source_response("""{"results": [{"name": "Homo", "matches": [
                {"is_synonym": true, "taxon": {"ott_id": 417957, "name": "Homo"}},
                {"is_synonym": false, "taxon": {"ott_id": 770309, "name": "Homo"}}]}]}""")
        end
        return fixture_http(method, url, headers, body)
    end
    with_sources(synonyms) do
        res = taxo_resolve(["Homo"]; source="ott", progress=false, ambiguity="first")
        @test res.status == ["ambiguous"]
        @test res.id == ["ott:770309"]
        @test res.candidates[1].id == ["ott:770309", "ott:417957"]
    end
end
