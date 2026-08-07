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
    end

    Taxodist.clear_cache()
end