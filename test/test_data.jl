@testset "Packaged taxobase structure" begin
    taxobase = load_taxobase()

    @test propertynames(taxobase) == (
        :taxa,
        :found_taxa,
        :coverage,
        :matrix,
        :pairwise,
        :lineage_homo,
        :lineage_tyrannosaurus,
        :closest,
        :filter,
        :search,
        :statistical_taxa,
        :statistical_matrix,
        :metadata,
    )

    @test length(taxobase.taxa) == 52
    @test length(taxobase.found_taxa) == 52
    @test size(taxobase.coverage) == (52, 2)
    @test all(taxobase.coverage.covered)

    @test size(taxobase.matrix) == (52, 52)
    @test taxobase.matrix.taxa == taxobase.found_taxa
    @test all(taxobase.matrix[i, i] == 0.0 for i in 1:52)
    @test all(isfinite, taxobase.matrix.values)
    @test taxobase.matrix.values == transpose(taxobase.matrix.values)

    @test length(taxobase.statistical_taxa) == 15
    @test size(taxobase.statistical_matrix) == (15, 15)
    @test taxobase.statistical_matrix.taxa == taxobase.statistical_taxa
    @test all(taxobase.statistical_matrix[i, i] == 0.0 for i in 1:15)
    @test all(isfinite, taxobase.statistical_matrix.values)

    @test taxobase.metadata.package_version == "0.6.0"
    @test taxobase.metadata.generated_on == "2026-08-06"
    @test taxobase.metadata.source == "The Taxonomicon"
end

@testset "Packaged taxobase examples" begin
    taxobase = load_taxobase()

    @test taxobase.lineage_homo[end] == "Homo"
    @test taxobase.lineage_tyrannosaurus[end] == "Tyrannosaurus"

    calculated = Taxodist._compute_distance(
        taxobase.lineage_tyrannosaurus,
        taxobase.lineage_homo;
        taxon_a="Tyrannosaurus",
        taxon_b="Homo",
    )
    @test calculated.distance ≈ taxobase.pairwise.distance
    @test calculated.mrca == taxobase.pairwise.mrca
    @test calculated.mrca_depth == taxobase.pairwise.mrca_depth

    @test names(taxobase.closest) == ["taxon", "distance"]
    @test size(taxobase.closest, 1) == 6
    @test taxobase.closest.taxon[1] == "Velociraptor"

    @test length(taxobase.filter) == 11
    @test "Tyrannosaurus" in taxobase.filter

    @test names(taxobase.search) == ["id", "name"]
    @test size(taxobase.search, 1) == 4
    @test all(occursin("Bacteria", name) for name in taxobase.search.name)
end