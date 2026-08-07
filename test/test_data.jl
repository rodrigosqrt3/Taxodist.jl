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

@testset "Taxobase conversion helpers" begin
    object_matrix = Taxodist.JSON3.read("{\"labels\":[\"A\",\"B\"],\"data\":[[0,0.5],[0.5,0]]}")
    matrix = Taxodist._matrix_from_json(object_matrix, String[])
    @test matrix.taxa == ["A", "B"]
    @test matrix[1, 2] == 0.5

    fallback = Taxodist._matrix_from_json([[0, 1], [1, 0]], ["A", "B"])
    @test fallback.taxa == ["A", "B"]
    @test_throws ArgumentError Taxodist._matrix_from_json([[0, 1]], ["A", "B"])
    @test_throws ArgumentError Taxodist._matrix_from_json([[0], [1]], ["A", "B"])

    empty_table = Taxodist._table_from_records(Any[], (:taxon, :distance))
    @test names(empty_table) == ["taxon", "distance"]
    @test size(empty_table, 1) == 0

    records = Taxodist.JSON3.read("[{\"taxon\":\"A\"},{}]")
    table = Taxodist._table_from_records(records, (:taxon, :distance))
    @test table.taxon == Any["A", nothing]
    @test table.distance == Any[nothing, nothing]

    vector_coverage = Taxodist._coverage_table(["A", "B"], [true, false])
    @test vector_coverage.taxon == ["A", "B"]
    @test vector_coverage.covered == [true, false]
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