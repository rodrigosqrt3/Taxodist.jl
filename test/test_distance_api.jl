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
    @test matrix.taxa == ["Alpha", "Beta", "Gamma"]
    @test [matrix[i, i] for i in 1:3] == zeros(3)
    @test matrix["Alpha", "Beta"] == 1 / 3
    @test matrix[1, 3] == 1 / 2
    @test Matrix(matrix) == transpose(Matrix(matrix))
    @test_throws KeyError matrix["Missing", "Beta"]
    @test size(distance_matrix(String[]; progress=false)) == (0, 0)
    @test_throws ArgumentError distance_matrix(["Alpha", "Alpha"]; progress=false)

    clear_cache()
end

@testset "Unresolved taxa remain explicit" begin
    clear_cache()
    seed_lineage("Alpha", "1", ["Biota", "Animalia", "Alpha"])
    Taxodist._taxodist_cache["id_Missing"] = nothing

    matrix = distance_matrix(["Alpha", "Missing"]; progress=false)
    @test isnan(matrix["Alpha", "Missing"])
    @test matrix["Missing", "Missing"] == 0.0

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

    focal = focal_distances("Alpha", ["Gamma", "Alpha", "Missing", "Beta"])
    @test focal.focal == "Alpha"
    @test focal.data.taxon == ["Alpha", "Beta", "Gamma", "Missing"]
    @test focal.data.distance[1] == 0.0
    @test focal.data.mrca[1] == "Alpha"
    @test focal.data.mrca_depth[1] == 4
    @test isnan(focal.data.distance[4])

    @test lineage_depth("Alpha") == 4
    @test lineage_depth("Missing") === nothing

    coverage = check_coverage(["Alpha", "Missing"])
    @test coverage.taxon == ["Alpha", "Missing"]
    @test coverage.covered == [true, false]

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

    singleton = TaxonomicDistanceMatrix(zeros(1, 1), ["A"])
    @test taxo_ordinate(singleton).points === nothing

    missing_matrix = TaxonomicDistanceMatrix([0.0 NaN; NaN 0.0], ["A", "B"])
    @test taxo_ordinate(missing_matrix).points === nothing

    infinite_matrix = TaxonomicDistanceMatrix([0.0 Inf; Inf 0.0], ["A", "B"])
    @test taxo_ordinate(infinite_matrix).points === nothing
end