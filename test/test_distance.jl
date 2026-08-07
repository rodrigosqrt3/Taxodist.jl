@testset "Core hierarchy distance" begin
    lineage_a = [
        "Biota", "Animalia", "Chordata", "Dinosauria", "Theropoda",
        "Tyrannosauridae", "Tyrannosaurus",
    ]
    lineage_b = [
        "Biota", "Animalia", "Chordata", "Dinosauria", "Theropoda",
        "Dromaeosauridae", "Velociraptor",
    ]

    result = Taxodist._compute_distance(
        lineage_a,
        lineage_b;
        taxon_a="Tyrannosaurus",
        taxon_b="Velociraptor",
    )

    @test result.distance == 1 / 5
    @test result.mrca == "Theropoda"
    @test result.mrca_depth == 5
    @test result.depth_a == 7
    @test result.depth_b == 7
end

@testset "Identity and ancestor-descendant pairs" begin
    ancestor = ["Biota", "Animalia", "Dinosauria"]
    descendant = ["Biota", "Animalia", "Dinosauria", "Theropoda", "A"]

    @test Taxodist._compute_distance(ancestor, ancestor).distance == 0.0
    @test Taxodist._compute_distance(ancestor, descendant).distance == 1 / 3
end

@testset "Continuous common prefix" begin
    lineage_a = ["Biota", "Animalia", "Alpha", "Repeated", "A"]
    lineage_b = ["Biota", "Animalia", "Beta", "Repeated", "B"]
    result = Taxodist._compute_distance(lineage_a, lineage_b)

    @test result.mrca == "Animalia"
    @test result.mrca_depth == 2
    @test result.distance == 1 / 2
end

@testset "Disconnected hierarchies" begin
    result = Taxodist._compute_distance(["Biota"], ["Natura"])

    @test isinf(result.distance)
    @test isnothing(result.mrca)
    @test result.mrca_depth == 0
end

@testset "Symmetry and ultrametricity" begin
    lineage_a = ["Biota", "Animalia", "Dinosauria", "Theropoda", "A"]
    lineage_b = ["Biota", "Animalia", "Dinosauria", "Theropoda", "B"]
    lineage_c = ["Biota", "Animalia", "Dinosauria", "Ornithischia", "C"]
    distance(a, b) = Taxodist._compute_distance(a, b).distance

    @test distance(lineage_a, lineage_b) == distance(lineage_b, lineage_a)
    @test distance(lineage_a, lineage_c) <=
          max(distance(lineage_a, lineage_b), distance(lineage_b, lineage_c))
end