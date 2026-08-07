@testset "Lineage comparison and shared clades" begin
    clear_cache()
    seed_lineage("Alpha", "1", ["Biota", "Animalia", "CladeA", "Alpha"])
    seed_lineage("Beta", "2", ["Biota", "Animalia", "CladeA", "Beta"])
    seed_lineage("Gamma", "3", ["Biota", "Animalia", "CladeB", "Gamma"])
    Taxodist._taxodist_cache["id_Missing"] = nothing

    comparison = compare_lineages("Alpha", "Beta")
    @test comparison.lineage_a[end] == "Alpha"
    @test comparison.lineage_b[end] == "Beta"
    @test comparison.mrca_depth == 3

    @test shared_clades("Alpha", "Beta") == ["Biota", "Animalia", "CladeA"]
    @test shared_clades("Alpha", "Gamma") == ["Biota", "Animalia"]
    @test shared_clades("Alpha", "Missing") === nothing
    @test compare_lineages("Alpha", "Missing") === nothing

    clear_cache()
end

@testset "Disconnected lineage utilities" begin
    clear_cache()
    seed_lineage("Animal", "1", ["Biota", "Animalia", "Animal"])
    seed_lineage("Mineral", "2", ["Mineralia", "Mineral"])

    @test shared_clades("Animal", "Mineral") == String[]
    @test taxo_path("Animal", "Mineral") === nothing
    disconnected = compare_lineages("Animal", "Mineral")
    @test disconnected.mrca_depth == 0

    clear_cache()
end

@testset "Clade membership and filtering" begin
    clear_cache()
    seed_lineage("Alpha", "1", ["Biota", "Animalia", "Clade (example)", "Species+", "Alpha"])
    seed_lineage("Beta", "2", ["Biota", "Animalia", "Other", "Beta"])
    Taxodist._taxodist_cache["id_Missing"] = nothing

    @test is_member("Alpha", "animalia") === true
    @test is_member("1", "  CLADE (EXAMPLE)  ") === true
    @test is_member("1", "Species+") === true
    @test is_member("1", "Clade (") === false
    @test is_member("Alpha", "Anim") === false
    @test is_member("Missing", "Animalia") === nothing
    @test filter_clade(["Alpha", "Beta", "Missing"], "Animalia") == ["Alpha", "Beta"]
    @test filter_clade(["1", "2", "Missing"], "Clade (example)") == ["1"]

    clear_cache()
end

@testset "Taxonomic paths" begin
    clear_cache()
    seed_lineage("Alpha", "1", ["Biota", "Animalia", "CladeA", "Alpha"])
    seed_lineage("Gamma", "2", ["Biota", "Animalia", "CladeB", "Gamma"])
    seed_lineage("CladeA", "3", ["Biota", "Animalia", "CladeA"])

    path = taxo_path("Alpha", "Gamma")
    @test path.taxon_a == "Alpha"
    @test path.taxon_b == "Gamma"
    @test path.data.node == ["Alpha", "CladeA", "Animalia", "CladeB", "Gamma"]
    @test path.data.depth == [4, 3, 2, 3, 4]
    @test path.data.direction == ["a", "a", "mrca", "b", "b"]
    @test count(==("mrca"), path.data.direction) == 1

    identical_path = taxo_path("Alpha", "Alpha")
    @test size(identical_path.data, 1) == 1
    @test identical_path.data.node == ["Alpha"]
    @test identical_path.data.direction == ["mrca"]

    ancestor_path = taxo_path("CladeA", "Alpha")
    @test ancestor_path.data.node == ["CladeA", "Alpha"]
    @test ancestor_path.data.direction == ["mrca", "b"]

    reverse_path = taxo_path("Alpha", "CladeA")
    @test reverse_path.data.node == ["Alpha", "CladeA"]
    @test reverse_path.data.direction == ["a", "mrca"]

    Taxodist._taxodist_cache["id_Missing"] = nothing
    @test taxo_path("Missing", "Alpha") === nothing
    @test taxo_path("Alpha", "Missing") === nothing

    clear_cache()
end

@testset "Analysis summaries and plots" begin
    matrix = example_analysis_matrix()
    cluster = taxo_cluster(matrix)
    ordination = taxo_ordinate(matrix; k=2)

    @test plot_taxodist_cluster(cluster) !== nothing
    @test plot_taxodist_ord(ordination) !== nothing

    summary = summary_taxodist_ord(ordination)
    @test names(summary) == ["Axis", "Eigenvalue", "Variance_Pct", "Cumulative_Pct"]
    @test summary.Axis == ["PC1", "PC2"]
    @test summary.Cumulative_Pct[end] ≈ 100.0

    @test plot_taxodist_cluster((hclust=nothing, dist=matrix)) === nothing
    @test plot_taxodist_ord((points=nothing, GOF=nothing)) === nothing
    @test summary_taxodist_ord((points=nothing, eig=nothing)) === nothing
    @test summary_taxodist_ord((points=Taxodist.DataFrame(taxon=["A"]), eig=nothing)) === nothing

    no_axes = (points=Taxodist.DataFrame(taxon=["A", "B"]), eig=[0.0, 0.0], GOF=nothing)
    @test size(summary_taxodist_ord(no_axes)) == (0, 4)
    @test plot_taxodist_ord(no_axes) === nothing

    one_axis = (
        points=Taxodist.DataFrame(taxon=["A", "B"], PC1=[-0.5, 0.5]),
        eig=[1.0, 0.0],
        GOF=nothing,
    )
    @test plot_taxodist_ord(one_axis) !== nothing

    zero_variance = (
        points=Taxodist.DataFrame(taxon=["A", "B"], PC1=[0.0, 0.0]),
        eig=[0.0, 0.0],
        GOF=[NaN, NaN],
    )
    zero_summary = summary_taxodist_ord(zero_variance)
    @test isempty(zero_summary.Axis)
end

@testset "Taxonomic heatmap" begin
    matrix = example_analysis_matrix()
    @test taxo_heatmap(matrix; display_plot=false) === matrix

    singleton = TaxonomicDistanceMatrix(zeros(1, 1), ["A"])
    @test taxo_heatmap(singleton; display_plot=false) === singleton

    missing_matrix = TaxonomicDistanceMatrix([0.0 NaN; NaN 0.0], ["A", "B"])
    @test taxo_heatmap(missing_matrix; display_plot=false) === missing_matrix

    infinite_matrix = TaxonomicDistanceMatrix([0.0 Inf; Inf 0.0], ["A", "B"])
    @test taxo_heatmap(infinite_matrix; display_plot=false) === infinite_matrix
end