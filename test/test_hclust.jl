# taxo_cluster() uses a port of the Fortran routines HCLUST and HCASS2 of R's
# stats::hclust(). hclust_r.json holds stats::hclust() results (R 4.6.2) for
# random dissimilarities - taxonomic distances with heavy ties, small
# integers, rounded and continuous values - under all eight methods; merges,
# heights and leaf order must be identical.

@testset "hclust port reproduces R's stats::hclust() exactly" begin
    cases = Taxodist.JSON3.read(read(joinpath(@__DIR__, "fixtures", "hclust_r.json"), String))
    as_vector(x) = x isa AbstractVector ? collect(x) : [x]
    @test length(cases) > 50
    n_checked = 0
    for case in cases
        n = case.n
        dist = Float64.(as_vector(case.dist))
        values = zeros(n, n)
        position = 0
        for j in 1:(n - 1), i in (j + 1):n   # R's dist storage, column by column
            position += 1
            values[i, j] = values[j, i] = dist[position]
        end
        for result in case.results
            hc = Taxodist._hclust(values, String(result.method))
            merges = reduce(vcat, [permutedims(Int.(collect(row))) for row in result.merge])
            @test hc.merges == merges
            @test hc.heights == Float64.(as_vector(result.height))
            @test hc.order == Int.(as_vector(result.order))
            n_checked += 1
        end
    end
    @test n_checked == 8 * length(cases)
end

@testset "taxo_cluster supports all of R's methods" begin
    values = [0 2 6 10 9; 2 0 5 9 8; 6 5 0 4 5; 10 9 4 0 3; 9 8 5 3 0]
    matrix = TaxonomicDistanceMatrix(values, ["a", "b", "c", "d", "e"])
    for (method, linkage) in (
        "ward.D" => :ward_presquared, "single" => :single, "complete" => :complete,
        "average" => :average, "mcquitty" => :mcquitty, "median" => :median,
        "centroid" => :centroid, "ward.D2" => :ward, "ward" => :ward,
        "ward_presquared" => :ward_presquared, "AVERAGE" => :average,
    )
        result = taxo_cluster(matrix; method=method)
        @test result.method === linkage
        @test result.hclust isa Taxodist.Clustering.Hclust
        @test result.hclust.linkage === linkage
        @test sort(result.hclust.order) == 1:5
        @test size(result.hclust.merges) == (4, 2)
    end
    @test taxo_cluster(matrix; method=:single).method === :single
    @test_throws ArgumentError taxo_cluster(matrix; method="nonsense")
    @test_throws ArgumentError taxo_cluster(matrix; method=1)
    @test_throws ArgumentError Taxodist._hclust(zeros(1, 1), "average")

    # Single linkage on this matrix, checked against R by hand:
    # hclust(as.dist(values), "single") gives merge
    # (-1,-2), (-4,-5), (-3,2), (1,3); heights 2, 3, 4, 5; order 1 2 3 4 5.
    single = taxo_cluster(matrix; method="single").hclust
    @test single.merges == [-1 -2; -4 -5; -3 2; 1 3]
    @test single.heights == [2.0, 3.0, 4.0, 5.0]
    @test single.order == [1, 2, 3, 4, 5]
end
