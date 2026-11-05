using Test
using Taxodist
using Plots  # activates the plotting extension

# Keep retry logic active but skip the real backoff pauses and the request
# throttle in tests.
Taxodist._sleep[] = seconds -> nothing
Taxodist.MIN_INTERVAL[] = 0.0

@testset "Taxodist" begin
    include("test_distance.jl")
    include("test_fetch.jl")
    include("test_distance_api.jl")
    include("test_utils.jl")
    include("test_data.jl")
    include("test_resolution.jl")
    include("test_matrix_tree.jl")
    include("test_parser.jl")
    include("test_regressions.jl")
    include("test_metrics.jl")
end
