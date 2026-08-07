module Taxodist

using Cascadia
using DataFrames
using Gumbo
using HTTP
using JSON3

export clear_cache,
       save_cache,
       load_cache,
       cache_info,
       get_taxonomicon_id,
       get_lineage_by_id,
       get_lineage,
       taxo_search,
       taxo_distance,
       mrca,
       TaxonomicDistanceMatrix,
       distance_matrix,
       closest_relative,
       focal_distances,
       lineage_depth,
       check_coverage,
       compare_lineages,
       shared_clades,
       is_member,
       filter_clade,
       taxo_path

include("fetch.jl")
include("distance.jl")
include("utilities.jl")

end