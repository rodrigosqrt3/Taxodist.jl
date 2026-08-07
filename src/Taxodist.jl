module Taxodist

using Cascadia
using Clustering
using DataFrames
using Gumbo
using HTTP
using JSON3
using LinearAlgebra
using Plots

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
       taxo_path,
       taxo_cluster,
       taxo_ordinate,
       taxo_heatmap,
       plot_taxodist_cluster,
       plot_taxodist_ord,
       summary_taxodist_ord,
       load_taxobase

include("fetch.jl")
include("distance.jl")
include("utils.jl")
include("data.jl")

end