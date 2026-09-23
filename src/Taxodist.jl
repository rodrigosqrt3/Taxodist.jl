module Taxodist

using Cascadia
using Clustering
using DataFrames
using Dates
using Gumbo
using HTTP
using JSON3
using LinearAlgebra
using Plots

const TAXODIST_VERSION = v"0.8.0"

export clear_cache,
       save_cache,
       load_cache,
       cache_info,
       get_taxonomicon_id,
       get_lineage_by_id,
       get_lineage,
       taxo_search,
       TaxodistResolution,
       taxo_resolve,
       taxo_from_lineages,
       summary_counts,
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
       load_taxobase,
       TaxodistBundle,
       taxo_bundle,
       validate_taxodist_bundle,
       write_taxodist_bundle,
       read_taxodist_bundle

include("fetch.jl")
include("distance.jl")
include("bundle.jl")
include("utils.jl")
include("data.jl")

end