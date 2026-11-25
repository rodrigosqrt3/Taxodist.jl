"""
    Taxodist

Taxonomic hierarchy distances and lineage analysis. Distances are computed
from the depth of the most recent common ancestor of ordered lineages, which
come from one of three taxonomic sources, chosen with
`taxo_resolve(taxa; source=...)` or `taxo_bundle(taxa; source=...)`:

- `"taxonomicon"` (default): The Taxonomicon; numeric identifiers (`"50841"`).
- `"ott"`: the Open Tree of Life reference taxonomy; identifiers `"ott:770309"`.
- `"ncbi"`: the NCBI Taxonomy; identifiers `"ncbi:9606"`. Requests are spaced
  by at least 0.4 s, or 0.12 s with an API key in the `NCBI_API_KEY`
  environment variable or `Taxodist.NCBI_API_KEY[]`.

Name-based functions (`get_lineage`, `taxo_distance`, `mrca`, ...) query The
Taxonomicon; resolutions from any source, and curated lineages from
`taxo_from_lineages`, can be passed to the matrix, tree, bundle, and community
functions. `taxo_sources()` lists the sources with their licences and
citations; please cite the source used as well as the package.

# Exported functions and types

- Retrieval and resolution: `get_lineage`, `get_lineage_by_id`,
  `get_taxonomicon_id`, `lineage_ranks`, `taxo_search`, `taxo_resolve`,
  `TaxodistResolution`, `summary_counts`, `taxo_from_lineages`,
  `taxonomicon_available`.
- Sources: `taxo_sources`, `compare_sources`, `TaxodistSourceComparison`.
- Distances: `taxo_distance`, `mrca`, `distance_matrix`,
  `TaxonomicDistanceMatrix`, `cross_distances`, `TaxonomicCrossDistances`,
  `closest_relative`, `focal_distances`, `lineage_depth`, `check_coverage`.
- Lineage utilities: `compare_lineages`, `shared_clades`, `is_member`,
  `filter_clade`, `taxo_path`.
- Trees, clustering, and ordination: `taxo_tree`, `TaxodistTree`, `newick`,
  `cophenetic`, `taxo_cluster`, `taxo_ordinate`, `summary_taxodist_ord`, and,
  with Plots.jl loaded, `taxo_heatmap`, `plot_taxodist_cluster`,
  `plot_taxodist_ord`.
- Community analyses: `taxo_distinctness`, `TaxodistDistinctness`,
  `taxo_distinctness_test`, `taxo_dissimilarity`, `taxo_ses`, `TaxodistSES`.
- Bundles: `taxo_bundle`, `TaxodistBundle`, `validate_taxodist_bundle`,
  `write_taxodist_bundle`, `read_taxodist_bundle`, `compare_bundles`,
  `TaxodistBundleComparison`.
- Cache and reference data: `clear_cache`, `save_cache`, `load_cache`,
  `cache_info`, `load_taxobase`.

# References

- Brands, S.J. (ed.) 1989-present. The Taxonomicon. Universal Taxonomic
  Services, Zwaag, The Netherlands. <http://taxonomicon.taxonomy.nl>
- Rees, J.A. & Cranston, K. (2017). Automated assembly of a reference taxonomy
  for phylogenetic data synthesis. Biodiversity Data Journal 5: e12581.
  <https://doi.org/10.3897/BDJ.5.e12581>
- Schoch, C.L. et al. (2020). NCBI Taxonomy: a comprehensive update on
  curation, resources and tools. Database 2020: baaa062.
  <https://doi.org/10.1093/database/baaa062>
"""
module Taxodist

using Cascadia
using Clustering
using DataFrames
using Dates
using Gumbo
using HTTP
using JSON3
using LinearAlgebra
using Random
using Random: AbstractRNG, randperm
using Statistics: cor, mean, median, quantile, std

# The version is read from Project.toml at precompilation; declaring the file
# as a dependency recompiles the package when it changes, so the version (and
# the user agent built from it) is never stale.
include_dependency(joinpath(dirname(@__DIR__), "Project.toml"))
const TAXODIST_VERSION = pkgversion(@__MODULE__)

export clear_cache,
       save_cache,
       load_cache,
       cache_info,
       get_taxonomicon_id,
       get_lineage_by_id,
       get_lineage,
       lineage_ranks,
       taxo_search,
       TaxodistResolution,
       taxo_resolve,
       taxo_from_lineages,
       taxonomicon_available,
       summary_counts,
       taxo_distance,
       mrca,
       TaxonomicDistanceMatrix,
       distance_matrix,
       TaxonomicCrossDistances,
       cross_distances,
       TaxodistTree,
       taxo_tree,
       newick,
       cophenetic,
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
       read_taxodist_bundle,
       TaxodistBundleComparison,
       compare_bundles,
       TaxodistDistinctness,
       taxo_distinctness,
       taxo_distinctness_test,
       taxo_sources,
       TaxodistSourceComparison,
       compare_sources,
       taxo_dissimilarity,
       TaxodistSES,
       taxo_ses

include("fetch.jl")
include("distance.jl")
include("sources.jl")
include("distinctness.jl")
include("community.jl")
include("bundle.jl")
include("tree.jl")
include("utils.jl")
include("data.jl")

end