# taxodist <picture><source media="(prefers-color-scheme: dark)" srcset="images/taxodist_dark.png"><source media="(prefers-color-scheme: light)" srcset="images/taxodist_sepia.png"><img alt="taxodist logo" src="images/taxodist_sepia.png" align="right" height="200"></picture>

[![Julia version](https://img.shields.io/badge/dynamic/toml?url=https%3A%2F%2Fraw.githubusercontent.com%2Frodrigosqrt3%2FTaxodist.jl%2Fmain%2FProject.toml&query=%24.version&prefix=v&label=julia&color=2f855a)](https://platform.juliahub.com/ui/Packages/General/Taxodist/) &nbsp; [![Julia Tests](https://github.com/rodrigosqrt3/Taxodist.jl/actions/workflows/julia.yml/badge.svg)](https://github.com/rodrigosqrt3/Taxodist.jl/actions/workflows/julia.yml) &nbsp; [![codecov](https://codecov.io/gh/rodrigosqrt3/Taxodist.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rodrigosqrt3/Taxodist.jl)

**Taxonomic hierarchy distance and lineage computation for any taxon on Earth.**

`Taxodist.jl` retrieves full hierarchical lineages and computes an ultrametric
distance between any two taxa: a pair of dinosaurs, a dinosaur and a fungus,
two species of fly, or an oak tree and a human. Lineages come from one of three
taxonomic sources: [The Taxonomicon](http://taxonomicon.taxonomy.nl) (the
default), the [Open Tree Taxonomy](https://tree.opentreeoflife.org/about/taxonomy-version),
and the [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy), or from
lineages you supply yourself. Besides pairwise distances and matrices, the
package provides classification trees, clustering and ordination, Clarke and
Warwick's taxonomic distinctness, taxonomic dissimilarities between samples,
and standardized effect sizes of community relatedness.

The same package is available for [R](https://github.com/rodrigosqrt3/taxodist)
and [Python](https://github.com/rodrigosqrt3/taxodist-py); the three
implementations give identical results and exchange analyses through a common
JSON bundle format. Documentation for all three is at
[taxodist.github.io](https://taxodist.github.io).

## Installation

```julia
using Pkg
Pkg.add("Taxodist")
```

The development version can be installed from GitHub:

```julia
Pkg.add(url="https://github.com/rodrigosqrt3/Taxodist.jl")
```

## Basic usage

Name-based functions query The Taxonomicon:

```julia
using Taxodist

# Get a full lineage
get_lineage("Tyrannosaurus")

# Distance between two taxa
taxo_distance("Tyrannosaurus", "Velociraptor")

# Most recent common ancestor
mrca("Tyrannosaurus", "Triceratops")  # "Dinosauria"
mrca("Tyrannosaurus", "Homo")         # "Amniota"

# Pairwise distance matrix
theropods = ["Tyrannosaurus", "Velociraptor", "Spinosaurus", "Allosaurus"]
distance_matrix(theropods)

# Filter taxa by clade
taxa = ["Tyrannosaurus", "Triceratops", "Homo", "Quercus"]
filter_clade(taxa, "Dinosauria")

# Get the path between two taxa
taxo_path("Tyrannosaurus", "Velociraptor")

# Save and restore the lineage cache across sessions
save_cache()             # persistent per-user file; or save_cache("my_cache.json")
load_cache()

# Distances between two sets of taxa (rows × columns)
cross_distances(["Tyrannosaurus", "Triceratops"], ["Velociraptor", "Homo"])

# Classification tree whose internal nodes are named clades
tree = taxo_tree(theropods)
newick(tree)             # export for Phylo.jl, ape, ete3, ...
cophenetic(tree)         # equals distance_matrix(theropods)
```

## Choosing a taxonomic source

`taxo_resolve()` and `taxo_bundle()` take a `source` keyword:

| `source` | classification | identifiers | root of every lineage |
|---|---|---|---|
| `"taxonomicon"` (default) | The Taxonomicon | numeric (`"50841"`) | `"Biota"` |
| `"ott"` | Open Tree of Life reference taxonomy | `"ott:770309"` | `"life"` |
| `"ncbi"` | NCBI Taxonomy | `"ncbi:9606"` | `"cellular organisms"` |

The resolution records its source (with the version of the Open Tree
Taxonomy) and works with every function that accepts one: `distance_matrix()`,
`cross_distances()`, `taxo_tree()`, `taxo_bundle()`, and, through the matrix,
the distinctness and community functions. When `source` is set, the prefix may
be omitted from identifiers in the input. Clade names differ between sources:
animals are `"Animalia"` in The Taxonomicon and `"Metazoa"` in the other two,
which matters for `within`.

```julia
taxo_sources()     # the sources, their identifier prefixes, licences, and citations

resolution = taxo_resolve(["Homo", "Pan", "Gorilla"]; source="ott")
distance_matrix(resolution)
taxo_resolve(["Homo", "Pan", "Gorilla"]; source="ncbi", within="Metazoa")
taxo_bundle(["Homo", "Pan", "Gorilla"]; source="ncbi")

# How much do the distances depend on the classification?
using Random
compare_sources(["Homo", "Pan", "Mus", "Gallus", "Danio", "Drosophila"];
                sources=["taxonomicon", "ott", "ncbi"], rng=Xoshiro(1))
```

`compare_sources()` resolves the same taxa in several sources and reports the
lineage depth of every taxon in each and the agreement between the distance
matrices (Spearman correlation and a Mantel permutation test). Distances are
comparable within one source, not across sources.

Requests to each source are spaced independently: at least
`Taxodist.MIN_INTERVAL[]` seconds (default 0.25) for The Taxonomicon and the
Open Tree of Life, and at least 0.4 s for NCBI, which allows three requests per
second without an API key. With an
[NCBI API key](https://www.ncbi.nlm.nih.gov/books/NBK25497/) (set it in the
`NCBI_API_KEY` environment variable or in `Taxodist.NCBI_API_KEY[]`) the limit
is ten per second and requests are spaced by 0.12 s. Refused requests (HTTP
429) and server errors are retried after a pause.

## Auditable and offline workflows

Resolve a batch once, inspect every decision, and reuse the stored lineages
without additional requests:

```julia
resolution = taxo_resolve(
    ["Tyrannosaurus", "Nereis", "50841"];
    ambiguity="warn",
    within="Animalia",   # optional: restrict homonyms to a clade
)

resolution.status
resolution.candidates
distance_matrix(resolution)
```

Local classifications can use the same workflow completely offline. Pass the
lineages as a vector of pairs to keep their order (and to give positional
`ids`); a `Dict` is accepted too, but its labels are sorted and `ids` must then
be a `Dict` keyed by label:

```julia
curated = taxo_from_lineages([
    "Alpha" => ["Biota", "Animalia", "Alpha"],
    "Beta" => ["Biota", "Animalia", "Beta"],
]; ids=["A1", "B1"], source="Curated study")
distance_matrix(curated)
```

Portable bundles combine the resolution table, matrix, metric definition, and
provenance. Their JSON schema is shared with the R and Python packages:

```julia
bundle = taxo_bundle(resolution; progress=false)
write_taxodist_bundle(bundle, "analysis.json")
restored = read_taxodist_bundle("analysis.json")

# Months later: what changed in the source classification?
compare_bundles("analysis.json", taxo_bundle(restored.resolution.input))
```

## The distance metric

`Taxodist.jl` measures relatedness by asking a single question: how deep is the most recent common ancestor (MRCA)?

$$
d(A,B) =
\begin{cases}
0, & A = B, \\
\dfrac{1}{\text{depth}(\text{MRCA}(A,B))}, & A \ne B.
\end{cases}
$$

The deeper the shared ancestor, the smaller the distance and the more related the two taxa are. A shallow MRCA means the two taxa diverged early; a deep MRCA means they share a long common history. Zero is reserved for identical hierarchy nodes. Consequently, a taxon and one of its descendants have a positive distance even though they are connected by ancestry. This distinction makes the measure a proper ultrametric on each connected hierarchy.

Distance and membership answer different questions. For example, *Tyrannosaurus* has a positive distance from *Dinosauria* because they are distinct nodes, while `is_member("Tyrannosaurus", "Dinosauria")` returns `true`. Use `is_member()` or `taxo_path()` when the relationship of interest is containment or ancestry.

The values represent classification depth in the chosen source, not evolutionary time or phylogenetic branch length. The Taxonomicon resolves lineages much more deeply than the Open Tree and NCBI taxonomies (*Tyrannosaurus* has over 70 nodes in its lineage), which makes the distances informative across all of life, including fossil groups; because lineage depths differ between sources, distances are comparable only within one source.

### Other metrics and ranks

`taxo_distance`, `distance_matrix`, `cross_distances`, `closest_relative`,
`focal_distances`, `taxo_heatmap`, and `taxo_bundle` accept a `metric`
keyword (also passed through by `taxo_cluster` and `taxo_ordinate`). With
`h` the MRCA depth and `d_A`, `d_B` the lineage depths:

| `metric` | distance | notes |
|---|---|---|
| `"reciprocal"` (default) | `1/h` | ultrametric |
| `"path"` | `d_A + d_B - 2h` | edges between the two nodes |
| `"wu_palmer"` | `1 - 2h/(d_A + d_B)` | normalized for local resolution; 1 for disconnected taxa |
| `"ranked"` | half the path through the standard Linnaean ranks | Clarke & Warwick step weights; needs ranks |

Every metric is 0 for identical lineages, and disconnected lineages are `Inf`
apart (except under `"wu_palmer"`). Matrices record their metric
(`distance_matrix(taxa; metric="path").metric`).

The rank of every node is recorded with the lineage, for every source.
`lineage_ranks("Passer")` returns a `DataFrame` of `node`, `rank`, and `depth`
(`rank` is `missing` for unranked nodes); `taxo_resolve` results gain a
`ranks` column, and `taxo_from_lineages` accepts
`ranks=["Alpha" => [missing, "Kingdom", "Genus"], ...]`. The `"ranked"` metric
uses `Taxodist.RANK_LEVELS[]` (kingdom, phylum, class, order, family, genus,
species by default).

Bundles use schema 1.1, which records the metric (and its rank levels) and
the lineage ranks; schema 1.0 bundles are still read.

## Community analyses

The functions below take a community as a vector of taxa, a `Dict` of
abundances, or a matrix / `DataFrame` with one row per sample, together with a
distance matrix whose labels include every taxon (from any source). The
example uses curated lineages, so it runs offline:

```julia
using DataFrames, Random

ranks = [missing, "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
birds = ["Biota", "Animalia", "Chordata", "Aves"]
mammals = ["Biota", "Animalia", "Chordata", "Mammalia"]
lineages = [
    "Passer" => [birds; "Passeriformes"; "Passeridae"; "Passer"],
    "Turdus" => [birds; "Passeriformes"; "Turdidae"; "Turdus"],
    "Corvus" => [birds; "Passeriformes"; "Corvidae"; "Corvus"],
    "Falco" => [birds; "Falconiformes"; "Falconidae"; "Falco"],
    "Homo" => [mammals; "Primates"; "Hominidae"; "Homo"],
    "Mus" => [mammals; "Rodentia"; "Muridae"; "Mus"],
]
example = taxo_from_lineages(lineages; ranks=[first(p) => ranks for p in lineages])
pool = distance_matrix(example; metric="ranked", progress=false)
community = DataFrame(sample=["site1", "site2", "site3"],
                      Passer=[3, 1, 0], Turdus=[1, 0, 0], Corvus=[0, 2, 0],
                      Falco=[1, 0, 1], Homo=[0, 0, 4], Mus=[0, 1, 2])

# Warwick and Clarke's taxonomic diversity and distinctness indices
taxo_distinctness(community, pool)
taxo_distinctness_test(["Passer", "Turdus", "Falco"], pool; n_iter=999, rng=Xoshiro(1))

# Taxonomic dissimilarity between samples (Γ+ or Θ+)
taxo_dissimilarity(community, pool; index="gamma_plus")

# Standardized effect sizes of MPD (NRI) and MNTD (NTI)
taxo_ses(community, pool; index="mpd", n_iter=999, rng=Xoshiro(1))
```

- `taxo_distinctness` computes Warwick and Clarke's Δ, Δ*, Δ+, Λ+, and SΔ+
  (results agree with `vegan::taxondive()`); `taxo_distinctness_test`
  compares Δ+ and Λ+ of a sample with random draws from a master list.
- `taxo_dissimilarity` computes the taxonomic dissimilarities Γ+ and Θ+ of
  Clarke, Somerfield and Chapman (2006) between samples, which credit closely
  related taxa as near matches, returned as a `TaxonomicDistanceMatrix` between
  samples for ordination, clustering, or PERMANOVA.
- `taxo_ses` computes the standardized effect sizes of the mean pairwise
  distance (NRI) and mean nearest-taxon distance (NTI) of Webb et al. (2002)
  under a species-pool (`null_model="pool"`) or independent-swap
  (`null_model="independent_swap"`) null model, optionally weighted by
  abundance as in picante.

## Caching

Lineages are cached in memory automatically during a session. To persist the
cache across sessions and avoid redundant network requests, use `save_cache()`
and `load_cache()`, which default to a per-user cache file
(`TAXODIST_CACHE_DIR` overrides the directory), or pass an explicit path.
Requests are retried automatically after transient network errors
(`Taxodist.MAX_TRIES[]`, default 3), and `check_coverage()` reports `missing`
when the server could not be reached. `taxonomicon_available()` checks quickly
whether The Taxonomicon can be reached.

## Plotting

Plots.jl is an optional dependency. Run `using Plots` to enable
`plot_taxodist_cluster()`, `plot_taxodist_ord()`, and `taxo_heatmap()`.
`taxo_cluster()` reproduces R's `hclust()` exactly (merges, heights, and leaf
order) for all its linkage methods.

## Data sources and citation

Published analyses should cite both the software and the classification the
lineages were taken from (`taxo_sources()` lists the citations, and every
resolution and bundle records its source):

- **The Taxonomicon** (`source="taxonomicon"`, the default):
  Brands, S.J. (ed.) 1989-present. *The Taxonomicon*. Universal Taxonomic
  Services, Zwaag, The Netherlands. http://taxonomicon.taxonomy.nl
- **Open Tree Taxonomy** (`source="ott"`), retrieved from the Open Tree of
  Life API (CC0; the version used is recorded in the resolution):
  Rees, J.A. & Cranston, K. (2017). Automated assembly of a reference
  taxonomy for phylogenetic data synthesis. *Biodiversity Data Journal* 5:
  e12581. https://doi.org/10.3897/BDJ.5.e12581
- **NCBI Taxonomy** (`source="ncbi"`), retrieved from the NCBI E-utilities:
  Schoch, C.L. et al. (2020). NCBI Taxonomy: a comprehensive update on
  curation, resources and tools. *Database* 2020: baaa062.
  https://doi.org/10.1093/database/baaa062

The NCBI E-utilities allow three requests per second without an API key and
ten with one; `Taxodist.jl` stays within these limits (see
[Choosing a taxonomic source](#choosing-a-taxonomic-source)). To cite
`Taxodist.jl` itself, use the metadata in `CITATION.cff` (GitHub's "Cite this
repository" button).

## Contributing

Found a taxon with an incorrect lineage? Please [open an issue](https://github.com/rodrigosqrt3/Taxodist.jl/issues)
with the queried names, the source, and the identifiers when available;
lineage corrections are the most valuable contribution to this package.
