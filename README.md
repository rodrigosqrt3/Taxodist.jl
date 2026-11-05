# taxodist <picture><source media="(prefers-color-scheme: dark)" srcset="images/taxodist_dark.png"><source media="(prefers-color-scheme: light)" srcset="images/taxodist_sepia.png"><img alt="taxodist logo" src="images/taxodist_sepia.png" align="right" height="200"></picture>

[![version](https://juliahub.com/docs/Taxodist/version.svg)](https://juliahub.com/ui/Packages/General/Taxodist/) &nbsp; [![Julia Tests](https://github.com/rodrigosqrt3/Taxodist.jl/actions/workflows/julia.yml/badge.svg)](https://github.com/rodrigosqrt3/Taxodist.jl/actions/workflows/julia.yml) &nbsp; [![codecov](https://codecov.io/gh/rodrigosqrt3/Taxodist.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rodrigosqrt3/Taxodist.jl)

**Taxonomic hierarchy distance and lineage computation for any taxon on Earth.**

`Taxodist.jl` retrieves full hierarchical lineages from [The Taxonomicon](http://taxonomicon.taxonomy.nl) and computes an ultrametric distance between any two taxa: a pair of dinosaurs, a dinosaur and a fungus, two species of fly, or an oak tree and a human.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/rodrigosqrt3/taxodist-jl")
```

## Basic usage

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
resolution = taxo_from_lineages([
    "Alpha" => ["Biota", "Animalia", "Alpha"],
    "Beta" => ["Biota", "Animalia", "Beta"],
]; ids=["A1", "B1"], source="Curated study")
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

The Taxonomicon provides substantially deeper lineage resolution than other programmatic sources, e.g., *Tyrannosaurus* has over 70 nodes in its lineage, which is what makes the distances meaningful across all of life.

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

The rank of every node is read from The Taxonomicon and cached with the
lineage. `lineage_ranks("Passer")` returns a `DataFrame` of `node`, `rank`,
and `depth` (`rank` is `missing` for unranked nodes); `taxo_resolve` results
gain a `ranks` column, and `taxo_from_lineages` accepts
`ranks=["Alpha" => [missing, "Kingdom", "Genus"], ...]`. The `"ranked"` metric
uses `Taxodist.RANK_LEVELS[]` (kingdom, phylum, class, order, family, genus,
species by default).

### Taxonomic distinctness

`taxo_distinctness` computes Warwick and Clarke's Δ, Δ*, Δ+, Λ+, and SΔ+ for
one or more samples (a vector of taxa, a `Dict` of abundances, or a matrix /
`DataFrame` with one row per sample); results agree with
`vegan::taxondive()`. `taxo_distinctness_test` compares Δ+ and Λ+ of a sample
with random draws from a master list:

```julia
pool = distance_matrix(resolution; metric="ranked")
taxo_distinctness(Dict("Alpha" => 5, "Beta" => 2), pool).data
taxo_distinctness_test(["Alpha", "Beta"], pool; n_iter=999, rng=Xoshiro(1))
```

Bundles use schema 1.1, which records the metric (and its rank levels) and
the lineage ranks; schema 1.0 bundles are still read.

## Caching

Lineages are cached in memory automatically during a session. To persist the
cache across sessions and avoid redundant network requests, use `save_cache()`
and `load_cache()`, which default to a per-user cache file
(`TAXODIST_CACHE_DIR` overrides the directory), or pass an explicit path.
Requests are retried automatically after transient network errors
(`Taxodist.MAX_TRIES[]`, default 3), consecutive requests are spaced by at
least `Taxodist.MIN_INTERVAL[]` seconds (default 0.25; 0 disables it), and
`check_coverage()` reports `missing` when the server could not be reached.
`taxonomicon_available()` checks quickly whether the server can be reached.

## Plotting

Plots.jl is an optional dependency. Run `using Plots` to enable
`plot_taxodist_cluster()`, `plot_taxodist_ord()`, and `taxo_heatmap()`.

## Data source

All lineage data is sourced from **The Taxonomicon** (taxonomy.nl), based on *Systema Naturae 2000* by Sheila J. Brands (1989 onwards). Please cite this resource in any published work using `Taxodist.jl`:

> Brands, S.J. (1989 onwards). *Systema Naturae 2000*. Amsterdam, The Netherlands. Retrieved from The Taxonomicon, http://taxonomicon.taxonomy.nl.

## Contributing

Found a taxon with an incorrect lineage? Please [open an issue](https://github.com/rodrigosqrt3/taxodist-jl/issues),
lineage corrections are the most valuable contribution to this package.