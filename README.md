# taxodist <picture><source media="(prefers-color-scheme: dark)" srcset="images/taxodist_dark.png"><source media="(prefers-color-scheme: light)" srcset="images/taxodist_sepia.png"><img alt="taxodist logo" src="images/taxodist_sepia.png" align="right" height="200"></picture>

[![version](https://juliahub.com/docs/Taxodist/version.svg)](https://juliahub.com/ui/Packages/General/Taxodist/) &nbsp; [![Julia Tests](https://github.com/rodrigosqrt3/taxodist-jl/actions/workflows/julia.yml/badge.svg)](https://github.com/rodrigosqrt3/taxodist-jl/actions/workflows/julia.yml) &nbsp; [![codecov](https://codecov.io/gh/rodrigosqrt3/Taxodist.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rodrigosqrt3/Taxodist.jl)

**Taxonomic hierarchy distance and lineage computation for any taxon on Earth.**

`Taxodist.jl` retrieves full hierarchical lineages from [The Taxonomicon](http://taxonomicon.taxonomy.nl) and computes an ultrametric distance between any two taxa: a pair of dinosaurs, a dinosaur and a fungus, two species of fly, or an oak tree and a human.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/rodrigosqrt3/Taxodist.jl")
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
save_cache("my_cache.json")
load_cache("my_cache.json")
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

## Caching

Lineages are cached in memory automatically during a session. To persist the
cache across sessions and avoid redundant network requests, use
`save_cache("file.json")` and `load_cache("file.json")`.

## Data source

All lineage data is sourced from **The Taxonomicon** (taxonomy.nl), based on *Systema Naturae 2000* by Sheila J. Brands (1989 onwards). Please cite this resource in any published work using `Taxodist.jl`:

> Brands, S.J. (1989 onwards). *Systema Naturae 2000*. Amsterdam, The Netherlands. Retrieved from The Taxonomicon, http://taxonomicon.taxonomy.nl.

## Contributing

Found a taxon with an incorrect lineage? Please [open an issue](https://github.com/rodrigosqrt3/Taxodist.jl/issues),
lineage corrections are the most valuable contribution to this package.
