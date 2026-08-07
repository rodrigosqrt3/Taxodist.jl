# Taxodist.jl <picture><source media="(prefers-color-scheme: dark)" srcset="images/taxodist_dark.png"><source media="(prefers-color-scheme: light)" srcset="images/taxodist_sepia.png"><img alt="taxodist logo" src="images/taxodist_sepia.png" align="right" height="200"></picture>

[![Julia](https://img.shields.io/badge/Julia-1.10%2B-9558B2?logo=julia)](https://julialang.org/) &nbsp; [![CI](https://github.com/rodrigosqrt3/Taxodist.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/rodrigosqrt3/Taxodist.jl/actions/workflows/CI.yml) &nbsp; [![codecov](https://codecov.io/gh/rodrigosqrt3/Taxodist.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rodrigosqrt3/Taxodist.jl) &nbsp; [![License: GPL v3+](https://img.shields.io/badge/License-GPL_v3%2B-blue.svg)](LICENSE.md)

**Taxonomic hierarchy distances derived from lineage classifications.**

Julia implementation of `taxodist`, a package for taxonomic hierarchy
distances and lineage analysis using ordered classifications retrieved from
The Taxonomicon. The Julia, R, and Python implementations share the same
distance definition and lineage semantics.

Development is currently focused on complete behavioral parity with the R and
Python implementations. The Julia package will not be registered until its
public API, edge-case behavior, reference data, and tests satisfy the shared
compatibility contract.

## Installation

Until registration in Julia's General Registry, install the development
version directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/rodrigosqrt3/Taxodist.jl")
```

## Distance definition

For two ordered lineages, the most recent common ancestor is the final node in
their continuous common prefix. The distance is defined as

```math
d(A,B)=
\begin{cases}
0, & L_A=L_B, \\
\dfrac{1}{\operatorname{depth}(\operatorname{MRCA}(A,B))}, & L_A\neq L_B.
\end{cases}
```

Lineages without a shared root have infinite distance. The values represent
classification depth, not evolutionary time or phylogenetic branch length.

## Current API

The current development version includes:

- taxon search and lineage retrieval from The Taxonomicon;
- validated in-memory and JSON cache management;
- pairwise distances and most recent common ancestors;
- labeled distance matrices;
- focal distances, closest-relative searches, lineage depth, and coverage
  checks;
- lineage comparison, shared-clade queries, membership filtering, and full
  paths through the MRCA.

```julia
using Taxodist

lineage = get_lineage("Carnotaurus")
comparison = taxo_distance("Carnotaurus", "Tyrannosaurus")

taxa = ["Carnotaurus", "Tyrannosaurus", "Triceratops"]
matrix = distance_matrix(taxa)
matrix["Carnotaurus", "Tyrannosaurus"]

closest_relative("Carnotaurus", ["Tyrannosaurus", "Triceratops"])
focal_distances("Carnotaurus", taxa)

shared_clades("Carnotaurus", "Tyrannosaurus")
is_member("Carnotaurus", "Theropoda")
taxo_path("Carnotaurus", "Triceratops").data
```

## Development

From a system terminal opened in the package directory:

```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

Alternatively, from inside the Julia REPL:

```julia
using Pkg
cd("path/to/taxodist-jl")
Pkg.activate(".")
Pkg.test()
```

Continuous integration runs the complete test suite on Linux, macOS, and
Windows with Julia 1.10 and the latest stable Julia release. Coverage from the
Linux job is reported to Codecov.

## Status

Retrieval, cache management, the distance kernel, labeled matrices,
multi-taxon helpers, and lineage utilities are implemented with offline
tests. Statistical analysis helpers, visualization, packaged reference data,
and release infrastructure are being ported incrementally.

## Related projects

- [R package](https://github.com/rodrigosqrt3/taxodist)
- [Python package](https://github.com/rodrigosqrt3/taxodist-py)
- [Documentation website](https://rodrigosqrt3.github.io/taxodist-site)

## Data source and citation

Retrieved lineage data originate from **The Taxonomicon**, based on *Systema
Naturae 2000*. Published analyses should cite both the software and the
underlying classification source. Citation metadata are provided in
`CITATION.cff`.

## License

Taxodist.jl is distributed under the GNU General Public License version 3 or
later.