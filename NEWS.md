# Taxodist.jl 0.6.0

## Initial Julia implementation

- Ported taxon search, identifier resolution, lineage retrieval, and validated
  JSON cache management from the R and Python packages.
- Implemented the shared continuous-prefix MRCA distance definition, including
  identity, ancestor--descendant, disconnected-hierarchy, symmetry, and
  ultrametric edge cases.
- Added labeled distance matrices, focal comparisons, closest relatives,
  lineage depth, coverage checks, clade membership, lineage comparison, and
  taxonomic paths.
- Added hierarchical clustering, PCoA ordination, dendrograms, ordination
  plots, summaries, and taxonomic-distance heatmaps.
- Added `load_taxobase()` using the same exported reference JSON distributed
  with the Python implementation and generated from the R `taxobase` object.
- Added cross-platform GitHub Actions, Codecov reporting, and citation
  metadata.