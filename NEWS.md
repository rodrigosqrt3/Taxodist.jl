# Taxodist.jl 0.8.0

- Added `taxo_resolve()` with explicit resolved, ambiguous, unresolved, and
  retrieval-error states, retained candidates, selected identifiers, stored
  lineages, duplicate-input reuse, and configurable ambiguity policies.
- Added `taxo_from_lineages()` for curated, frozen, and fully offline
  classifications with source provenance.
- Extended `distance_matrix()` to consume `TaxodistResolution` objects without
  further retrieval and to return matrices stored in `TaxodistBundle` objects.
- Added `taxo_bundle()`, `validate_taxodist_bundle()`,
  `write_taxodist_bundle()`, and `read_taxodist_bundle()` using the same JSON
  schema version 1.0 as the R and Python implementations.
- Added portable handling for missing distances and positive/negative infinity,
  plus cross-language bundle reconstruction and scientific consistency checks.
- Updated the package version, user agent, citation metadata, documentation,
  and deterministic tests for the complete 0.8 API.

# Taxodist.jl 0.7.0

- Aligned live Taxonomicon lineage parsing with the R and Python
  implementations by reading the complete hierarchy content, including
  intermediate nodes that are rendered as plain text rather than hyperlinks.
- Aligned ambiguous-name candidate matching with the case-insensitive
  whole-word behavior used by R and Python.
- Added deterministic tests for structured hierarchy pages and exclusion of
  descendant taxa after the requested page subject.
- Added a numeric-only pairwise path for distance-matrix construction and a
  cache for final name-resolved lineages, avoiding repeated result allocation,
  filtering, and regular-expression work in warm calls.
- Preserved labeled matrix output and the continuous-prefix MRCA distance
  semantics shared with R and Python.
- Expanded deterministic coverage of cache validation and Gumbo node handling;
  the package test suite now reports complete source-line coverage.

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