# Taxodist.jl 0.9.0

This release mirrors taxodist 0.9.0 for R. See `PORTING-0.9.0.md` for the
correspondence between the three implementations.

## Bug fixes

- `taxo_resolve()` now truncates each selected lineage at the queried taxon
  exactly as `get_lineage()` does, so `distance_matrix(taxo_resolve(x))`
  equals `distance_matrix(x)` and `resolved_name` reports the queried taxon
  (previously the family for genera such as *Periplaneta*).
- Cleaned and uncleaned lineages are cached separately (`lin_<id>` /
  `rawlin_<id>`); `get_lineage_by_id(id; clean=false)` no longer poisons later
  cleaned lookups. One request fills both entries.
- Taxon names are normalized before lookup and located case-insensitively when
  truncating the lineage; `get_lineage("tyrannosaurus")` no longer appends a
  spurious node.
- `taxo_heatmap()` passes plot keyword arguments to `Plots.heatmap` (where they
  override the defaults) instead of to `distance_matrix`, which rejected them;
  `verbose` and `progress` are explicit and resolutions are accepted.
- `cache_info()` reports the names of cached taxa and their IDs (`ids`).
- `get_lineage_by_id()` accepts integer IDs and `nothing`, and raises an
  `ArgumentError` for a vector of IDs; `get_lineage()` likewise rejects
  vectors and returns `nothing` for `nothing`/`missing`.
- `get_taxonomicon_id()` now warns when The Taxonomicon cannot be reached, as
  in R and Python.

## Behavior changes

- `check_coverage()` reports `missing` in `covered` when the server could not
  be reached, instead of `false`.
- `save_cache()` returns the path of the written file.
- Displaying a bundle or passing it to `distance_matrix()` validates its
  structure only; `write_taxodist_bundle()` and `read_taxodist_bundle()` still
  recompute every distance (`validate_taxodist_bundle(bundle;
  check_distances=true)`).

## Performance

- Distance matrices compare integer prefix identifiers instead of strings,
  with type-stable loops.
- HTTP requests are retried with exponential backoff after network errors and
  server errors (5xx); `Taxodist.MAX_TRIES[]` sets the number of attempts
  (default 3).
- Each candidate lineage is retrieved once while resolving a name.

## New features

- `taxo_tree()` builds a `TaxodistTree` directly from lineages, with named
  clades as internal nodes; `newick(tree)` exports it and `cophenetic(tree)`
  reproduces `distance_matrix()` exactly.
- `cross_distances(x, y)` returns a labeled `TaxonomicCrossDistances` matrix
  between two sets of taxa.
- `compare_bundles(old, new)` reports which taxa changed status, identifier,
  or lineage between two bundles (or bundle files) and which distances changed.
- `taxo_resolve()` and `taxo_bundle()` gain `within` to restrict candidates to
  a clade.
- `save_cache()` and `load_cache()` default to a persistent per-user cache
  file (`TAXODIST_CACHE_DIR` overrides the directory).

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