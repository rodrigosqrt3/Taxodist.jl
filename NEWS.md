# Taxodist.jl 1.0.0

First stable release, matching taxodist 1.0.0 for R. The default distance
metric, the other metrics, and the bundle format (schema 1.1) are now stable:
later 1.x releases will keep their definitions and continue to read bundles
written by this version.

## New features

- **Ranks.** The rank of every lineage node (e.g. `"Kingdom"`, `"Clade"`,
  `"Superfamily"`, `"Genus"`) is read from The Taxonomicon (structured and
  text parsers) and cached with the lineage (`ranks_<id>` / `rawranks_<id>`;
  missing ranks are `null` in cache files). The new `lineage_ranks()` returns
  a `DataFrame` of `node`, `rank` (`missing` when unknown), and `depth`.
  `taxo_resolve()` results gain a `ranks` column, and `taxo_from_lineages()`
  accepts a `ranks` keyword.
- **Distance metrics.** `taxo_distance()`, `distance_matrix()`,
  `cross_distances()`, `closest_relative()`, `focal_distances()`,
  `taxo_heatmap()`, and `taxo_bundle()` gain a `metric` keyword, also passed
  through by `taxo_cluster()` and `taxo_ordinate()`: `"reciprocal"` (default,
  `1/depth(MRCA)`), `"path"` (edges between the nodes), `"wu_palmer"`
  (`1 - 2 depth(MRCA) / (depth(A) + depth(B))`), and `"ranked"` (half the
  path through the standard Linnaean ranks, configurable with
  `Taxodist.RANK_LEVELS[]`; the step weights of Clarke and Warwick).
  `TaxonomicDistanceMatrix` gains a `metric` field (R's `"method"`
  attribute), and `taxo_distance()` results a `metric` element.
- **Taxonomic distinctness.** `taxo_distinctness()` computes Warwick and
  Clarke's Δ, Δ*, Δ+, Λ+, and SΔ+ from presence or abundance data (agreeing
  with `vegan::taxondive()`); `taxo_distinctness_test()` compares Δ+ and Λ+
  of a sample with random samples from a master list (`rng` keyword).
- `taxonomicon_available()` reports whether The Taxonomicon can be reached,
  with a short timeout, as in R.
- **Bundle schema 1.1.** Bundles record the metric (and, for `"ranked"`, its
  rank levels) and the lineage ranks (`lineage_ranks`, `null` for missing
  ranks). Distances are verified with the stored metric. Schema 1.0 bundles
  are still read but must use the reciprocal metric, and a bundle's `metric`
  must name a known metric.
- New dependencies: the `Random` and `Statistics` standard libraries.

## Breaking changes: lineage parsing

- Lineages are read from the structure of the Taxonomicon page
  (`#ctl00_cphTabContent_divHierarchy`) instead of its text: each node is the
  bold scientific name of a row, rows after the page's own taxon are dropped,
  and the page subject is appended when its own row is missing. Author
  strings, rank labels, `[plesion]`/`Anamorph genus` labels, and the `‡`
  marker no longer leak into or remove nodes (e.g. `"Passeres (C."`), and
  clades such as `Terrestrornithes` are no longer lost. **Distances can
  differ from earlier versions.**
- Species are kept as nodes of a queried subspecies, and `Genus (Subgenus)`
  nodes of a queried species of the same genus; other multi-word nodes and
  `[...]` nodes are dropped, and the lineage is truncated at the last match of
  the queried name.
- Pages without the hierarchy block are parsed from text exactly as in R:
  the scientific name is extracted from the start of each line (handling
  `de`, `van der`, `d'`, `Özdikmen`, `R.Br. ex`, ...), qualified prefixes
  (`[crown]`, `[plesion]`, `[stem]`) and rank prefixes require whitespace
  after the rank, the cutoff at the page subject is case-sensitive, and the
  page text is split into lines like R/libxml2 (the extra line breaks Julia
  used to insert at block tags were removed).
- Search results take the ID from the first TaxonTree link whose class is
  exactly `Valid`, and "astronomy" is no longer an astronomical keyword
  (matching R).
- The packaged `taxobase` data was rebuilt with the new parser.

## Bug fixes

- Names with regular-expression metacharacters (e.g. `"Microtus (Agricola)"`)
  are matched literally, and whole-word matching uses Unicode-aware
  lookarounds instead of `\b`.
- `taxo_resolve()` reports names whose candidates were all retrieved but none
  is biological as `"unresolved"` instead of `"retrieval_error"`.
- Among several biological candidates, names marked `nom. inval.`,
  `nom. rej.`, `nom. nud.`, or `nom. illeg.` are considered last (in
  `get_taxonomicon_id()` and `taxo_resolve()`).
- `save_cache()` records the parser version (`meta_parser_version`);
  `load_cache()` keeps only taxon IDs from caches written by older versions,
  with a warning.

## Julia-specific fixes

- The 30 s HTTP timeout was ignored (`request_timeout` exists only in
  HTTP.jl 2); requests now use `readtimeout` and `connect_timeout`.
- `_request_html` retries only network errors; interrupts and programming
  errors are rethrown.
- `taxo_from_lineages()` accepts an ordered vector of `label => lineage`
  pairs (or any ordered dictionary). A plain `Dict` is sorted by label and
  positional `ids` are rejected, because they could previously be assigned to
  the wrong taxa.
- Numeric taxa (`50841`, `50841.0`) no longer crash `taxo_distance()`,
  `mrca()`, `compare_lineages()`, `shared_clades()`, `taxo_path()`,
  `closest_relative()`, `focal_distances()`, `filter_clade()`, or
  `distance_matrix()`; integral floats are treated as integer IDs.
- `closest_relative()` and `focal_distances()` sort unresolved taxa (`NaN`)
  after taxa without a common ancestor (`Inf`), like R's `na.last`.
- `distance_matrix()` accepts `missing` entries (NaN rows, like R's `NA`).
- `taxo_resolve()` accepts a single string; `focal_distances()` gains
  `progress`.
- `taxo_cluster()` accepts R's `"ward.D2"`/`"ward.D"` names and lists the
  supported methods when given an unsupported one.
- Bundle reading rejects booleans as distances or counts, reports
  non-integral counts as `ArgumentError`, and checks that `candidates` (and
  `taxa`) are arrays.
- `summary_taxodist_ord()` reports `eig[1:k]` like R.
- The session cache is guarded by a lock, so lineages can be retrieved from
  several tasks or threads.
- `cophenetic()` no longer rebuilds a set per pair of tips.

## Network

- Consecutive requests are spaced by at least `Taxodist.MIN_INTERVAL[]`
  seconds (default 0.25; 0 disables it).

## Packaging

- Plots.jl is now a weak dependency loaded through the `TaxodistPlotsExt`
  extension: `using Taxodist` no longer loads Plots; run `using Plots` to
  enable plotting.
- Compat bounds for every dependency (including stdlibs and the test-only
  Test); `Plots = "1"`.
- `TAXODIST_VERSION` is taken from `Project.toml` via `pkgversion`.
- Added `LICENSE` (GPL-3) and `.gitignore`; removed stray coverage files.
- Shared parser fixtures (`test/fixtures`) with tests mirroring R's
  `test-parser.R`, plus `expected_text_fallback.json` generated with R to
  check the text fallback parser against R.

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