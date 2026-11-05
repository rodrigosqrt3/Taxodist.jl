const _taxodist_cache = Dict{String,Any}()

# Guards every access to `_taxodist_cache` made by the package, so lineages
# can be retrieved from several tasks or threads at once.
const _cache_lock = ReentrantLock()

const _absent = Symbol("taxodist.absent")

"""
    _cache_lookup(key)

Return the cached value for `key`, or `_absent` when there is none. A cached
`nothing` (a name known to be absent) is returned as `nothing`.
"""
_cache_lookup(key::AbstractString) = lock(() -> get(_taxodist_cache, key, _absent), _cache_lock)
_cache_set!(key::AbstractString, value) = lock(() -> (_taxodist_cache[key] = value), _cache_lock)
_cache_snapshot() = lock(() -> copy(_taxodist_cache), _cache_lock)

const _bare_ranks = Set([
    "Go to", "Superphylum", "Subfamily", "Suborder", "Epifamily",
    "Infraorder", "Superclass", "Subclass", "Superfamily", "Subgenus",
    "Section", "Division", "Candidatus", "Parvphylum", "Branch",
    "Supercohort", "Infracohort", "Subdivision", "Subsection", "Grade",
    "[unranked]", "(Supercluster)", "(Region)", "[crown]", "",
])

const _astronomical_pattern =
    r"\bastronomical\b|\bplanet\b|\bMinor planet\b|\bcomet\b|\basteroid\b"i

const _user_agent = "taxodist Julia package/$(TAXODIST_VERSION)"
const _http_get = Ref{Function}(HTTP.get)
const _sleep = Ref{Function}(sleep)

"""
    MAX_TRIES

Number of attempts made for each request to The Taxonomicon (default 3).
Network errors and server errors (5xx) are retried with exponential backoff;
client errors (4xx) are not. Change it with `Taxodist.MAX_TRIES[] = 5`.
"""
const MAX_TRIES = Ref{Int}(3)

"""
    MIN_INTERVAL

Minimum interval, in seconds, between consecutive requests to The Taxonomicon
(default 0.25), so that large batches do not overload the server. Change it
with `Taxodist.MIN_INTERVAL[] = 1.0`; a value of zero or less disables it.
"""
const MIN_INTERVAL = Ref{Float64}(0.25)

const _last_request = Ref{Float64}(-Inf)
const _throttle_lock = ReentrantLock()

function _throttle()
    interval = MIN_INTERVAL[]
    interval > 0 || return nothing
    lock(_throttle_lock) do
        wait = interval - (time() - _last_request[])
        wait > 0 && sleep(wait)
        _last_request[] = time()
    end
    return nothing
end

_normalise_space(text) = strip(replace(String(text), r"\s+" => " "))

# R's trimws(): removes ASCII spaces, tabs, and line breaks only.
_trimws(text) = strip(String(text), (' ', '\t', '\r', '\n'))

"""
    _taxon_label(taxon)

Convert a taxon name or numeric ID to a string. Integral numbers, including
floating-point values such as `50841.0`, become integer IDs as in R.
"""
_taxon_label(taxon::AbstractString) = String(taxon)
_taxon_label(taxon::Integer) = string(taxon)
_taxon_label(taxon::AbstractFloat) =
    isfinite(taxon) && isinteger(taxon) ? string(BigInt(taxon)) : string(taxon)
_taxon_label(taxon) = string(taxon)

const _regex_special_characters = Set(
    ['\\', '.', '^', '$', '*', '+', '?', '{', '}', '[', ']', '|', '(', ')'],
)

function _regex_escape(text::AbstractString)
    output = IOBuffer()
    for character in String(text)
        character in _regex_special_characters && write(output, '\\')
        write(output, character)
    end
    return String(take!(output))
end

# PCRE pattern matching `word` as a whole word. Unicode-aware lookarounds are
# used instead of `\b`, which fails next to non-word characters (e.g. a name
# ending in ")") and ignores non-ASCII letters.
_word_pattern(word::AbstractString) = "(*UCP)(?<!\\w)" * _regex_escape(word) * "(?!\\w)"

function _contains_word(text::AbstractString, word::AbstractString; ignore_case::Bool=false)
    return occursin(Regex(_word_pattern(word), ignore_case ? "i" : ""), text)
end

"""
    _finalize_lineage(lineage, taxon)

Truncate a name-resolved lineage at the queried taxon. Nodes starting with
`[` are dropped, and nodes containing spaces are kept only when they are the
queried taxon or one of its name-ancestors (see `_keep_lineage_node`). The
taxon is matched exactly first and case-insensitively as a fallback; the
lineage is truncated at its last match, and the taxon is appended only when it
is absent.
"""
function _finalize_lineage(lineage::AbstractVector, taxon::AbstractString)
    plan = _finalize_plan(lineage, taxon)
    nodes = String.(lineage[plan.keep])
    plan.append && push!(nodes, String(taxon))
    return nodes
end

"""
    _finalize_ranks(ranks, lineage, taxon)

Apply the truncation of `_finalize_lineage(lineage, taxon)` to the ranks
aligned with `lineage`. An appended taxon receives the rank implied by the
form of its name.
"""
function _finalize_ranks(ranks::AbstractVector, lineage::AbstractVector, taxon::AbstractString)
    plan = _finalize_plan(lineage, taxon)
    out = _rank_vector(ranks[plan.keep])
    plan.append && push!(out, _rank_from_name(taxon))
    return out
end

# Indices of `lineage` retained by `_finalize_lineage`, and whether the
# queried taxon must be appended.
function _finalize_plan(lineage::AbstractVector, taxon::AbstractString)
    candidates = [index for (index, node) in enumerate(lineage) if !startswith(node, "[")]
    candidates = candidates[_keep_lineage_node(lineage[candidates], taxon)]
    kept = lineage[candidates]
    matches = findall(==(taxon), kept)
    if isempty(matches)
        folded = lowercase(taxon)
        matches = findall(node -> lowercase(node) == folded, kept)
    end
    isempty(matches) || return (keep=candidates[1:last(matches)], append=false)
    return (keep=candidates, append=true)
end

"""
    _keep_lineage_node(lineage, taxon)

Nodes without spaces are always kept. A node with spaces (a species,
subspecies, or subgenus name) is kept only when it is the queried taxon or one
of its name-ancestors: the species of a queried subspecies, or the subgenus
`Genus (Subgenus)` of a queried species of the same genus.
"""
function _keep_lineage_node(lineage::AbstractVector, taxon::AbstractString)
    plain = [!occursin(' ', node) for node in lineage]
    occursin(' ', taxon) || return plain
    target = lowercase(taxon)
    genus = replace(target, r" .*$" => "")
    return [
        plain[index] ||
        lowercase(node) == target ||
        startswith(target, lowercase(node) * " ") ||
        (occursin(r"^\S+ \([^()]+\)$", node) && replace(lowercase(node), r" .*$" => "") == genus)
        for (index, node) in enumerate(lineage)
    ]
end

function clear_cache()
    lock(() -> empty!(_taxodist_cache), _cache_lock)
    return nothing
end

"""
    _default_cache_file()

Location of the persistent per-user cache file. `TAXODIST_CACHE_DIR`
overrides the directory; otherwise `%LOCALAPPDATA%\\taxodist\\Cache` on
Windows, `~/Library/Caches/taxodist` on macOS, and `\$XDG_CACHE_HOME/taxodist`
(default `~/.cache/taxodist`) elsewhere.
"""
_platform() = Sys.iswindows() ? :windows : Sys.isapple() ? :apple : :unix

function _default_cache_file(platform::Symbol=_platform())
    directory = get(ENV, "TAXODIST_CACHE_DIR", "")
    if isempty(directory)
        directory = if platform === :windows
            joinpath(get(ENV, "LOCALAPPDATA", homedir()), "taxodist", "Cache")
        elseif platform === :apple
            joinpath(homedir(), "Library", "Caches", "taxodist")
        else
            xdg = get(ENV, "XDG_CACHE_HOME", "")
            joinpath(isempty(xdg) ? joinpath(homedir(), ".cache") : xdg, "taxodist")
        end
    end
    return joinpath(directory, "lineage-cache.json")
end

# Version of the lineage parser. Saved caches record it so that lineages
# produced by an older parser are not reused.
const _parser_version = 2
const _parser_version_key = "meta_parser_version"

"""
    save_cache(file=nothing)

Write the session cache to a JSON file and return its path. Without a file,
a persistent per-user cache file is used (see `TAXODIST_CACHE_DIR`). The
version of the lineage parser is recorded in the file.
"""
function save_cache(file::Union{Nothing,AbstractString}=nothing)
    if file === nothing
        file = _default_cache_file()
        mkpath(dirname(file))
    end
    data = _cache_snapshot()
    data[_parser_version_key] = _parser_version
    open(file, "w") do stream
        JSON3.write(stream, data)
    end
    println("Cache saved to '$(file)' ($(length(data)) entries).")
    return String(file)
end

function _validated_cache(raw)
    data = Dict{String,Any}()

    for (raw_key, raw_value) in pairs(raw)
        key = String(raw_key)
        isempty(key) && throw(ArgumentError(
            "Invalid cache file: cache keys must be non-empty strings.",
        ))

        if startswith(key, "id_")
            value = raw_value isa AbstractString ? String(raw_value) : nothing
            (value === nothing || isempty(value)) && throw(ArgumentError(
                "Invalid cache file: taxon ID entries must be non-empty strings.",
            ))
            data[key] = value
        elseif startswith(key, "lin_") || startswith(key, "rawlin_")
            raw_value isa AbstractVector || throw(ArgumentError(
                "Invalid cache file: lineage entries must be sequences of strings.",
            ))
            all(node -> node isa AbstractString, raw_value) || throw(ArgumentError(
                "Invalid cache file: lineage entries must be sequences of strings.",
            ))
            data[key] = String.(collect(raw_value))
        elseif startswith(key, "ranks_") || startswith(key, "rawranks_")
            # Missing ranks are stored as JSON null.
            raw_value isa AbstractVector && all(
                rank -> rank === nothing || rank isa AbstractString, raw_value,
            ) || throw(ArgumentError(
                "Invalid cache file: rank entries must be sequences of strings or nulls.",
            ))
            data[key] = _rank_vector(raw_value)
        elseif startswith(key, "resolved_lineage_")
            raw_value isa AbstractVector || throw(ArgumentError(
                "Invalid cache file: resolved lineage entries must be sequences of strings.",
            ))
            all(node -> node isa AbstractString, raw_value) || throw(ArgumentError(
                "Invalid cache file: resolved lineage entries must be sequences of strings.",
            ))
            data[key] = String.(collect(raw_value))
        else
            data[key] = raw_value
        end
    end

    return data
end

"""
    load_cache(file=nothing)

Merge a cache written by `save_cache()` into the session cache after
validating it. Without a file, the persistent per-user cache is used. Lineages
from a cache written by an older parser are discarded (with a warning); their
taxon IDs are kept.
"""
function load_cache(file::Union{Nothing,AbstractString}=nothing)
    file === nothing && (file = _default_cache_file())
    isfile(file) || throw(ArgumentError("Cache file not found: '$(file)'"))
    raw = JSON3.read(read(file, String))
    raw isa JSON3.Object || throw(ArgumentError(
        "Invalid cache file: expected an object created by save_cache().",
    ))

    data = _validated_cache(raw)
    version = pop!(data, _parser_version_key, nothing)
    current = version isa Real && !(version isa Bool) && !isnan(version) &&
              version >= _parser_version
    if !current
        n_entries = length(data)
        filter!(entry -> startswith(first(entry), "id_"), data)
        n_dropped = n_entries - length(data)
        if n_dropped > 0
            @warn "The cache was written by an older version of taxodist. " *
                  "$(n_dropped) lineage $(n_dropped == 1 ? "entry was" : "entries were") " *
                  "discarded and will be retrieved again; taxon IDs were kept."
        end
    end
    lock(() -> merge!(_taxodist_cache, data), _cache_lock)
    println("Cache loaded from '$(file)' ($(length(data)) entries).")
    return nothing
end

"""
    cache_info()

Print and return cache statistics: `n_lineages`, `n_ids`, `taxa` (the
terminal node of each cached lineage, `nothing` for an empty lineage), `ids`
(the matching Taxonomicon IDs), and `size_bytes`.
"""
function cache_info()
    cache = _cache_snapshot()
    lineage_keys = sort([key for key in keys(cache) if startswith(key, "lin_")])
    id_keys = sort([key for key in keys(cache) if startswith(key, "id_")])
    ids = String[replace(key, r"^lin_" => "") for key in lineage_keys]
    taxa = Union{Nothing,String}[
        isempty(cache[key]) ? nothing : String(last(cache[key]))
        for key in lineage_keys
    ]
    size_bytes = Base.summarysize(cache)

    println("taxodist Cache")
    println("* Lineages cached : $(length(lineage_keys))")
    println("* IDs cached      : $(length(id_keys))")
    println("* Memory used     : $(size_bytes) bytes")

    if isempty(taxa)
        println("\nNo lineages cached yet.")
    else
        println("\nCached taxa:")
        for (taxon, id) in zip(taxa, ids)
            println("  $(something(taxon, "nothing")) (ID $(id))")
        end
    end

    return (
        n_lineages=length(lineage_keys),
        n_ids=length(id_keys),
        taxa=taxa,
        ids=ids,
        size_bytes=size_bytes,
    )
end

"""
    _is_network_error(error)

Whether `error` is a transient HTTP, socket, or I/O failure that a retry may
fix. Interrupts and programming errors are not.
"""
function _is_network_error(error)
    if error isa HTTP.Exceptions.RequestError || error isa HTTP.Exceptions.ConnectError
        error.error isa InterruptException && return false
    end
    return error isa HTTP.Exceptions.HTTPError || error isa Base.IOError || error isa EOFError
end

"""
    taxonomicon_available(; timeout=3) -> Bool

Send a single quick request to The Taxonomicon and report whether it answered
with a successful response. Returns `false` when there is no internet
connection or the server is down. Useful to skip code that needs the server;
lineage retrieval itself already handles failed requests.

`timeout` is the maximum time to wait for an answer, in seconds.

# Examples
```julia
taxonomicon_available()
```
"""
function taxonomicon_available(; timeout::Real=3)
    (timeout isa Bool || isnan(timeout) || timeout <= 0) && throw(ArgumentError(
        "timeout must be a single positive number of seconds.",
    ))
    response = try
        _http_get[](
            "http://taxonomicon.taxonomy.nl/",
            ["User-Agent" => _user_agent];
            status_exception=false,
            retry=false,
            readtimeout=ceil(Int, timeout),
            connect_timeout=ceil(Int, timeout),
        )
    catch error
        _is_network_error(error) || rethrow()
        return false
    end
    return response.status == 200
end

function _request_html(url::AbstractString; verbose::Bool=false)
    tries = max(MAX_TRIES[], 1)
    for attempt in 1:tries
        _throttle()
        response = try
            _http_get[](
                url,
                ["User-Agent" => _user_agent];
                status_exception=false,
                retry=false,
                readtimeout=30,
                connect_timeout=30,
            )
        catch error
            _is_network_error(error) || rethrow()
            if attempt == tries
                verbose && println("Could not reach Taxonomicon: $(sprint(showerror, error))")
                return nothing
            end
            _sleep[](min(2^(attempt - 1), 5))
            continue
        end
        if response.status >= 500 && attempt < tries
            _sleep[](min(2^(attempt - 1), 5))
            continue
        end
        if response.status != 200
            verbose && println("Could not reach Taxonomicon (HTTP $(response.status)).")
            return nothing
        end
        return String(response.body)
    end
end

_parse_document(html::AbstractString) = parsehtml(String(html); preserve_whitespace=true)

# Elements after which (or as first child of which) libxml2 keeps
# whitespace-only text; see `_xml_text`.
const _libxml_pcdata_tags = Set([
    :a, :abbr, :acronym, :address, :applet, :b, :bdo, :big, :blockquote,
    :body, :button, :caption, :center, :cite, :code, :dd, :del, :dfn, :div,
    :dt, :em, :font, :form, :h1, :h2, :h3, :h4, :h5, :h6, :i, :iframe, :ins,
    :kbd, :label, :legend, :li, :map, :menu, :object, :ol, :p, :pre, :q, :s,
    :samp, :small, :span, :strike, :strong, :td, :th, :tt, :u, :ul, :var,
])

_is_blank_text(text::AbstractString) = all(in((' ', '\t', '\n', '\r')), text)

# Whether libxml2's HTML parser, as used by R's xml2/rvest (`read_html()`
# with the NOBLANKS option), drops a whitespace-only text node: it does so
# inside <html> or <head>, and when the node follows an element (or opens an
# element) that does not normally contain text, such as <img> or <br>.
function _libxml_drops_blank(parent::Gumbo.HTMLElement, children, index::Int)
    Gumbo.tag(parent) in (:html, :head) && return true
    index == 1 && return !(Gumbo.tag(parent) in _libxml_pcdata_tags)
    previous = children[index - 1]
    # Gumbo siblings are either text or elements.
    previous isa HTMLText && return false
    return !(Gumbo.tag(previous) in _libxml_pcdata_tags)
end

function _write_xml_text(io::IO, node::Gumbo.HTMLElement)
    children = Gumbo.children(node)
    for (index, child) in enumerate(children)
        if child isa HTMLText
            _is_blank_text(child.text) && _libxml_drops_blank(node, children, index) && continue
            write(io, child.text)
        elseif child isa Gumbo.HTMLElement
            _write_xml_text(io, child)
        end
    end
end

"""
    _xml_text(node)

Text content of an element as returned by R's `rvest::html_text()`: the
concatenated text of its descendants, without the whitespace-only text nodes
that libxml2 removes while parsing. Matching this matters for the text
fallback parser, which splits the text into lines.
"""
_xml_text(node::Gumbo.HTMLElement) = sprint(_write_xml_text, node)

"""
    _parse_search_results(html)

Parse the rows of a TaxonList search page into unique candidate IDs and their
descriptions, as R's `.parse_search_rows()`: astronomical homonyms are
skipped, and the ID is taken from the first TaxonTree link whose class is
exactly `Valid`; rows whose link has no numeric ID are skipped.
"""
function _parse_search_results(html::AbstractString)
    document = _parse_document(html)
    results = NamedTuple{(:id, :name),Tuple{String,String}}[]
    seen = Set{String}()

    for row in eachmatch(Selector("tr"), document.root)
        text_entry = _normalise_space(_xml_text(row))
        occursin(_astronomical_pattern, text_entry) && continue
        links = [
            link for link in eachmatch(Selector("a"), row)
            if occursin("TaxonTree", Gumbo.getattr(link, "href", ""))
        ]
        isempty(links) && continue
        valid = findfirst(link -> Gumbo.getattr(link, "class", nothing) == "Valid", links)
        valid === nothing && continue
        match_id = match(r"id=([0-9]+)", Gumbo.getattr(links[valid], "href", ""))
        match_id === nothing && continue
        id = String(match_id.captures[1])

        text_entry = replace(text_entry, r"^N\|T\|P\|R\|B\|L\s*" => "")
        id in seen && continue
        push!(results, (id=id, name=String(text_entry)))
        push!(seen, id)
    end

    return results
end

const _search_parser = Ref{Function}(_parse_search_results)

function _taxo_search_details(taxon; verbose::Bool=false)
    taxon = _normalise_space(taxon)
    verbose && println("Searching Taxonomicon for '$(taxon)'...")
    safe_taxon = HTTP.URIs.escapeuri(String(taxon))
    url = "http://taxonomicon.taxonomy.nl/TaxonList.aspx?" *
          "subject=Entity&by=ScientificName&search=$(safe_taxon)"
    html = _request_html(url; verbose=verbose)
    html === nothing && return (status="retrieval_error", results=nothing)

    parsed = try
        _search_parser[](html)
    catch error
        error isa InterruptException && rethrow()
        verbose && println("Could not parse the Taxonomicon response: $(sprint(showerror, error))")
        return (status="retrieval_error", results=nothing)
    end
    if isempty(parsed)
        verbose && println("No matches found.")
        return (status="not_found", results=nothing)
    end

    results = DataFrame(parsed)
    verbose && println("Found $(nrow(results)) entries.")
    return (status="ok", results=results)
end

function taxo_search(taxon; verbose::Bool=false)
    return _taxo_search_details(taxon; verbose=verbose).results
end

"""
    _parse_hierarchy(document)

Parse the classification tree of a TaxonTree page. Each row of
`#ctl00_cphTabContent_divHierarchy` is an `<a>` (an ancestor) or a `<span>`
(the page's own taxon) that is not nested in another `a`/`span`, followed by a
`<br>`. The indentation of a row is the number of non-breaking spaces in the
text since the last `<br>`. The name is the text of the row's `<b>` element;
ranks, markers, and author strings live outside it. Rows after the page's own
taxon are its children and are excluded.

Returns `(nodes, ranks, self_missing)`, or `nothing` when the page has no such
block (or no named rows), so that the caller can fall back to text parsing.
The rank of a row is the text before its name (see `_rank_before_name`). When
the row of the page's own taxon is missing, the leading run of strictly
increasing indentation is kept and `self_missing` is `true`.
"""
function _parse_hierarchy(document)
    containers = eachmatch(Selector("#ctl00_cphTabContent_divHierarchy"), document.root)
    isempty(containers) && return nothing

    names = String[]
    ranks = Union{Missing,String}[]
    indents = Int[]
    is_self = Bool[]
    indent = 0

    function visit(node)
        for child in Gumbo.children(node)
            if child isa HTMLText
                indent += count(==(' '), child.text)
            elseif child isa HTMLElement
                node_tag = Gumbo.tag(child)
                if node_tag === :br
                    indent = 0
                elseif node_tag === :a || node_tag === :span
                    class = Gumbo.getattr(child, "class", nothing)
                    class in ("Ordinal", "Small") && continue
                    bold = eachmatch(Selector("b"), child)
                    isempty(bold) && continue
                    name = replace(_trimws(_xml_text(first(bold))), r"\s+" => " ")
                    isempty(name) && continue
                    push!(names, name)
                    push!(ranks, _rank_before_name(_xml_text(child), name))
                    push!(indents, indent)
                    push!(is_self, node_tag === :span)
                else
                    visit(child)
                end
            end
        end
    end
    visit(first(containers))
    isempty(names) && return nothing

    self_index = findfirst(is_self)
    if self_index !== nothing
        return (nodes=names[1:self_index], ranks=ranks[1:self_index], self_missing=false)
    end

    # Some pages (e.g. unverified taxa) omit the row of the page's own taxon.
    # Keep the leading run of strictly increasing indentation, i.e. the path
    # from the root, and let the caller append the subject name.
    keep = 1
    while keep < length(indents) && indents[keep + 1] > indents[keep]
        keep += 1
    end
    return (nodes=names[1:keep], ranks=ranks[1:keep], self_missing=true)
end

# -- Ranks --------------------------------------------------------------------

"""
    _rank_before_name(row_text, name)

Rank of a tree row: the text that precedes the taxon name (e.g.
`"[crown] Clade"`, `"Superfamily ‡?"`, `"Anamorph genus"`), without crown,
plesion, or stem qualifiers, extinction and type markers, quotes, or a
`"Taxon:"` caption. Rows without such text are species, subspecies, or
subgenera, whose rank is implied by the form of the name (see
`_rank_from_name`). Returns `missing` when no rank can be determined.
"""
function _rank_before_name(row_text::AbstractString, name::AbstractString)
    text = replace(String(row_text), r"\s+" => " ")
    position = findfirst(String(name), text)
    prefix = position === nothing ? "" : text[1:prevind(text, first(position))]
    return _normalize_rank(prefix, name)
end

function _normalize_rank(prefix::AbstractString, name::AbstractString)
    prefix = replace(String(prefix), r"^.*:" => ""; count=1)
    prefix = replace(prefix, r"\[(crown|plesion|stem)\]" => "")
    prefix = replace(prefix, r"[†‡ᵀ?\"]" => "")
    prefix = String(_trimws(replace(prefix, r"\s+" => " ")))
    return isempty(prefix) ? _rank_from_name(name) : prefix
end

"""
    _rank_from_name(name)

Rank implied by the form of a scientific name: `"Genus (Subgenus)"` is a
subgenus, a binomial a species, and a trinomial a subspecies. Uninomials carry
no rank information (`missing`).
"""
function _rank_from_name(name::AbstractString)
    occursin(r"^\S+ \([^()]+\)$", name) && return "Subgenus"
    words = length(split(_trimws(name), r"\s+"; keepempty=false))
    words == 2 && return "Species"
    words == 3 && return "Subspecies"
    return missing
end

const RankVector = Vector{Union{Missing,String}}

"""
    _rank_vector(ranks)

Convert ranks to a `Vector{Union{Missing,String}}`; `nothing` and `missing`
entries (JSON `null`, R's `NA`) become `missing`.
"""
_rank_vector(ranks) = Union{Missing,String}[
    (rank === nothing || ismissing(rank)) ? missing : String(rank) for rank in ranks
]

"""
    _lineage_from_page(html_or_document, taxon_id)

Extract the raw lineage (before trimming at `Biota`) of a TaxonTree page.
Uses the structured hierarchy block when present and falls back to text
parsing otherwise. Repeated names are collapsed to their first occurrence.
See `_lineage_and_ranks_from_page` for the ranks of the nodes.
"""
_lineage_from_page(page, taxon_id) = _lineage_and_ranks_from_page(page, taxon_id).nodes

"""
    _lineage_and_ranks_from_page(html_or_document, taxon_id)

Like `_lineage_from_page`, returning `(nodes, ranks)` with the rank of every
node (`missing` when unknown), as R's `.lineage_from_page()` with its
`"ranks"` attribute.
"""
_lineage_and_ranks_from_page(html::AbstractString, taxon_id) =
    _lineage_and_ranks_from_page(_parse_document(html), taxon_id)

function _lineage_and_ranks_from_page(document::Gumbo.HTMLDocument, taxon_id)
    subject = eachmatch(Selector("#ctl00_divSubject b"), document.root)
    current_name = if isempty(subject)
        nothing
    else
        String(_trimws(replace(_trimws(_xml_text(first(subject))), r"[†‡ᵀ?]" => "")))
    end

    structured = _parse_hierarchy(document)
    nodes, ranks = if structured !== nothing
        nodes = copy(structured.nodes)
        ranks = copy(structured.ranks)
        if structured.self_missing && current_name !== nothing &&
           !isempty(current_name) && last(nodes) != current_name
            subject_block = eachmatch(Selector("#ctl00_divSubject"), document.root)
            subject_text = isempty(subject_block) ? "" : _xml_text(first(subject_block))
            push!(nodes, current_name)
            push!(ranks, _rank_before_name(subject_text, current_name))
        end
        nodes, ranks
    else
        parsed = _parse_lineage_text_ranked(document, current_name, string(taxon_id))
        parsed.nodes, parsed.ranks
    end
    seen = Set{String}()
    first_occurrence = BitVector(undef, length(nodes))
    for (index, node) in enumerate(nodes)
        first_occurrence[index] = !(node in seen)
        push!(seen, node)
    end
    return (nodes=nodes[first_occurrence], ranks=RankVector(ranks[first_occurrence]))
end

const _rank_alternatives = join([
    "Superdomain", "Domain", "Superkingdom", "Kingdom", "Subkingdom",
    "Infrakingdom", "Superphylum", "Phylum", "Subphylum", "Infraphylum",
    "Superclass", "Class", "Subclass", "Infraclass", "Supercohort", "Cohort",
    "Subcohort", "Magnorder", "Grandorder", "Superorder", "Order", "Suborder",
    "Infraorder", "Parvorder", "Legion", "Superfamily", "Family", "Subfamily",
    "Tribe", "Subtribe", "Genus", "Species", "Subspecies", "Clade",
    "Supergrade", "Grade", "Subgrade",
], "|")

const _qualified_prefix = Regex(
    "^\\[(crown|plesion|stem)\\]\\s+(?:(?:" * _rank_alternatives * ")\\s+)?",
)
const _rank_prefix = Regex("^(?:Anamorph\\s+\\w+|" * _rank_alternatives * ")\\s+")

# Lower-case words that start an author string rather than continue a
# scientific name (e.g. "de Blainville", "van der Hoeven", "R.Br. ex").
const _author_particles = [
    "de", "del", "della", "des", "di", "do", "dos", "da", "das", "du", "d",
    "van", "von", "der", "den", "ter", "zu", "zur", "la", "le", "les", "ex",
    "in", "et", "al", "auct", "sensu", "nom", "non", "emend", "corrig", "ms",
    "pro", "sp", "spp", "var", "nec", "fide",
]

const _scientific_name_pattern = Regex(
    "(*UCP)^[^\\s\",\\[(]+" *
    "(?:\\s(?!(?:" * join(_author_particles, "|") * ")(?![\\w'’-]))" *
    "[a-zà-ÿ][\\w-]*(?![\\w.'’]))" *
    "{0,2}",
)

"""
    _extract_scientific_name(text)

Keep only the scientific name at the start of a cleaned tree line: the first
word (a uninomial, possibly hyphenated) followed by up to two lower-case
epithets that are not author particles. Quoted and bracketed lines, and
`missing`, are returned unchanged.
"""
_extract_scientific_name(::Missing) = missing

function _extract_scientific_name(text::AbstractString)
    text = String(text)
    (startswith(text, "\"") || startswith(text, "[") || startswith(text, "(")) && return text
    found = match(_scientific_name_pattern, text)
    return found === nothing ? text : String(found.match)
end

_extract_scientific_name(texts::AbstractVector) = [_extract_scientific_name(text) for text in texts]

_marked_label(raw_text) = replace(replace(String(raw_text), r"[†‡ᵀ]" => ""), r"^\?\s*" => ""; count=1)

const _rank_label_pattern = Regex(
    "^(?:\\[(?:crown|plesion|stem)\\]\\s+)?(Anamorph\\s+\\w+|" * _rank_alternatives * ")\\s+",
)

"""
    _label_rank(raw_text)

Rank prefix of a text line of the fallback parser (R's `rank_match`), or
`missing` when the line does not start with a known rank.
"""
function _label_rank(raw_text)
    found = match(_rank_label_pattern, _marked_label(raw_text))
    return found === nothing ? missing : String(found.captures[1])
end

function _clean_lineage_label(raw_text)
    text = _marked_label(raw_text)
    text = replace(text, _qualified_prefix => ""; count=1)
    text = replace(text, _rank_prefix => ""; count=1)
    text = replace(text, r"[†‡?]" => "")
    return _extract_scientific_name(strip(text))
end

"""
    _parse_lineage_text(document, current_name, taxon_id)

Fallback parser for pages without the structured hierarchy block, mirroring
R's `.parse_lineage_text()`: the text of `#divPageContent` is split into lines
(one node per line) and cut at the first line containing the page's subject as
a whole word (case-sensitive). Without that block the TaxonTree links up to
the last link to the page itself are used. Rank prefixes and markers are then
removed and the scientific name is extracted from each line.
"""
_parse_lineage_text(document, current_name, taxon_id::AbstractString) =
    _parse_lineage_text_ranked(document, current_name, taxon_id).nodes

"""
    _parse_lineage_text_ranked(document, current_name, taxon_id)

`_parse_lineage_text` returning `(nodes, ranks)`. The rank of a line is its
rank prefix; lines without one get the rank implied by the form of the name.
"""
function _parse_lineage_text_ranked(document, current_name, taxon_id::AbstractString)
    content = eachmatch(Selector("#divPageContent"), document.root)

    texts = if current_name !== nothing && !isempty(content)
        tree_text = _trimws(_xml_text(first(content)))
        raw_lines = [String(_trimws(line)) for line in split(tree_text, '\n')]
        filter!(!isempty, raw_lines)
        tree_start = findfirst(line -> startswith(line, "Natura"), raw_lines)
        tree_start === nothing || (raw_lines = raw_lines[tree_start:end])
        searchable = [replace(line, r"[†ᵀ]" => "") for line in raw_lines]
        cutoff = findfirst(line -> _contains_word(line, current_name), searchable)
        cutoff === nothing ? raw_lines : raw_lines[1:cutoff]
    else
        links = [
            link for link in eachmatch(Selector("a"), document.root)
            if occursin("TaxonTree", Gumbo.getattr(link, "href", "")) &&
               occursin(r"id=[0-9]", Gumbo.getattr(link, "href", ""))
        ]
        own = Regex("id=" * _regex_escape(taxon_id) * "(\$|&)")
        own_indices = findall(link -> occursin(own, Gumbo.getattr(link, "href", "")), links)
        isempty(own_indices) || (links = links[1:maximum(own_indices)])
        [String(_trimws(_xml_text(link))) for link in links]
    end

    lineage = String[]
    ranks = RankVector()
    for raw_text in texts
        label = _clean_lineage_label(raw_text)
        label in _bare_ranks && continue
        occursin(r"^\s*$", label) && continue
        startswith(label, "\"") && continue
        startswith(label, "Population") && continue
        push!(lineage, label)
        rank = _label_rank(raw_text)
        push!(ranks, ismissing(rank) ? _rank_from_name(label) : rank)
    end
    return (nodes=lineage, ranks=ranks)
end

function _parse_lineage_html(html::AbstractString, taxon_id::AbstractString; clean::Bool=true)
    lineage = _lineage_from_page(html, taxon_id)
    if clean
        biota_index = findfirst(==("Biota"), lineage)
        biota_index === nothing || (lineage = lineage[biota_index:end])
    end
    return lineage
end

"""
    get_lineage_by_id(taxon_id; clean=true, verbose=false)

Retrieve the lineage of one Taxonomicon ID (a string or an integral number).
Cleaned and uncleaned lineages are cached separately (`lin_<id>` /
`rawlin_<id>`), so the value of `clean` used in an earlier call never affects
later calls.
"""
function get_lineage_by_id(taxon_id; clean::Bool=true, verbose::Bool=false)
    (taxon_id === nothing || ismissing(taxon_id)) && return nothing
    taxon_id isa AbstractVector && throw(ArgumentError(
        "taxon_id must be a single ID, not a vector of $(length(taxon_id)) IDs.",
    ))
    id = String(strip(_taxon_label(taxon_id)))
    isempty(id) && return nothing
    occursin(r"^[0-9]+$", id) || return nothing

    clean_key = "lin_$(id)"
    raw_key = "rawlin_$(id)"
    cache_key = clean ? clean_key : raw_key
    cached = _cache_lookup(cache_key)
    if cached !== _absent
        verbose && println("Using cached lineage for ID $(id)")
        return copy(cached)
    end

    url = "http://taxonomicon.taxonomy.nl/TaxonTree.aspx?id=$(id)&src=0"
    html = _request_html(url; verbose=verbose)
    html === nothing && return nothing

    # One request serves both variants, so cache both, with their ranks.
    parsed = _lineage_and_ranks_from_page(html, id)
    raw = parsed.nodes
    isempty(raw) && return nothing
    biota_index = findfirst(==("Biota"), raw)
    kept = biota_index === nothing ? (1:length(raw)) : (biota_index:length(raw))
    cleaned = raw[kept]
    lock(_cache_lock) do
        _taxodist_cache[raw_key] = raw
        _taxodist_cache[clean_key] = cleaned
        _taxodist_cache["rawranks_$(id)"] = parsed.ranks
        _taxodist_cache["ranks_$(id)"] = parsed.ranks[kept]
    end
    return copy(clean ? cleaned : raw)
end

_format_taxon_id(taxon_id) = String(_trimws(_taxon_label(taxon_id)))

"""
    _ranks_by_id(taxon_id; clean=true, verbose=false)

Ranks aligned with `get_lineage_by_id(taxon_id; clean)`. Lineages cached
without ranks (e.g. restored from a cache file of an older version) are
retrieved again.
"""
function _ranks_by_id(taxon_id; clean::Bool=true, verbose::Bool=false)
    lineage = get_lineage_by_id(taxon_id; clean=clean, verbose=verbose)
    lineage === nothing && return nothing
    id = _format_taxon_id(taxon_id)
    key = (clean ? "ranks_" : "rawranks_") * id
    if _cache_lookup(key) === _absent
        lock(_cache_lock) do
            delete!(_taxodist_cache, "lin_$(id)")
            delete!(_taxodist_cache, "rawlin_$(id)")
        end
        lineage = get_lineage_by_id(id; clean=clean, verbose=verbose)
        (lineage === nothing || _cache_lookup(key) === _absent) && return nothing
    end
    ranks = _cache_lookup(key)
    length(ranks) == length(lineage) || return RankVector(fill(missing, length(lineage)))
    return _rank_vector(ranks)
end

"""
    _cached_ranks(taxon_id, lineage)

Cached ranks of a cleaned lineage, without any network request. Returns
`missing` ranks when none are cached or they do not align with `lineage`.
"""
function _cached_ranks(taxon_id, lineage)
    ranks = _cache_lookup("ranks_" * _format_taxon_id(taxon_id))
    if ranks !== _absent && ranks !== nothing && length(ranks) == length(lineage)
        return _rank_vector(ranks)
    end
    return RankVector(fill(missing, length(lineage)))
end

"""
    lineage_ranks(taxon; clean=true, verbose=false)

Return each node of the lineage of a taxon (name or numeric ID) with its rank
in The Taxonomicon classification (e.g. `"Kingdom"`, `"Clade"`,
`"Superfamily"`, `"Genus"`), as a `DataFrame` with columns `node`, `rank`, and
`depth`, or `nothing` when the taxon cannot be found. The nodes are exactly
those returned by `get_lineage` for the same arguments. `rank` is `missing`
for nodes without a stated rank (such as the philosophical root nodes above
`Biota`).

Ranks are read from the label that precedes each name in the classification,
without the `[crown]`, `[plesion]`, and `[stem]` qualifiers or extinction
markers. Rows without a label are species, subspecies, or subgenera, whose
rank is inferred from the form of the name. Ranks are used by the `"ranked"`
distance metric (see `taxo_distance`).
"""
function lineage_ranks(taxon; clean::Bool=true, verbose::Bool=false)
    (taxon === nothing || ismissing(taxon)) && return nothing
    taxon isa AbstractVector && throw(ArgumentError(
        "taxon must be a single taxon name or ID.",
    ))
    taxon_name = String(_normalise_space(_taxon_label(taxon)))
    isempty(taxon_name) && return nothing
    is_id = occursin(r"^[0-9]+$", taxon_name)
    id = is_id ? taxon_name : get_taxonomicon_id(taxon_name; verbose=verbose)
    id === nothing && return nothing

    ranks = _ranks_by_id(id; clean=clean, verbose=verbose)
    lineage = get_lineage_by_id(id; clean=clean, verbose=verbose)
    (ranks === nothing || lineage === nothing) && return nothing
    if !is_id
        ranks = _finalize_ranks(ranks, lineage, taxon_name)
        lineage = _finalize_lineage(lineage, taxon_name)
    end
    return DataFrame(node=String.(lineage), rank=ranks, depth=collect(1:length(lineage)))
end

"""
    _candidate_order(descriptions)

Order of homonymous candidates: entries whose search description marks the
name as invalid, rejected, nude, or illegitimate (`nom. inval.`, `nom. rej.`,
`nom. nud.`, `nom. illeg.`) are moved after the others, keeping the search
order within each group.
"""
function _candidate_order(descriptions)
    invalid = [
        occursin(r"\bnom\.\s*(inval|rej|nud|illeg)\."i, String(description))
        for description in descriptions
    ]
    return vcat(findall(!, invalid), findall(invalid))
end

"""
    get_taxonomicon_id(taxon; verbose=false)

Return the Taxonomicon ID of a taxon name, or `nothing` when it is absent or
the server cannot be reached. Whitespace in the name is normalized. Only
candidates descending from `Biota` are considered; when several remain, those
whose lineage contains the name are preferred, invalid or rejected names are
considered last, and a warning is emitted.
"""
function get_taxonomicon_id(taxon; verbose::Bool=false)
    return _lookup_taxonomicon_id(taxon; verbose=verbose).id
end

"""
    _lookup_taxonomicon_id(taxon; verbose=false)

Workhorse for `get_taxonomicon_id()`. Returns `(id, status)` where status is
`"found"`, `"not_found"`, or `"retrieval_error"`.
"""
function _lookup_taxonomicon_id(taxon; verbose::Bool=false)
    taxon_name = String(_normalise_space(_taxon_label(taxon)))
    cache_key = "id_$(taxon_name)"
    cached = _cache_lookup(cache_key)
    if cached !== _absent
        verbose && println("Using cached ID for $(taxon_name)")
        return (id=cached, status=cached === nothing ? "not_found" : "found")
    end

    search = _taxo_search_details(taxon_name; verbose=verbose)
    if search.status == "retrieval_error"
        @warn "Cannot reach The Taxonomicon server. The website (taxonomy.nl) appears to be offline or unreachable. Please try again later."
        return (id=nothing, status="retrieval_error")
    end
    search.status == "not_found" && return (id=nothing, status="not_found")

    biological = NamedTuple[]
    lineage_failed = false
    for candidate in eachrow(search.results)
        lineage = get_lineage_by_id(candidate.id; clean=true, verbose=false)
        if lineage === nothing
            lineage_failed = true
            continue
        end
        "Biota" in lineage || continue
        push!(biological, (id=String(candidate.id), name=String(candidate.name), lineage=lineage))
    end

    if isempty(biological)
        return (id=nothing, status=lineage_failed ? "retrieval_error" : "not_found")
    end

    if length(biological) > 1
        exact = [
            candidate for candidate in biological
            if any(node -> _contains_word(node, taxon_name; ignore_case=true), candidate.lineage)
        ]
        isempty(exact) || (biological = exact)
        biological = biological[_candidate_order([candidate.name for candidate in biological])]
    end

    if length(biological) > 1
        others = join(["ID $(candidate.id): $(candidate.name)" for candidate in biological[2:end]], "; ")
        @warn "Multiple valid biological entries found for '$(taxon_name)'. " *
              "Using: $(biological[1].name) (ID: $(biological[1].id)). " *
              "Pass a numeric ID directly to select another entry. Other available IDs: $(others)"
    end

    final_id = biological[1].id
    _cache_set!(cache_key, final_id)
    verbose && println("Found $(taxon_name) with ID $(final_id)")
    return (id=final_id, status="found")
end

"""
    get_lineage(taxon; clean=true, verbose=false)

Retrieve the lineage of a single taxon name or numeric ID. Whitespace is
normalized and the queried name is located case-insensitively when the
lineage is truncated.
"""
function get_lineage(taxon; clean::Bool=true, verbose::Bool=false)
    (taxon === nothing || ismissing(taxon)) && return nothing
    taxon isa AbstractVector && throw(ArgumentError(
        "taxon must be a single taxon name or ID.",
    ))
    taxon_name = String(_normalise_space(_taxon_label(taxon)))
    isempty(taxon_name) && return nothing
    is_id = occursin(r"^[0-9]+$", taxon_name)
    resolved_key = "resolved_lineage_$(clean ? 1 : 0)_$(taxon_name)"
    if !is_id
        cached = _cache_lookup(resolved_key)
        if cached !== _absent
            verbose && println("Using cached resolved lineage for $(taxon_name)")
            return copy(cached)
        end
    end
    id = is_id ? taxon_name : get_taxonomicon_id(taxon_name; verbose=verbose)
    id === nothing && return nothing

    lineage = get_lineage_by_id(id; clean=clean, verbose=verbose)
    lineage === nothing && return nothing
    is_id && return lineage

    lineage = _finalize_lineage(lineage, taxon_name)
    _cache_set!(resolved_key, copy(lineage))
    return lineage
end

const RESOLUTION_COLUMNS = [
    "input",
    "resolved_name",
    "id",
    "status",
    "n_candidates",
    "lineage_depth",
    "lineage",
    "candidates",
]

_utc_timestamp() = Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-dd HH:MM:SS") * " UTC"

function _empty_candidates()
    return DataFrame(id=String[], name=String[])
end

"""
    TaxodistResolution

An auditable table of taxon-resolution decisions. Columns are available as
properties (for example, `resolution.status`) and the complete `DataFrame` is
stored in `resolution.data`. Source provenance is retained alongside the
table.
"""
struct TaxodistResolution
    data::DataFrame
    source::String
    source_url::Union{Nothing,String}
    retrieved_at::Union{Nothing,String}
end

function Base.getproperty(resolution::TaxodistResolution, name::Symbol)
    if name in fieldnames(TaxodistResolution)
        return getfield(resolution, name)
    end
    data = getfield(resolution, :data)
    name in propertynames(data) && return data[!, name]
    return getproperty(data, name)
end

function Base.propertynames(resolution::TaxodistResolution; private::Bool=false)
    fields = collect(fieldnames(TaxodistResolution))
    columns = Symbol.(names(getfield(resolution, :data)))
    return Tuple(unique(vcat(fields, columns)))
end

Base.size(resolution::TaxodistResolution) = size(resolution.data)
Base.size(resolution::TaxodistResolution, dimension::Integer) =
    size(resolution.data, dimension)
Base.length(resolution::TaxodistResolution) = nrow(resolution.data)
Base.getindex(resolution::TaxodistResolution, args...) = getindex(resolution.data, args...)

function Base.show(io::IO, ::MIME"text/plain", resolution::TaxodistResolution)
    counts = summary_counts(resolution)
    println(io, "$(size(resolution, 1))×$(size(resolution, 2)) TaxodistResolution")
    println(
        io,
        "resolved=$(counts.resolved), ambiguous=$(counts.ambiguous), " *
        "unresolved=$(counts.unresolved), retrieval_error=$(counts.retrieval_error)",
    )
    show(io, MIME"text/plain"(), resolution.data)
end

function summary_counts(resolution::TaxodistResolution)
    statuses = resolution.status
    return (
        resolved=count(==("resolved"), statuses),
        ambiguous=count(==("ambiguous"), statuses),
        unresolved=count(==("unresolved"), statuses),
        retrieval_error=count(==("retrieval_error"), statuses),
    )
end

function _empty_resolution_table(; ranks::Bool=false)
    table = DataFrame(
        input=String[],
        resolved_name=Union{Nothing,String}[],
        id=Union{Nothing,String}[],
        status=String[],
        n_candidates=Int[],
        lineage_depth=Union{Nothing,Int}[],
        lineage=Union{Nothing,Vector{String}}[],
        candidates=DataFrame[],
    )
    ranks && (table.ranks = Union{Nothing,RankVector}[])
    return table
end

_has_rank_field(row) = row isa DataFrameRow ? hasproperty(row, :ranks) : haskey(row, :ranks)

"""
    _make_resolution(rows; source, source_url=nothing, retrieved_at=missing, with_ranks=nothing)

Build a `TaxodistResolution` from row records. The optional `ranks` column
(one rank vector, or `nothing`, per row) is added when `with_ranks` is `true`,
or by default when any row has a `ranks` field.
"""
function _make_resolution(
    rows;
    source::AbstractString,
    source_url::Union{Nothing,AbstractString}=nothing,
    retrieved_at=missing,
    with_ranks::Union{Nothing,Bool}=nothing,
)
    with_ranks === nothing && (with_ranks = any(_has_rank_field, rows))
    table = _empty_resolution_table(; ranks=with_ranks)
    for row in rows
        record = (
            input=String(row.input),
            resolved_name=row.resolved_name === nothing ? nothing : String(row.resolved_name),
            id=row.id === nothing ? nothing : String(row.id),
            status=String(row.status),
            n_candidates=Int(row.n_candidates),
            lineage_depth=row.lineage_depth === nothing ? nothing : Int(row.lineage_depth),
            lineage=row.lineage === nothing ? nothing : String.(collect(row.lineage)),
            candidates=copy(row.candidates),
        )
        if with_ranks
            ranks = _has_rank_field(row) ? row.ranks : nothing
            record = merge(record, (ranks=ranks === nothing ? nothing : _rank_vector(ranks),))
        end
        push!(table, record)
    end
    return TaxodistResolution(
        table,
        String(source),
        source_url === nothing ? nothing : String(source_url),
        ismissing(retrieved_at) ? _utc_timestamp() :
            retrieved_at === nothing ? nothing : String(retrieved_at),
    )
end

function _copy_resolution_row(row)
    copied = (
        input=row.input,
        resolved_name=row.resolved_name,
        id=row.id,
        status=row.status,
        n_candidates=row.n_candidates,
        lineage_depth=row.lineage_depth,
        lineage=row.lineage === nothing ? nothing : copy(row.lineage),
        candidates=copy(row.candidates),
    )
    _has_rank_field(row) || return copied
    return merge(copied, (ranks=row.ranks === nothing ? nothing : copy(row.ranks),))
end

"""
    taxo_resolve(taxa; ambiguity="warn", verbose=false, progress=true, within=nothing)

Resolve names or numeric Taxonomicon IDs (a vector, or a single string) into
an auditable batch result. Every input is assigned one of `resolved`,
`ambiguous`, `unresolved`, or `retrieval_error`, while candidate tables and
selected lineages are retained. A name whose candidates were all retrieved
but none is biological is `unresolved`; `retrieval_error` means a search or a
candidate lineage could not be retrieved.

Names are matched after normalizing whitespace, and the stored lineage is
truncated at the queried taxon exactly as in `get_lineage()`, so
`distance_matrix(taxo_resolve(x))` equals `distance_matrix(x)`. When several
candidates remain, invalid or rejected names (`nom. inval.`, `nom. rej.`, ...)
are ranked last. `within` optionally restricts candidates to those whose
lineage contains a clade (exact, case-insensitive); names without such a
candidate are `unresolved`.
"""
taxo_resolve(taxon::AbstractString; kwargs...) = taxo_resolve([taxon]; kwargs...)

function taxo_resolve(
    taxa::AbstractVector;
    ambiguity="warn",
    verbose::Bool=false,
    progress::Bool=true,
    within::Union{Nothing,AbstractString}=nothing,
)
    policy = lowercase(String(ambiguity))
    policy in ("warn", "first", "error") || throw(ArgumentError(
        "ambiguity must be \"warn\", \"first\", or \"error\".",
    ))
    all(taxon -> taxon isa AbstractString, taxa) || throw(ArgumentError(
        "taxa must be a vector of strings.",
    ))
    labels = String.(taxa)
    all(taxon -> !isempty(strip(taxon)), labels) || throw(ArgumentError(
        "taxa cannot contain missing or empty values.",
    ))
    within !== nothing && isempty(strip(within)) && throw(ArgumentError(
        "within must be nothing or one non-empty clade name.",
    ))
    within_clade = within === nothing ? nothing : lowercase(_normalise_space(within))

    function failed_record(taxon, status, candidates=_empty_candidates())
        return (
            input=taxon,
            resolved_name=nothing,
            id=nothing,
            status=status,
            n_candidates=nrow(candidates),
            lineage_depth=nothing,
            lineage=nothing,
            candidates=candidates,
            ranks=nothing,
        )
    end

    function resolve_one(taxon::String)
        if occursin(r"^[0-9]+$", taxon)
            lineage = get_lineage_by_id(taxon; clean=true, verbose=verbose)
            lineage === nothing && return failed_record(taxon, "retrieval_error")
            resolved_lineage = String.(lineage)
            name = last(resolved_lineage)
            return (
                input=taxon,
                resolved_name=name,
                id=taxon,
                status="resolved",
                n_candidates=1,
                lineage_depth=length(resolved_lineage),
                lineage=resolved_lineage,
                candidates=DataFrame(id=[taxon], name=[name]),
                ranks=_cached_ranks(taxon, resolved_lineage),
            )
        end

        search = _taxo_search_details(taxon; verbose=verbose)
        search.status == "retrieval_error" && return failed_record(taxon, "retrieval_error")

        candidates = search.results
        if search.status == "not_found" || candidates === nothing || nrow(candidates) == 0
            return failed_record(taxon, "unresolved")
        end

        searched_candidates = copy(candidates)
        candidate_lineages = [
            get_lineage_by_id(identifier; clean=true, verbose=verbose)
            for identifier in candidates.id
        ]
        retrieved = [lineage !== nothing for lineage in candidate_lineages]
        valid = [lineage !== nothing && "Biota" in lineage for lineage in candidate_lineages]
        candidates = candidates[valid, :]
        candidate_lineages = candidate_lineages[valid]

        if nrow(candidates) == 0
            # Only a failed lineage retrieval is a retrieval error; candidates
            # that were retrieved but are not biological (e.g. astronomical
            # objects) mean the name is simply not a taxon.
            # As for names outside `within`, an unresolved record carries no
            # candidates (bundle validation requires this).
            all(retrieved) && return failed_record(taxon, "unresolved")
            return failed_record(taxon, "retrieval_error", searched_candidates)
        end

        if within_clade !== nothing
            inside = [within_clade in lowercase.(lineage) for lineage in candidate_lineages]
            any(inside) || return failed_record(taxon, "unresolved")
            candidates = candidates[inside, :]
            candidate_lineages = candidate_lineages[inside]
        end

        # When a broad search returns several entries, prefer candidates whose
        # lineage contains the complete queried name, then rank invalid or
        # rejected names last.
        if nrow(candidates) > 1
            exact = [
                any(node -> _contains_word(node, taxon; ignore_case=true), lineage)
                for lineage in candidate_lineages
            ]
            if any(exact)
                candidates = candidates[exact, :]
                candidate_lineages = candidate_lineages[exact]
            end
            preferred = _candidate_order(candidates.name)
            candidates = candidates[preferred, :]
            candidate_lineages = candidate_lineages[preferred]
        end

        # Apply the same truncation as get_lineage() so that name-based and
        # resolution-based matrices are identical.
        selected_lineage = _finalize_lineage(String.(candidate_lineages[1]), taxon)
        selected_ranks = _finalize_ranks(
            _cached_ranks(candidates.id[1], candidate_lineages[1]),
            candidate_lineages[1],
            taxon,
        )
        candidate_count = nrow(candidates)
        return (
            input=taxon,
            resolved_name=last(selected_lineage),
            id=String(candidates.id[1]),
            status=candidate_count > 1 ? "ambiguous" : "resolved",
            n_candidates=candidate_count,
            lineage_depth=length(selected_lineage),
            lineage=selected_lineage,
            candidates=copy(candidates),
            ranks=selected_ranks,
        )
    end

    keys_ = String[_normalise_space(label) for label in labels]
    unique_keys = unique(keys_)
    resolved_unique = Dict{String,Any}()
    for (index, key) in pairs(unique_keys)
        resolved_unique[key] = resolve_one(key)
        progress && println("Resolving taxa [$(index)/$(length(unique_keys))]: $(key)")
    end
    rows = [
        merge(_copy_resolution_row(resolved_unique[key]), (input=label,))
        for (label, key) in zip(labels, keys_)
    ]
    result = _make_resolution(
        rows;
        source="The Taxonomicon",
        source_url="http://taxonomicon.taxonomy.nl",
        with_ranks=true,
    )

    ambiguous_rows = findall(==("ambiguous"), result.status)
    if !isempty(ambiguous_rows)
        details = join(
            ["$(result.input[i]) ($(result.n_candidates[i]) candidates)" for i in ambiguous_rows],
            ", ",
        )
        message = "Ambiguous taxon names were resolved using the first candidate: " *
                  details * ". Inspect taxo_search(), use `within`, or pass numeric IDs explicitly."
        policy == "error" && throw(ArgumentError(message))
        policy == "warn" && @warn message
    end
    return result
end

"""
    taxo_from_lineages(lineages; ids=nothing, source="user-supplied", ranks=nothing)

Create an auditable offline resolution from root-to-tip lineages without
contacting an online taxonomy service.

`lineages` maps input labels to lineages. Pass an ordered collection (a
vector of `label => lineage` pairs, or an ordered dictionary such as
`OrderedDict`) to keep the input order. A plain `Dict` has no defined order,
so its labels are sorted; positional (vector) `ids` are then rejected because
they could not be matched to labels reliably; pass `ids` as a dictionary keyed
by label instead. Without `ids`, identifiers `custom:<label>` are created.

`ranks` optionally gives the rank of every node of every lineage (use
`missing` or `nothing` for unranked nodes), keyed by label like `ids` (a
dictionary or a vector of `label => ranks` pairs) or, for ordered lineages,
by position. Ranks are required only by the `"ranked"` distance metric.
Without them every rank is `missing`.
"""
function taxo_from_lineages(
    lineages::AbstractDict;
    ids=nothing,
    source::AbstractString="user-supplied",
    ranks=nothing,
)
    raw_labels = collect(keys(lineages))
    all(label -> label isa AbstractString, raw_labels) || throw(ArgumentError(
        "lineages must have unique, non-empty names.",
    ))
    if lineages isa Union{Dict,IdDict,WeakKeyDict}
        ids isa AbstractVector && throw(ArgumentError(
            "Positional ids cannot be used with an unordered Dict of lineages, whose " *
            "iteration order is undefined. Pass ids as a Dict keyed by lineage name, or " *
            "pass lineages as an ordered vector of `label => lineage` pairs.",
        ))
        ranks isa AbstractVector && !(ranks isa AbstractVector{<:Pair}) && throw(ArgumentError(
            "Positional ranks cannot be used with an unordered Dict of lineages. Pass " *
            "ranks as a Dict keyed by lineage name.",
        ))
        raw_labels = sort(raw_labels; by=String)
    end
    return _taxo_from_ordered(
        raw_labels, [lineages[label] for label in raw_labels];
        ids=ids, source=source, ranks=ranks,
    )
end

function taxo_from_lineages(
    lineages::AbstractVector{<:Pair};
    ids=nothing,
    source::AbstractString="user-supplied",
    ranks=nothing,
)
    raw_labels = [first(entry) for entry in lineages]
    all(label -> label isa AbstractString, raw_labels) || throw(ArgumentError(
        "lineages must have unique, non-empty names.",
    ))
    return _taxo_from_ordered(
        raw_labels, [last(entry) for entry in lineages];
        ids=ids, source=source, ranks=ranks,
    )
end

function _taxo_from_ordered(raw_labels, raw_lineages; ids, source, ranks=nothing)
    labels = String.(raw_labels)
    length(unique(labels)) == length(labels) || throw(ArgumentError(
        "lineages must have unique, non-empty names.",
    ))
    all(label -> !isempty(strip(label)), labels) || throw(ArgumentError(
        "lineages must have unique, non-empty names.",
    ))
    isempty(strip(source)) && throw(ArgumentError("source must be one non-empty string."))

    ordered_lineages = Vector{Vector{String}}()
    for lineage in raw_lineages
        lineage isa AbstractVector || throw(ArgumentError(
            "Every lineage must be a non-empty vector without missing or empty nodes.",
        ))
        all(node -> node isa AbstractString && !isempty(strip(String(node))), lineage) ||
            throw(ArgumentError(
                "Every lineage must be a non-empty vector without missing or empty nodes.",
            ))
        isempty(lineage) && throw(ArgumentError(
            "Every lineage must be a non-empty vector without missing or empty nodes.",
        ))
        push!(ordered_lineages, String.(lineage))
    end

    identifiers = if ids === nothing
        ["custom:$(label)" for label in labels]
    elseif ids isa AbstractDict
        all(label -> haskey(ids, label), raw_labels) || throw(ArgumentError(
            "Named ids must contain every lineage name.",
        ))
        [ids[label] for label in raw_labels]
    elseif ids isa AbstractVector
        collect(ids)
    else
        throw(ArgumentError("ids must be a vector or dictionary of strings."))
    end
    length(identifiers) == length(labels) || throw(ArgumentError(
        "ids must contain one unique, non-empty identifier per lineage.",
    ))
    all(identifier -> identifier isa AbstractString && !isempty(strip(String(identifier))), identifiers) ||
        throw(ArgumentError(
            "ids must contain one unique, non-empty identifier per lineage.",
        ))
    identifiers = String.(identifiers)
    length(unique(identifiers)) == length(identifiers) || throw(ArgumentError(
        "ids must contain one unique, non-empty identifier per lineage.",
    ))
    rank_vectors = _lineage_rank_vectors(ranks, raw_labels, ordered_lineages)

    rows = [begin
        lineage = ordered_lineages[index]
        name = last(lineage)
        (
            input=label,
            resolved_name=name,
            id=identifiers[index],
            status="resolved",
            n_candidates=1,
            lineage_depth=length(lineage),
            lineage=lineage,
            candidates=DataFrame(id=[identifiers[index]], name=[name]),
            ranks=rank_vectors[index],
        )
    end for (index, label) in pairs(labels)]

    return _make_resolution(rows; source=source, source_url=nothing, with_ranks=true)
end

# Rank vectors aligned with `lineages` from the `ranks` argument of
# taxo_from_lineages(): `nothing`, a dictionary or vector of pairs keyed by
# label, or a positional vector.
function _lineage_rank_vectors(ranks, raw_labels, lineages)
    ranks === nothing && return [RankVector(fill(missing, length(lineage))) for lineage in lineages]
    selected = if ranks isa AbstractDict || ranks isa AbstractVector{<:Pair}
        lookup = ranks isa AbstractDict ? ranks : Dict(first(entry) => last(entry) for entry in ranks)
        all(label -> haskey(lookup, label), raw_labels) || throw(ArgumentError(
            "Named ranks must contain every lineage name.",
        ))
        [lookup[label] for label in raw_labels]
    elseif ranks isa AbstractVector
        collect(ranks)
    else
        throw(ArgumentError(
            "ranks must be a dictionary or vector of rank vectors (one per lineage).",
        ))
    end
    valid = length(selected) == length(lineages) && all(
        index -> selected[index] isa AbstractVector &&
            length(selected[index]) == length(lineages[index]) &&
            all(rank -> rank === nothing || ismissing(rank) || rank isa AbstractString, selected[index]),
        eachindex(lineages),
    )
    valid || throw(ArgumentError(
        "ranks must give one rank (or missing) per node of every lineage.",
    ))
    return [_rank_vector(rank_vector) for rank_vector in selected]
end
