# Taxonomic sources other than The Taxonomicon.
#
# Every source resolves names to records with the same structure as
# taxo_resolve() uses for The Taxonomicon, so the resulting
# `TaxodistResolution` works with every matrix, tree, bundle, and
# distinctness function. Identifiers carry a source prefix ("ott:770309",
# "ncbi:9606") so they are never confused with Taxonomicon IDs. Mirrors R's
# R/sources.R.

const _SOURCE_IDS = ("taxonomicon", "ott", "ncbi")

const _OTT_API = "https://api.opentreeoflife.org/v3"
const _OTT_TAXONOMY_URL = "https://tree.opentreeoflife.org/about/taxonomy-version"
const _NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

"""
    NCBI_API_KEY

API key sent with requests to the NCBI E-utilities, which raises their rate
limit from three to ten requests per second. With the default `nothing`, the
`NCBI_API_KEY` environment variable is used. Set it with
`Taxodist.NCBI_API_KEY[] = "your key"`.
"""
const NCBI_API_KEY = Ref{Any}(nothing)

# Performs an HTTP request: `(method, url, headers, body; kwargs...)`.
# Replaced in tests.
const _http_request = Ref{Function}(HTTP.request)

# Name candidate of an external source; `name` is `nothing` for an NCBI
# identifier whose lineage could not be retrieved.
const _SourceCandidate = NamedTuple{(:id, :name, :synonym),Tuple{String,Union{Nothing,String},Bool}}

"""
    taxo_sources()

List the taxonomic sources from which `taxo_resolve` can retrieve lineages,
as a `DataFrame` with one row per source and columns `id` (the value of the
`source` keyword), `name`, `url`, `id_prefix` (prefix of the identifiers
stored in resolutions), `ranks`, `licence`, and `citation`.

- `"taxonomicon"`: The Taxonomicon (Brands 1989-present), a single curated
  classification with very deep lineages, many unranked clades, and extensive
  fossil coverage. Retrieved from its web pages.
- `"ott"`: the Open Tree of Life reference taxonomy (Rees and Cranston
  2017), a synthesis of NCBI, GBIF, WoRMS, IRMNG, and other taxonomies,
  retrieved from the Open Tree of Life API. Its version is recorded in the
  resolution.
- `"ncbi"`: the NCBI Taxonomy database (Schoch et al. 2020), covering
  organisms with sequence data, retrieved from the NCBI E-utilities. Requests
  are limited to three per second, or ten with an API key (`NCBI_API_KEY`
  environment variable or `Taxodist.NCBI_API_KEY[]`); requests are spaced by
  at least 0.4 s (0.12 s with a key), and refused requests (HTTP 429) are
  retried.

Lineages are ordered from the root of each source (`"Biota"`, `"life"`, and
`"cellular organisms"` respectively) to the taxon. Clade names differ between
sources: animals are `"Animalia"` in The Taxonomicon but `"Metazoa"` in the
other two, which matters for the `within` keyword of `taxo_resolve`.
Distances are only comparable within one source; see `compare_sources`.
"""
function taxo_sources()
    return DataFrame(
        id=collect(_SOURCE_IDS),
        name=["The Taxonomicon", "Open Tree Taxonomy", "NCBI Taxonomy"],
        url=[
            "http://taxonomicon.taxonomy.nl",
            _OTT_TAXONOMY_URL,
            "https://www.ncbi.nlm.nih.gov/taxonomy",
        ],
        id_prefix=["", "ott:", "ncbi:"],
        ranks=fill("ranked and unranked clades", 3),
        licence=[
            "All rights reserved; cite the source",
            "CC0 (Open Tree of Life); underlying sources keep their terms",
            "Public domain (US government work)",
        ],
        citation=[
            "Brands, S.J. (ed.) 1989-present. The Taxonomicon. Universal Taxonomic Services, Zwaag, The Netherlands.",
            "Rees, J.A. & Cranston, K. (2017). Automated assembly of a reference taxonomy for phylogenetic data synthesis. Biodiversity Data Journal 5: e12581.",
            "Schoch, C.L. et al. (2020). NCBI Taxonomy: a comprehensive update on curation, resources and tools. Database 2020: baaa062.",
        ],
    )
end

const _SOURCE_CHOICES = join(repr.(_SOURCE_IDS), ", ", ", or ")

function _match_source(source)
    source isa Union{AbstractString,Symbol} || throw(ArgumentError(
        "source must be a single string, one of $(_SOURCE_CHOICES).",
    ))
    name = String(source)
    name in _SOURCE_IDS || throw(ArgumentError(
        "Unknown taxonomic source $(repr(name)). Use one of $(_SOURCE_CHOICES); " *
        "see taxo_sources().",
    ))
    return name
end

# -- HTTP ------------------------------------------------------------------

# Minimum interval between requests to each source. NCBI allows at most
# three requests per second without an API key and ten with one.
function _source_interval(source::AbstractString)
    interval = MIN_INTERVAL[]
    (isnan(interval) || interval < 0) && (interval = 0.0)
    if source == "ncbi"
        # NCBI counts requests in a sliding window, and bursts at exactly the
        # limit are refused with HTTP 429, so the floor is above 1/3 (1/10) s.
        interval = max(interval, isempty(_ncbi_api_key()) ? 0.4 : 0.12)
    end
    return interval
end

function _ncbi_api_key()
    key = NCBI_API_KEY[]
    key === nothing && (key = get(ENV, "NCBI_API_KEY", ""))
    return key isa AbstractString ? String(key) : ""
end

"""
    _source_request(source, method, url; body=nothing)

Perform one request to an external source and return `(status, text)`.
`status` is `nothing` when the request failed without a response. Requests
to each source are throttled independently, and transient failures (network
errors, 5xx responses, and HTTP 429) are retried with backoff; other client
errors are final. A
`body` is sent as JSON.
"""
function _source_request(source::AbstractString, method::AbstractString, url::AbstractString; body=nothing)
    headers = ["User-Agent" => _user_agent]
    payload = ""
    if body !== nothing
        payload = JSON3.write(body)
        push!(headers, "Content-Type" => "application/json")
    end
    interval = _source_interval(source)
    tries = max(MAX_TRIES[], 1)
    for attempt in 1:tries
        _throttle(source, interval)
        response = try
            _http_request[](
                method,
                url,
                headers,
                payload;
                status_exception=false,
                retry=false,
                readtimeout=30,
                connect_timeout=30,
            )
        catch error
            _is_network_error(error) || rethrow()
            attempt == tries && return (status=nothing, text=nothing)
            _sleep[](min(2^(attempt - 1), 5))
            continue
        end
        # Client errors are final, except 429 (too many requests), which is
        # retried like server errors.
        if (response.status >= 500 || response.status == 429) && attempt < tries
            _sleep[](min(2^(attempt - 1), 5))
            continue
        end
        return (status=Int(response.status), text=String(response.body))
    end
end

function _parse_json(text)
    text === nothing && return nothing
    try
        return JSON3.read(text)
    catch error
        error isa InterruptException && rethrow()
        return nothing
    end
end

# JSON body of a successful response, or `nothing`.
_response_json(response) = response.status == 200 ? _parse_json(response.text) : nothing

_json_get(data, key::Symbol) = data isa AbstractDict ? get(data, key, nothing) : nothing

const _UNRANKED_SOURCE_RANKS = Set(["no rank", "no rank - terminal", "cellular root", "acellular root"])

"""
    _normalize_source_rank(rank)

Rank labels as used by The Taxonomicon: capitalised, with unranked nodes
(other than explicit clades) as `missing`.
"""
function _normalize_source_rank(rank)
    (rank === nothing || ismissing(rank)) && return missing
    text = lowercase(String(_trimws(string(rank))))
    (isempty(text) || text in _UNRANKED_SOURCE_RANKS) && return missing
    return uppercasefirst(text)
end

function _source_cache_get(key::AbstractString)
    value = _cache_lookup(key)
    return value === _absent ? nothing : value
end

"""
    _source_cache_value(key, value)

Validate and convert a source entry of a cache file read by `load_cache` to
the form stored in the session cache.
"""
function _source_cache_value(key::AbstractString, value)
    invalid = ArgumentError("Invalid cache file: malformed source entry '$(key)'.")
    if occursin(r"^src_(ott|ncbi)_lin_", key)
        nodes = _json_get(value, :nodes)
        ranks = _json_get(value, :ranks)
        valid = nodes isa AbstractVector && all(node -> node isa AbstractString, nodes) &&
                ranks isa AbstractVector && length(ranks) == length(nodes) &&
                all(rank -> rank === nothing || rank isa AbstractString, ranks)
        valid || throw(invalid)
        return (nodes=String[String(node) for node in nodes], ranks=_rank_vector(ranks))
    elseif startswith(key, "src_ncbi_search_")
        (value isa AbstractVector && all(id -> id isa AbstractString, value)) || throw(invalid)
        return String[String(id) for id in value]
    elseif startswith(key, "src_ott_match_")
        valid = value isa AbstractVector && all(
            candidate -> _json_get(candidate, :id) isa AbstractString &&
                _json_get(candidate, :name) isa AbstractString &&
                _json_get(candidate, :synonym) isa Bool,
            value,
        )
        valid || throw(invalid)
        return _SourceCandidate[
            (id=String(candidate[:id]), name=String(candidate[:name]), synonym=candidate[:synonym])
            for candidate in value
        ]
    elseif key == "src_ott_about"
        (_json_get(value, :name) isa AbstractString && _json_get(value, :url) isa AbstractString) ||
            throw(invalid)
        return (name=String(value[:name]), url=String(value[:url]))
    end
    return value
end

# -- Open Tree Taxonomy ----------------------------------------------------

# Version of the Open Tree Taxonomy served by the API, as
# `(name="Open Tree Taxonomy 3.7 (ott3.7draft3)", url=...)`. Not cached when
# the request fails.
function _ott_about()
    cached = _source_cache_get("src_ott_about")
    cached === nothing || return cached
    data = _response_json(
        _source_request("ott", "POST", "$(_OTT_API)/taxonomy/about"; body=Dict{String,Any}()),
    )
    data isa AbstractDict || return (name="Open Tree Taxonomy", url=_OTT_TAXONOMY_URL)
    text(value) = value === nothing ? "" : string(value)
    weburl = _json_get(data, :weburl)
    about = (
        name="Open Tree Taxonomy $(text(_json_get(data, :version))) ($(text(_json_get(data, :source))))",
        url=weburl === nothing ? _OTT_TAXONOMY_URL : string(weburl),
    )
    _cache_set!("src_ott_about", about)
    return about
end

_plural(n::Integer, word::AbstractString) = n == 1 ? word : word * "s"

"""
    _ott_match_names(names; verbose=false)

Candidates of each name in the Open Tree Taxonomy, as a dictionary from name
to a vector of `(id, name, synonym)` (empty when the name is unknown); names
whose request failed are absent. Names are sent in batches of 250, without
approximate matching.
"""
function _ott_match_names(names; verbose::Bool=false)
    out = Dict{String,Vector{_SourceCandidate}}()
    todo = String[]
    for name in names
        cached = _source_cache_get("src_ott_match_" * lowercase(name))
        cached === nothing ? push!(todo, name) : (out[name] = cached)
    end
    for batch in Iterators.partition(todo, 250)
        verbose && println(
            "Matching $(length(batch)) $(_plural(length(batch), "name")) in the Open Tree Taxonomy...",
        )
        data = _response_json(_source_request(
            "ott", "POST", "$(_OTT_API)/tnrs/match_names";
            body=(names=collect(batch), do_approximate_matching=false),
        ))
        data === nothing && continue
        by_query = Dict{String,Any}()
        for result in something(_json_get(data, :results), [])
            by_query[lowercase(string(result[:name]))] = _json_get(result, :matches)
        end
        for name in batch
            matches = something(get(by_query, lowercase(name), nothing), [])
            candidates = _SourceCandidate[
                (
                    id="ott:" * string(hit[:taxon][:ott_id]),
                    name=String(hit[:taxon][:name]),
                    synonym=_json_get(hit, :is_synonym) === true,
                )
                for hit in matches
            ]
            _cache_set!("src_ott_match_" * lowercase(name), candidates)
            out[name] = candidates
        end
    end
    return out
end

"""
    _ott_lineage(id; verbose=false)

Lineage of an OTT taxon from the root (`"life"`) to the taxon, as
`(nodes, ranks)`, or `nothing` when it cannot be retrieved.
"""
function _ott_lineage(id::AbstractString; verbose::Bool=false)
    number = replace(id, r"^ott:" => "")
    key = "src_ott_lin_" * number
    cached = _source_cache_get(key)
    cached === nothing || return cached
    verbose && println("Retrieving lineage of ott:$(number)...")
    ott_id = tryparse(Int, number)
    ott_id === nothing && return nothing
    data = _response_json(_source_request(
        "ott", "POST", "$(_OTT_API)/taxonomy/taxon_info";
        body=(ott_id=ott_id, include_lineage=true),
    ))
    _json_get(data, :name) === nothing && return nothing
    ancestors = reverse(collect(something(_json_get(data, :lineage), [])))
    nodes = String[[String(node[:name]) for node in ancestors]; String(data[:name])]
    ranks = RankVector([_normalize_source_rank(_json_get(node, :rank)) for node in ancestors])
    push!(ranks, _normalize_source_rank(_json_get(data, :rank)))
    lineage = (nodes=nodes, ranks=ranks)
    _cache_set!(key, lineage)
    return lineage
end

# -- NCBI Taxonomy ---------------------------------------------------------

function _ncbi_url(endpoint::AbstractString, query::Vector{Pair{String,String}})
    query = vcat(query, ["tool" => "taxodist"])
    key = _ncbi_api_key()
    isempty(key) || push!(query, "api_key" => key)
    return "$(_NCBI_EUTILS)/$(endpoint)?" *
           join(["$(name)=$(HTTP.URIs.escapeuri(value))" for (name, value) in query], "&")
end

"""
    _ncbi_search(name; verbose=false)

NCBI identifiers (`"ncbi:<taxid>"`) matching a name: by scientific name,
then by any name (synonyms, common names). Empty when the name is unknown;
`nothing` when a request failed.
"""
function _ncbi_search(name::AbstractString; verbose::Bool=false)
    key = "src_ncbi_search_" * lowercase(name)
    cached = _source_cache_get(key)
    cached === nothing || return cached
    verbose && println("Searching NCBI Taxonomy for $(name)...")
    for field in ("Scientific Name", "All Names")
        term = "\"" * replace(name, "\"" => "") * "\"[" * field * "]"
        data = _response_json(_source_request(
            "ncbi", "GET",
            _ncbi_url("esearch.fcgi", [
                "db" => "taxonomy", "term" => term, "retmode" => "json", "retmax" => "50",
            ]),
        ))
        result = _json_get(data, :esearchresult)
        result isa AbstractDict || return nothing
        ids = String["ncbi:" * string(id) for id in something(_json_get(result, :idlist), [])]
        if !isempty(ids)
            _cache_set!(key, ids)
            return ids
        end
    end
    _cache_set!(key, String[])
    return String[]
end

# Retrieves and caches the lineages of NCBI taxa, in batches of 200.
function _ncbi_fetch(ids; verbose::Bool=false)
    numbers = unique([replace(id, r"^ncbi:" => "") for id in ids])
    todo = [number for number in numbers if _cache_lookup("src_ncbi_lin_" * number) === _absent]
    for batch in Iterators.partition(todo, 200)
        verbose && println(
            "Retrieving $(length(batch)) $(_plural(length(batch), "lineage")) from NCBI Taxonomy...",
        )
        response = _source_request(
            "ncbi", "GET",
            _ncbi_url("efetch.fcgi", [
                "db" => "taxonomy", "id" => join(batch, ","), "retmode" => "xml",
            ]),
        )
        response.status == 200 || continue
        _ncbi_parse_taxa(response.text)
    end
    return nothing
end

"""
    _ncbi_parse_taxa(text)

Parse an efetch `TaxaSet` and cache the lineage of every taxon in it, from
its `LineageEx` and itself. Merged identifiers (`AkaTaxIds`) answer with the
taxon they were merged into. Malformed XML is ignored.
"""
function _ncbi_parse_taxa(text::AbstractString)
    root = try
        _parse_xml(text)
    catch error
        error isa ArgumentError || rethrow()
        return nothing
    end
    root.name == "TaxaSet" || return nothing
    for taxon in _xml_children(root, "Taxon")
        ancestors = [
            ancestor for container in _xml_children(taxon, "LineageEx")
            for ancestor in _xml_children(container, "Taxon")
        ]
        nodes = String[
            [_xml_child_text(ancestor, "ScientificName") for ancestor in ancestors];
            _xml_child_text(taxon, "ScientificName")
        ]
        ranks = RankVector([_normalize_source_rank(_xml_child_text(ancestor, "Rank")) for ancestor in ancestors])
        push!(ranks, _normalize_source_rank(_xml_child_text(taxon, "Rank")))
        lineage = (nodes=nodes, ranks=ranks)
        _cache_set!("src_ncbi_lin_" * _xml_child_text(taxon, "TaxId"), lineage)
        for container in _xml_children(taxon, "AkaTaxIds"), old in _xml_children(container, "TaxId")
            _cache_set!("src_ncbi_lin_" * _xml_text(old), lineage)
        end
    end
    return nothing
end

function _ncbi_lineage(id::AbstractString; verbose::Bool=false)
    number = replace(id, r"^ncbi:" => "")
    key = "src_ncbi_lin_" * number
    _cache_lookup(key) === _absent && _ncbi_fetch([number]; verbose=verbose)
    return _source_cache_get(key)
end

# -- Minimal XML reader ----------------------------------------------------
#
# Enough of XML 1.0 for the NCBI E-utilities: elements, text, character and
# predefined entity references, CDATA, comments, processing instructions, and
# a document type declaration (skipped). Attributes are ignored. Malformed
# documents throw an ArgumentError.

struct _XmlElement
    name::String
    children::Vector{Any}  # _XmlElement or String
end

const _XML_ENTITIES = Dict("lt" => "<", "gt" => ">", "amp" => "&", "quot" => "\"", "apos" => "'")

function _xml_entity(reference::AbstractString)
    body = reference[2:(end - 1)]
    startswith(body, "#x") && return string(Char(parse(UInt32, body[3:end]; base=16)))
    startswith(body, "#") && return string(Char(parse(UInt32, body[2:end])))
    haskey(_XML_ENTITIES, body) || throw(ArgumentError("Unknown XML entity &$(body);"))
    return _XML_ENTITIES[body]
end

_xml_unescape(text::AbstractString) =
    replace(text, r"&(#x[0-9A-Fa-f]+|#[0-9]+|[A-Za-z]+);" => _xml_entity)

# Index just after the next `delimiter` at or after `i`.
function _xml_after(s::String, i::Int, delimiter::AbstractString)
    found = findnext(delimiter, s, i)
    found === nothing && throw(ArgumentError("Malformed XML: missing $(repr(delimiter))."))
    return last(found) + 1
end

_xml_at(s::String, i::Int, prefix::AbstractString) = startswith(SubString(s, i), prefix)

# Skips whitespace, comments, processing instructions, and document type
# declarations outside the root element.
function _xml_skip_misc(s::String, i::Int)
    while i <= ncodeunits(s)
        if codeunit(s, i) in (0x20, 0x09, 0x0a, 0x0d)
            i += 1
        elseif _xml_at(s, i, "<?")
            i = _xml_after(s, i, "?>")
        elseif _xml_at(s, i, "<!--")
            i = _xml_after(s, i, "-->")
        elseif _xml_at(s, i, "<!DOCTYPE")
            i = _xml_after(s, i, ">")
        else
            break
        end
    end
    return i
end

function _parse_xml(text::AbstractString)
    s = String(text)
    i = _xml_skip_misc(s, 1)
    (i <= ncodeunits(s) && codeunit(s, i) == UInt8('<')) ||
        throw(ArgumentError("Malformed XML: no root element."))
    root, i = _xml_element(s, i)
    _xml_skip_misc(s, i) > ncodeunits(s) ||
        throw(ArgumentError("Malformed XML: content after the root element."))
    return root
end

function _xml_flush!(element::_XmlElement, buffer::IOBuffer)
    text = String(take!(buffer))
    isempty(text) || push!(element.children, text)
    return nothing
end

# Parses the element starting at `s[i] == '<'`; returns it and the index
# after its end.
function _xml_element(s::String, i::Int)
    close = _xml_after(s, i, ">")
    tag = SubString(s, i + 1, prevind(s, close - 1))
    self_closing = endswith(tag, "/")
    found = match(r"^[^\s/>]+", tag)
    found === nothing && throw(ArgumentError("Malformed XML: invalid tag <$(tag)>."))
    element = _XmlElement(String(found.match), Any[])
    self_closing && return element, close
    i = close
    buffer = IOBuffer()
    while true
        i > ncodeunits(s) && throw(ArgumentError("Malformed XML: unclosed <$(element.name)>."))
        if codeunit(s, i) == UInt8('<')
            if _xml_at(s, i, "</")
                close = _xml_after(s, i, ">")
                strip(SubString(s, i + 2, prevind(s, close - 1))) == element.name ||
                    throw(ArgumentError("Malformed XML: mismatched </$(element.name)>."))
                _xml_flush!(element, buffer)
                return element, close
            elseif _xml_at(s, i, "<!--")
                i = _xml_after(s, i, "-->")
            elseif _xml_at(s, i, "<![CDATA[")
                stop = _xml_after(s, i, "]]>")
                write(buffer, SubString(s, i + 9, prevind(s, stop - 3)))
                i = stop
            elseif _xml_at(s, i, "<?")
                i = _xml_after(s, i, "?>")
            else
                _xml_flush!(element, buffer)
                child, i = _xml_element(s, i)
                push!(element.children, child)
            end
        else
            stop = something(findnext('<', s, i), ncodeunits(s) + 1)
            write(buffer, _xml_unescape(SubString(s, i, prevind(s, stop))))
            i = stop
        end
    end
end

_xml_children(element::_XmlElement, name::AbstractString) =
    [child for child in element.children if child isa _XmlElement && child.name == name]

# Text of an element and its descendants, as xml2's xml_text().
_xml_text(element::_XmlElement) =
    join(child isa _XmlElement ? _xml_text(child) : child for child in element.children)

# Text of the first child element named `name`; "" when there is none.
function _xml_child_text(element::_XmlElement, name::AbstractString)
    children = _xml_children(element, name)
    return isempty(children) ? "" : _xml_text(first(children))
end

# -- Resolution ------------------------------------------------------------

_source_lineage(source::AbstractString, id::AbstractString; verbose::Bool=false) =
    source == "ott" ? _ott_lineage(id; verbose=verbose) : _ncbi_lineage(id; verbose=verbose)

"""
    _resolve_external(names, source, within_clade, verbose, progress, failed)

Resolve unique, normalized names (or identifiers) in an external source.
Returns one record per name, with the fields used by `taxo_resolve`.
"""
function _resolve_external(names, source, within_clade, verbose, progress, failed)
    id_pattern = Regex("^($(source):)?[0-9]+\$")
    is_id = [occursin(id_pattern, name) for name in names]
    search_names = names[.!is_id]

    candidates_by_name = if source == "ott"
        _ott_match_names(search_names; verbose=verbose)
    else
        found = [_ncbi_search(name; verbose=verbose) for name in search_names]
        ids = reduce(vcat, [ids for ids in found if ids !== nothing]; init=String[])
        isempty(ids) || _ncbi_fetch(ids; verbose=verbose)
        by_name = Dict{String,Vector{_SourceCandidate}}()
        for (name, ids) in zip(search_names, found)
            ids === nothing && continue
            by_name[name] = _SourceCandidate[
                (id=id, name=_ncbi_candidate_name(id), synonym=false) for id in ids
            ]
        end
        by_name
    end

    records = Vector{Any}(undef, length(names))
    for (index, name) in enumerate(names)
        records[index] = if is_id[index]
            _resolve_external_id(name, source, verbose, failed)
        else
            _resolve_external_name(
                name, get(candidates_by_name, name, nothing), source, within_clade, verbose, failed,
            )
        end
        progress && println("Resolving taxa [$(index)/$(length(names))]: $(name)")
    end
    return records
end

function _ncbi_candidate_name(id::AbstractString)
    lineage = _source_cache_get("src_ncbi_lin_" * replace(id, r"^ncbi:" => ""))
    return lineage === nothing ? nothing : last(lineage.nodes)
end

function _external_record(input, lineage, candidates::DataFrame)
    n_candidates = nrow(candidates)
    return (
        input=input,
        resolved_name=last(lineage.nodes),
        id=candidates.id[1],
        status=n_candidates > 1 ? "ambiguous" : "resolved",
        n_candidates=n_candidates,
        lineage_depth=length(lineage.nodes),
        lineage=copy(lineage.nodes),
        candidates=candidates,
        ranks=copy(lineage.ranks),
    )
end

function _resolve_external_id(input, source, verbose, failed)
    id = startswith(input, source * ":") ? input : source * ":" * input
    lineage = _source_lineage(source, id; verbose=verbose)
    lineage === nothing && return failed(input, "retrieval_error")
    return _external_record(input, lineage, DataFrame(id=[id], name=[last(lineage.nodes)]))
end

_candidate_table(candidates) = DataFrame(
    id=String[candidate.id for candidate in candidates],
    name=String[candidate.name for candidate in candidates],
)

function _resolve_external_name(input, candidates, source, within_clade, verbose, failed)
    candidates === nothing && return failed(input, "retrieval_error")
    # NCBI candidates have no name when their lineage could not be retrieved.
    named = filter(candidate -> candidate.name !== nothing, candidates)
    !isempty(candidates) && isempty(named) && return failed(input, "retrieval_error")
    candidates = named
    isempty(candidates) && return failed(input, "unresolved")

    # Exact matches of the queried name take precedence over matches through
    # synonyms or other names, as the lineage match does for The Taxonomicon.
    exact = filter(candidate -> lowercase(candidate.name) == lowercase(input), candidates)
    isempty(exact) || (candidates = exact)
    # NCBI returns identifiers in no particular order; lower (older)
    # identifiers come first, which is deterministic and usually puts the
    # best-known taxon first.
    if source == "ncbi"
        candidates = sort(
            candidates;
            by=candidate -> parse(Float64, replace(candidate.id, r"^ncbi:" => "")),
            alg=Base.Sort.DEFAULT_STABLE,
        )
    end
    candidates = sort(candidates; by=candidate -> candidate.synonym, alg=Base.Sort.DEFAULT_STABLE)

    needs_all = within_clade !== nothing || length(candidates) > 1
    lineages = Any[
        (needs_all || index == 1) ? _source_lineage(source, candidate.id; verbose=verbose) : nothing
        for (index, candidate) in enumerate(candidates)
    ]
    if length(candidates) > 1
        keep = _outermost_candidates([lineage === nothing ? nothing : lineage.nodes for lineage in lineages])
        candidates = candidates[keep]
        lineages = lineages[keep]
    end
    if within_clade !== nothing
        any(lineage -> lineage !== nothing, lineages) ||
            return failed(input, "retrieval_error", _candidate_table(candidates))
        inside = [
            lineage !== nothing && within_clade in lowercase.(lineage.nodes)
            for lineage in lineages
        ]
        any(inside) || return failed(input, "unresolved")
        candidates = candidates[inside]
        lineages = lineages[inside]
    end
    lineage = first(lineages)
    lineage === nothing && return failed(input, "retrieval_error", _candidate_table(candidates))
    return _external_record(input, lineage, _candidate_table(candidates))
end

"""
    _outermost_candidates(lineages)

Candidates nested in another candidate of the same name, such as the
nominotypical subgenus *Mus (Mus)* inside the genus *Mus*, are not homonyms:
only the outermost is kept. A candidate is dropped when the lineage of
another candidate is a proper prefix of its own. Candidates without a
lineage (`nothing`) are kept. Returns a `Bool` mask.
"""
function _outermost_candidates(lineages)
    n = length(lineages)
    return Bool[
        lineages[i] === nothing || !any(
            j -> j != i && lineages[j] !== nothing &&
                length(lineages[j]) < length(lineages[i]) &&
                lineages[i][1:length(lineages[j])] == lineages[j],
            1:n,
        )
        for i in 1:n
    ]
end

function _source_metadata(source::AbstractString)
    source == "taxonomicon" && return (name="The Taxonomicon", url="http://taxonomicon.taxonomy.nl")
    source == "ott" && return _ott_about()
    return (name="NCBI Taxonomy", url="https://www.ncbi.nlm.nih.gov/taxonomy")
end

# -- Comparison of sources -------------------------------------------------

"""
    TaxodistSourceComparison

Result of `compare_sources`, with fields `sources`, `resolutions` and
`matrices` (dictionaries keyed by source), `shared`, `taxa`, `agreement`,
`metric`, and `permutations`.
"""
struct TaxodistSourceComparison
    sources::Vector{String}
    resolutions::Dict{String,TaxodistResolution}
    matrices::Dict{String,TaxonomicDistanceMatrix}
    shared::Vector{String}
    taxa::DataFrame
    agreement::DataFrame
    metric::String
    permutations::Int
end

"""
    compare_sources(taxa; sources=["taxonomicon", "ott", "ncbi"], metric="reciprocal",
                    ambiguity="first", within=nothing, permutations=999,
                    verbose=false, progress=true, rng=Random.default_rng())

Resolve the same taxa (a vector of at least three names; identifiers are
specific to one source) in several taxonomic sources (see `taxo_sources`)
and compare the resulting classifications: the depth of every lineage in each
source and the agreement between the distance matrices, measured by rank
correlation and a Mantel permutation test. This shows how much an analysis
depends on the classification it is based on.

`ambiguity`, `verbose`, and `progress` are passed to `taxo_resolve`
(`ambiguity="first"` reports ambiguity in the result rather than as
warnings). `within` is one clade name used for every source, or a dictionary
(or named tuple, or vector of pairs) with one clade name per source, since
clade names differ between sources (e.g. `"Animalia"` and `"Metazoa"`).
`permutations=0` skips the Mantel test. Pass `rng` (e.g. `Random.Xoshiro(1)`)
for reproducible p-values.

Returns a `TaxodistSourceComparison` with:

- `resolutions`: the `TaxodistResolution` of every source;
- `matrices`: distance matrices restricted to `shared`, the taxa resolved
  (or ambiguous) in every source, each input once;
- `taxa`: a `DataFrame` with one row per input and, for every source, its
  status (`status_<source>`), identifier (`id_<source>`), and lineage depth
  (`depth_<source>`);
- `agreement`: a `DataFrame` with one row per pair of sources (`source_a`,
  `source_b`), the number of taxa (`n_taxa`), the Spearman correlation
  between their distances (`spearman`), and the Mantel p-value (`p_value`,
  `missing` when not computed).

The Mantel test permutes the taxa of one matrix and reports the proportion
of permutations whose correlation is at least the observed one, counting the
observed arrangement, i.e. `(1 + count) / (1 + permutations)`. Lineage depths
are not comparable across sources in absolute terms; correlations of
distances are.
"""
function compare_sources(
    taxa;
    sources=["taxonomicon", "ott", "ncbi"],
    metric="reciprocal",
    ambiguity="first",
    within=nothing,
    permutations=999,
    verbose::Bool=false,
    progress::Bool=true,
    rng::AbstractRNG=Random.default_rng(),
)
    metric = _match_metric(metric)
    valid_sources = sources isa AbstractVector && length(sources) >= 2 &&
                    all(source -> source isa Union{AbstractString,Symbol}, sources) &&
                    allunique(String.(sources))
    valid_sources || throw(ArgumentError("sources must name at least two different sources."))
    sources = String[_match_source(source) for source in sources]
    (taxa isa AbstractVector && length(taxa) >= 3 && all(taxon -> taxon isa AbstractString, taxa)) ||
        throw(ArgumentError("taxa must be a character vector of at least three names."))
    (permutations isa Real && !(permutations isa Bool) && isfinite(permutations) && permutations >= 0) ||
        throw(ArgumentError("permutations must be a single non-negative number."))
    permutations = trunc(Int, permutations)
    within_by_source = _within_by_source(within, sources)
    taxa = String.(taxa)

    resolutions = Dict{String,TaxodistResolution}()
    for source in sources
        resolutions[source] = taxo_resolve(
            taxa;
            ambiguity=ambiguity,
            verbose=verbose,
            progress=progress,
            within=within_by_source[source],
            source=source,
        )
    end

    usable(resolution) = [status in ("resolved", "ambiguous") for status in resolution.status]
    shared_rows = reduce(.&, [usable(resolutions[source]) for source in sources])
    # Identical inputs are compared once.
    seen = Set{String}()
    keep = falses(length(taxa))
    for (index, taxon) in enumerate(taxa)
        keep[index] = shared_rows[index] && !(taxon in seen)
        push!(seen, taxon)
    end
    shared = taxa[keep]

    matrices = Dict{String,TaxonomicDistanceMatrix}()
    for source in sources
        resolution = resolutions[source]
        subset = TaxodistResolution(
            resolution.data[keep, :], resolution.source, resolution.source_url, resolution.retrieved_at,
        )
        matrices[source] = distance_matrix(subset; progress=false, metric=metric)
    end

    table = DataFrame(input=taxa)
    for source in sources
        resolution = resolutions[source]
        table[!, "status_$(source)"] = String.(resolution.status)
        table[!, "id_$(source)"] = [id === nothing ? missing : id for id in resolution.id]
        table[!, "depth_$(source)"] = [depth === nothing ? missing : depth for depth in resolution.lineage_depth]
    end

    agreement = DataFrame(
        source_a=String[],
        source_b=String[],
        n_taxa=Int[],
        spearman=Union{Missing,Float64}[],
        p_value=Union{Missing,Float64}[],
    )
    for i in 1:(length(sources) - 1), j in (i + 1):length(sources)
        stats = _mantel(matrices[sources[i]], matrices[sources[j]], permutations, rng)
        push!(agreement, (sources[i], sources[j], length(shared), stats.r, stats.p))
    end

    return TaxodistSourceComparison(
        sources, resolutions, matrices, shared, table, agreement, metric, permutations,
    )
end

# The clade name passed as `within` to every source.
function _within_by_source(within, sources)
    within === nothing && return Dict{String,Any}(source => nothing for source in sources)
    within isa AbstractString && return Dict{String,Any}(source => String(within) for source in sources)
    entries = if within isa Union{AbstractDict,NamedTuple}
        collect(pairs(within))
    elseif within isa AbstractVector{<:Pair}
        within
    elseif within isa AbstractVector && length(within) == 1
        return _within_by_source(only(within), sources)
    else
        throw(ArgumentError("within must be one clade name or a vector named by source."))
    end
    lookup = Dict(string(first(entry)) => last(entry) for entry in entries)
    return Dict{String,Any}(source => get(lookup, source, nothing) for source in sources)
end

# Ranks with ties averaged, as R's rank().
function _average_ranks(values::AbstractVector)
    order = sortperm(values)
    ranks = Vector{Float64}(undef, length(values))
    i = 1
    while i <= length(values)
        j = i
        while j < length(values) && values[order[j + 1]] == values[order[i]]
            j += 1
        end
        ranks[order[i:j]] .= (i + j) / 2
        i = j + 1
    end
    return ranks
end

function _spearman(x::AbstractVector, y::AbstractVector)
    r = cor(_average_ranks(x), _average_ranks(y))
    return isnan(r) ? missing : r
end

"""
    _mantel(a, b, permutations, rng)

Spearman Mantel statistic between two distance matrices over the same taxa,
with a one-sided permutation p-value. Pairs with a non-finite distance in
either matrix are excluded. Returns `(r, p)`, `missing` when not computed.
"""
function _mantel(a::TaxonomicDistanceMatrix, b::TaxonomicDistanceMatrix, permutations::Integer, rng::AbstractRNG)
    x = _lower_triangle(a.values)
    y = _lower_triangle(b.values)
    ok = isfinite.(x) .& isfinite.(y)
    count(ok) < 3 && return (r=missing, p=missing)
    r = _spearman(x[ok], y[ok])
    (r === missing || permutations == 0) && return (r=r, p=missing)
    n = size(b.values, 1)
    hits = 0
    for _ in 1:permutations
        perm = randperm(rng, n)
        yp = _lower_triangle(b.values[perm, perm])
        okp = isfinite.(x) .& isfinite.(yp)
        rp = _spearman(x[okp], yp[okp])
        rp !== missing && rp >= r && (hits += 1)
    end
    return (r=r, p=(1 + hits) / (1 + permutations))
end

_format_median(value) = value === missing ? "NA" : isinteger(value) ? string(Int(value)) : string(value)

function Base.show(io::IO, ::MIME"text/plain", comparison::TaxodistSourceComparison)
    sources = comparison.sources
    println(io, "Comparison of taxonomic sources")
    println(io, "* Sources: ", join(sources, ", "))
    println(io, "* Metric: ", comparison.metric)
    println(
        io,
        "* Taxa resolved in every source: $(length(comparison.shared)) of $(nrow(comparison.taxa))",
    )
    medians = map(sources) do source
        depths = collect(skipmissing(comparison.taxa[!, "depth_$(source)"]))
        "$(source) $(_format_median(isempty(depths) ? missing : median(depths)))"
    end
    println(io, "Median lineage depth: ", join(medians, ", "))
    println(io)
    out = copy(comparison.agreement)
    out.spearman = [value === missing ? missing : round(value; digits=3) for value in out.spearman]
    show(io, MIME"text/plain"(), out; summary=false, eltypes=false, show_row_number=false)
end
