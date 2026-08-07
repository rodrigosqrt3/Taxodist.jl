const _taxodist_cache = Dict{String,Any}()

const _rank_prefixes =
    "Clade|Kingdom|Phylum|Superphylum|Subphylum|Infraphylum|Class|Order|" *
    "Suborder|Infraorder|Parvorder|Grandorder|Magnorder|Cohort|Subcohort|" *
    "Legion|Family|Subfamily|Tribe|Subtribe|Genus|Species|Subkingdom|" *
    "Infrakingdom|Superclass|Subclass|Infraclass|Superorder|Superfamily|" *
    "Domain|Superkingdom|Grade|Subgrade|Supergrade"

const _bare_ranks = Set([
    "Go to", "Superphylum", "Subfamily", "Suborder", "Epifamily",
    "Infraorder", "Superclass", "Subclass", "Superfamily", "Subgenus",
    "Section", "Division", "Candidatus", "Parvphylum", "Branch",
    "Supercohort", "Infracohort", "Subdivision", "Subsection", "Grade",
    "[unranked]", "(Supercluster)", "(Region)", "[crown]", "",
])

const _astronomical_pattern =
    r"\bastronomical\b|\bplanet\b|\bMinor planet\b|\bcomet\b|\bastronomy\b|\basteroid\b"i

const _user_agent = "taxodist Julia package/0.6.0"

_normalise_space(text) = strip(replace(String(text), r"\s+" => " "))

function _clean_lineage_label(raw_text)
    text = replace(String(raw_text), r"[†ᵀ]" => "")
    text = _normalise_space(text)
    text = replace(text, Regex("^\\[crown\\]\\s+(" * _rank_prefixes * ")?\\s*") => "")
    text = replace(text, Regex("^(" * _rank_prefixes * ") ") => "")
    text = replace(text, r"\s+[A-Z][a-záàâãéèêíïóôõöúüç].*$" => "")
    text = replace(text, r"\s+[A-Z]\.[A-Z]\..*$" => "")
    text = replace(text, r"\s+auct\..*$" => "")
    text = replace(text, r"\s+von.*$" => "")
    text = replace(text, r"\s+\([A-Z][a-z].*$" => "")
    text = replace(text, r"\s+\(\d{4}\).*$" => "")
    text = replace(text, r"\s+\[.*$" => "")
    text = replace(text, r"\s+[A-Z]\.$" => "")
    text = replace(text, r"\s+\([a-z].*$" => "")
    text = replace(text, r"\s+\".*$" => "")
    text = replace(text, r"^\".*" => "")
    return strip(text)
end

function clear_cache()
    empty!(_taxodist_cache)
    return nothing
end

function save_cache(file::AbstractString)
    open(file, "w") do stream
        JSON3.write(stream, _taxodist_cache)
    end
    println("Cache saved to '$(file)' ($(length(_taxodist_cache)) entries).")
    return nothing
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
        elseif startswith(key, "lin_")
            raw_value isa AbstractVector || throw(ArgumentError(
                "Invalid cache file: lineage entries must be sequences of strings.",
            ))
            all(node -> node isa AbstractString, raw_value) || throw(ArgumentError(
                "Invalid cache file: lineage entries must be sequences of strings.",
            ))
            data[key] = String.(collect(raw_value))
        else
            data[key] = raw_value
        end
    end

    return data
end

function load_cache(file::AbstractString)
    isfile(file) || throw(ArgumentError("Cache file not found: '$(file)'"))
    raw = JSON3.read(read(file, String))
    raw isa JSON3.Object || throw(ArgumentError(
        "Invalid cache file: expected an object created by save_cache().",
    ))

    data = _validated_cache(raw)
    merge!(_taxodist_cache, data)
    println("Cache loaded from '$(file)' ($(length(data)) entries).")
    return nothing
end

function cache_info()
    lineage_keys = sort([key for key in keys(_taxodist_cache) if startswith(key, "lin_")])
    id_keys = sort([key for key in keys(_taxodist_cache) if startswith(key, "id_")])
    taxa = replace.(lineage_keys, r"^lin_" => "")
    size_bytes = Base.summarysize(_taxodist_cache)

    println("taxodist Cache")
    println("* Lineages cached : $(length(lineage_keys))")
    println("* IDs cached      : $(length(id_keys))")
    println("* Memory used     : $(size_bytes) bytes")

    if isempty(taxa)
        println("\nNo lineages cached yet.")
    else
        println("\nCached taxa:")
        foreach(taxon -> println("  $(taxon)"), taxa)
    end

    return (
        n_lineages=length(lineage_keys),
        n_ids=length(id_keys),
        taxa=taxa,
        size_bytes=size_bytes,
    )
end

function _request_html(url::AbstractString; verbose::Bool=false)
    try
        response = HTTP.get(
            url,
            ["User-Agent" => _user_agent];
            status_exception=false,
            retry=false,
            readtimeout=30,
        )
        if response.status != 200
            verbose && println("Could not reach Taxonomicon (HTTP $(response.status)).")
            return nothing
        end
        return String(response.body)
    catch error
        verbose && println("Could not reach Taxonomicon: $(sprint(showerror, error))")
        return nothing
    end
end

function _valid_tree_link(row)
    for link in eachmatch(Selector("a.Valid"), row)
        href = Gumbo.getattr(link, "href", "")
        occursin("TaxonTree", href) || continue
        match_id = match(r"id=([0-9]+)", href)
        match_id === nothing || return (id=match_id.captures[1], href=href)
    end
    return nothing
end

function _parse_search_results(html::AbstractString)
    document = parsehtml(String(html))
    results = NamedTuple{(:id, :name),Tuple{String,String}}[]
    seen = Set{String}()

    for row in eachmatch(Selector("tr"), document.root)
        text_entry = _normalise_space(nodeText(row))
        occursin(_astronomical_pattern, text_entry) && continue
        tree_link = _valid_tree_link(row)
        tree_link === nothing && continue
        tree_link.id in seen && continue

        text_entry = replace(
            text_entry,
            r"^N\s*\|\s*T\s*\|\s*P\s*\|\s*R\s*\|\s*B\s*\|\s*L\s*" => "",
        )
        push!(results, (id=tree_link.id, name=text_entry))
        push!(seen, tree_link.id)
    end

    return results
end

function taxo_search(taxon; verbose::Bool=false)
    verbose && println("Searching Taxonomicon for '$(taxon)'...")
    safe_taxon = HTTP.URIs.escapeuri(String(taxon))
    url = "http://taxonomicon.taxonomy.nl/TaxonList.aspx?" *
          "subject=Entity&by=ScientificName&search=$(safe_taxon)"
    html = _request_html(url; verbose=verbose)
    html === nothing && return nothing

    parsed = _parse_search_results(html)
    if isempty(parsed)
        verbose && println("No matches found.")
        return nothing
    end

    results = DataFrame(parsed)
    verbose && println("Found $(nrow(results)) entries.")
    return results
end

function _parse_lineage_html(html::AbstractString, taxon_id::AbstractString; clean::Bool=true)
    document = parsehtml(String(html))
    links = eachmatch(Selector("a"), document.root)
    texts = String[]

    for link in links
        href = Gumbo.getattr(link, "href", "")
        occursin("TaxonTree", href) || continue
        occursin(r"id=[0-9]+", href) || continue
        push!(texts, _normalise_space(nodeText(link)))
        link_id = match(r"id=([0-9]+)", href)
        link_id !== nothing && link_id.captures[1] == taxon_id && break
    end

    lineage = String[]
    seen = Set{String}()
    for raw_text in texts
        label = _clean_lineage_label(raw_text)
        label in _bare_ranks && continue
        isempty(label) && continue
        startswith(label, "\"") && continue
        startswith(label, "Population") && continue
        label in seen && continue
        push!(lineage, label)
        push!(seen, label)
    end

    if clean
        biota_index = findfirst(==("Biota"), lineage)
        biota_index === nothing || (lineage = lineage[biota_index:end])
    end

    return lineage
end

function get_lineage_by_id(taxon_id; clean::Bool=true, verbose::Bool=false)
    id = strip(String(taxon_id))
    isempty(id) && return nothing
    occursin(r"^[0-9]+$", id) || return nothing

    cache_key = "lin_$(id)"
    if haskey(_taxodist_cache, cache_key)
        verbose && println("Using cached lineage for ID $(id)")
        return copy(_taxodist_cache[cache_key])
    end

    url = "http://taxonomicon.taxonomy.nl/TaxonTree.aspx?id=$(id)&src=0"
    html = _request_html(url; verbose=verbose)
    html === nothing && return nothing

    lineage = _parse_lineage_html(html, id; clean=clean)
    isempty(lineage) && return nothing
    _taxodist_cache[cache_key] = lineage
    return copy(lineage)
end

function get_taxonomicon_id(taxon; verbose::Bool=false)
    taxon_name = String(taxon)
    cache_key = "id_$(taxon_name)"
    if haskey(_taxodist_cache, cache_key)
        verbose && println("Using cached ID for $(taxon_name)")
        return _taxodist_cache[cache_key]
    end

    candidates = taxo_search(taxon_name; verbose=verbose)
    candidates === nothing && return nothing

    biological = NamedTuple[]
    for candidate in eachrow(candidates)
        lineage = get_lineage_by_id(candidate.id; clean=true, verbose=false)
        lineage === nothing && continue
        "Biota" in lineage || continue
        push!(biological, (id=String(candidate.id), name=String(candidate.name)))
    end

    isempty(biological) && return nothing

    if length(biological) > 1
        exact = [candidate for candidate in biological if begin
            lineage = get_lineage_by_id(candidate.id; clean=true, verbose=false)
            lineage !== nothing && any(
                node -> lowercase(node) == lowercase(taxon_name),
                lineage,
            )
        end]
        isempty(exact) || (biological = exact)
    end

    if length(biological) > 1
        @warn "Multiple valid biological entries found for '$(taxon_name)'. Using ID $(biological[1].id). Pass a numeric ID directly to select another entry."
    end

    final_id = biological[1].id
    _taxodist_cache[cache_key] = final_id
    verbose && println("Found $(taxon_name) with ID $(final_id)")
    return final_id
end

function get_lineage(taxon; clean::Bool=true, verbose::Bool=false)
    taxon_name = String(taxon)
    is_id = occursin(r"^[0-9]+$", taxon_name)
    id = is_id ? taxon_name : get_taxonomicon_id(taxon_name; verbose=verbose)
    id === nothing && return nothing

    lineage = get_lineage_by_id(id; clean=clean, verbose=verbose)
    lineage === nothing && return nothing
    is_id && return lineage

    if !occursin(r"\s", taxon_name)
        lineage = [node for node in lineage if !occursin(" ", node) && !startswith(node, "[")]
        target_index = findlast(==(taxon_name), lineage)
    else
        lineage = [node for node in lineage if !occursin(" ", node) || node == taxon_name]
        target_index = findfirst(==(taxon_name), lineage)
    end

    if target_index === nothing
        push!(lineage, taxon_name)
    else
        lineage = lineage[1:target_index]
    end

    return isempty(lineage) ? nothing : lineage
end