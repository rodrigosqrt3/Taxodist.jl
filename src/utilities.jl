"""
    compare_lineages(taxon_a, taxon_b; verbose=false)

Print the shared trunk and the two lineage branches after their MRCA. Returns
the complete lineages and MRCA depth for programmatic use.
"""
function compare_lineages(taxon_a, taxon_b; verbose::Bool=false)
    lineage_a = get_lineage(taxon_a; verbose=verbose)
    lineage_b = get_lineage(taxon_b; verbose=verbose)

    if lineage_a === nothing || lineage_b === nothing
        println("Error: Could not retrieve one or both lineages")
        return nothing
    end

    result = _compute_distance(
        lineage_a,
        lineage_b;
        taxon_a=String(taxon_a),
        taxon_b=String(taxon_b),
    )
    mrca_depth = result.mrca_depth

    println("Lineage Comparison")
    println("MRCA: $(result.mrca) at depth $(mrca_depth)")

    if mrca_depth > 0
        println("\nShared lineage ($(mrca_depth) nodes):")
        for node in lineage_a[1:mrca_depth]
            println("  $(node)")
        end
    end

    println("\n$(taxon_a) only ($(length(lineage_a) - mrca_depth) nodes):")
    if length(lineage_a) > mrca_depth
        for node in lineage_a[(mrca_depth + 1):end]
            println("> $(node)")
        end
    end

    println("\n$(taxon_b) only ($(length(lineage_b) - mrca_depth) nodes):")
    if length(lineage_b) > mrca_depth
        for node in lineage_b[(mrca_depth + 1):end]
            println("> $(node)")
        end
    end

    return (
        lineage_a=String.(lineage_a),
        lineage_b=String.(lineage_b),
        mrca_depth=mrca_depth,
    )
end

"""
    shared_clades(taxon_a, taxon_b; verbose=false)

Return the continuous shared lineage from the root through the MRCA. Returns
an empty vector for disconnected hierarchies and `nothing` when either taxon
cannot be resolved.
"""
function shared_clades(taxon_a, taxon_b; verbose::Bool=false)
    lineage_a = get_lineage(taxon_a; verbose=verbose)
    lineage_b = get_lineage(taxon_b; verbose=verbose)
    (lineage_a === nothing || lineage_b === nothing) && return nothing

    result = _compute_distance(
        lineage_a,
        lineage_b;
        taxon_a=String(taxon_a),
        taxon_b=String(taxon_b),
    )
    result.mrca_depth == 0 && return String[]
    return String.(lineage_a[1:result.mrca_depth])
end

"""
    is_member(taxon, clade; verbose=false)

Test whether the complete clade name occurs in a taxon's lineage. Matching is
case-insensitive and ignores outer whitespace, but does not perform partial or
regular-expression matching.
"""
function is_member(taxon, clade; verbose::Bool=false)
    lineage = get_lineage(taxon; verbose=verbose)
    lineage === nothing && return nothing

    normalized_clade = lowercase(strip(String(clade)))
    normalized_lineage = lowercase.(strip.(String.(lineage)))
    return normalized_clade in normalized_lineage
end

"""
    filter_clade(taxa, clade; verbose=false)

Return only the taxa whose resolved lineage contains the requested clade.
Unresolved taxa are omitted.
"""
function filter_clade(
    taxa::AbstractVector,
    clade;
    verbose::Bool=false,
)
    kept = String[]
    for taxon in String.(taxa)
        membership = is_member(taxon, clade; verbose=verbose)
        membership === true && push!(kept, taxon)
    end
    return kept
end

"""
    taxo_path(taxon_a, taxon_b; verbose=false)

Return the node-by-node path ascending from taxon A to the MRCA and descending
to taxon B. The result contains `taxon_a`, `taxon_b`, and a `data` DataFrame
with columns `node`, `depth`, and `direction`.
"""
function taxo_path(taxon_a, taxon_b; verbose::Bool=false)
    name_a = String(taxon_a)
    name_b = String(taxon_b)
    lineage_a = get_lineage(name_a; verbose=verbose)
    lineage_b = get_lineage(name_b; verbose=verbose)

    if lineage_a === nothing
        @warn "Could not retrieve lineage for $(name_a)"
        return nothing
    end
    if lineage_b === nothing
        @warn "Could not retrieve lineage for $(name_b)"
        return nothing
    end

    result = _compute_distance(
        lineage_a,
        lineage_b;
        taxon_a=name_a,
        taxon_b=name_b,
    )
    mrca_depth = result.mrca_depth
    if mrca_depth == 0
        @warn "No common ancestor found for $(repr(name_a)) and $(repr(name_b))."
        return nothing
    end

    path = DataFrame(node=String[], depth=Int[], direction=String[])

    if mrca_depth < length(lineage_a)
        for index in length(lineage_a):-1:(mrca_depth + 1)
            push!(path, (
                node=String(lineage_a[index]),
                depth=index,
                direction="a",
            ))
        end
    end

    push!(path, (
        node=String(result.mrca),
        depth=mrca_depth,
        direction="mrca",
    ))

    if mrca_depth < length(lineage_b)
        for index in (mrca_depth + 1):length(lineage_b)
            push!(path, (
                node=String(lineage_b[index]),
                depth=index,
                direction="b",
            ))
        end
    end

    return (taxon_a=name_a, taxon_b=name_b, data=path)
end