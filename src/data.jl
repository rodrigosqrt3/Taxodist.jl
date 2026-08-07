_json_dict(value) = Dict(String(key) => item for (key, item) in pairs(value))

function _matrix_from_json(value, fallback_labels)
    labels = String.(fallback_labels)
    rows = value

    if value isa JSON3.Object
        object = _json_dict(value)
        labels = String.(collect(get(object, "labels", fallback_labels)))
        rows = get(object, "data", Any[])
    end

    n_labels = length(labels)
    matrix = Matrix{Float64}(undef, n_labels, n_labels)
    collected_rows = collect(rows)
    length(collected_rows) == n_labels || throw(ArgumentError(
        "Invalid taxobase matrix: row count does not match its labels.",
    ))

    for row_index in 1:n_labels
        row = collect(collected_rows[row_index])
        length(row) == n_labels || throw(ArgumentError(
            "Invalid taxobase matrix: column count does not match its labels.",
        ))
        matrix[row_index, :] = Float64.(row)
    end

    return TaxonomicDistanceMatrix(matrix, labels)
end

function _table_from_records(records, columns)
    table = DataFrame()
    for column in columns
        table[!, column] = Any[]
    end

    for record in records
        values = _json_dict(record)
        push!(table, (; (column => get(values, String(column), nothing) for column in columns)...))
    end
    return table
end

function _coverage_table(taxa, coverage_raw)
    coverage_values = if coverage_raw isa AbstractVector
        Dict(taxon => value for (taxon, value) in zip(taxa, coverage_raw))
    else
        _json_dict(coverage_raw)
    end
    return DataFrame(
        taxon=copy(taxa),
        covered=[Bool(get(coverage_values, taxon, false)) for taxon in taxa],
    )
end

"""
    load_taxobase()

Load the packaged offline reference dataset exported from the R `taxobase`
object. The result contains the same taxa, matrices, examples, search output,
and provenance metadata used by the R and Python implementations.
"""
function load_taxobase()
    path = joinpath(@__DIR__, "data", "taxobase.json")
    isfile(path) || throw(ArgumentError(
        "Packaged taxobase file not found at $(repr(path)).",
    ))

    raw = JSON3.read(read(path, String))
    raw isa JSON3.Object || throw(ArgumentError(
        "Invalid taxobase file: expected a JSON object.",
    ))
    data = _json_dict(raw)

    taxa = String.(collect(get(data, "taxa", Any[])))
    found_taxa = String.(collect(get(data, "found_taxa", taxa)))

    coverage_raw = get(data, "coverage", Dict{String,Any}())
    coverage = _coverage_table(taxa, coverage_raw)

    matrix = _matrix_from_json(get(data, "matrix", Any[]), found_taxa)

    pairwise_values = _json_dict(get(
        data,
        "pairwise",
        Dict{String,Any}(),
    ))
    pairwise = (
        distance=Float64(get(pairwise_values, "distance", NaN)),
        mrca=String(get(pairwise_values, "mrca", "")),
        mrca_depth=Int(get(pairwise_values, "mrca_depth", 0)),
        depth_a=Int(get(pairwise_values, "depth_a", 0)),
        depth_b=Int(get(pairwise_values, "depth_b", 0)),
        taxon_a=String(get(pairwise_values, "taxon_a", "")),
        taxon_b=String(get(pairwise_values, "taxon_b", "")),
    )

    lineage_homo = String.(collect(get(data, "lineage_homo", Any[])))
    lineage_tyrannosaurus = String.(collect(
        get(data, "lineage_tyrannosaurus", Any[]),
    ))

    closest = _table_from_records(
        get(data, "closest", Any[]),
        (:taxon, :distance),
    )
    closest.taxon = String.(closest.taxon)
    closest.distance = Float64.(closest.distance)

    filtered = String.(collect(get(data, "filter", Any[])))

    search = _table_from_records(
        get(data, "search", Any[]),
        (:id, :name),
    )
    search.id = String.(search.id)
    search.name = String.(search.name)

    statistical_taxa = String.(collect(get(data, "statistical_taxa", Any[])))
    statistical_matrix = _matrix_from_json(
        get(data, "statistical_matrix", Any[]),
        statistical_taxa,
    )

    metadata_values = _json_dict(get(
        data,
        "metadata",
        Dict{String,Any}(),
    ))
    metadata = (
        source=String(get(metadata_values, "source", "The Taxonomicon")),
        source_url=String(get(
            metadata_values,
            "source_url",
            "http://taxonomicon.taxonomy.nl",
        )),
        generated_on=String(get(metadata_values, "generated_on", "")),
        package_version=String(get(metadata_values, "package_version", "legacy")),
        distance_definition=String(get(
            metadata_values,
            "distance_definition",
            "legacy packaged data",
        )),
    )

    return (
        taxa=taxa,
        found_taxa=found_taxa,
        coverage=coverage,
        matrix=matrix,
        pairwise=pairwise,
        lineage_homo=lineage_homo,
        lineage_tyrannosaurus=lineage_tyrannosaurus,
        closest=closest,
        filter=filtered,
        search=search,
        statistical_taxa=statistical_taxa,
        statistical_matrix=statistical_matrix,
        metadata=metadata,
    )
end