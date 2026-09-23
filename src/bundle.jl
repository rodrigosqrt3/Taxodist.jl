"""
    TaxodistBundle

A dictionary-like portable analysis bundle containing taxon resolution,
distances, metric definition, and provenance.
"""
mutable struct TaxodistBundle
    data::Dict{String,Any}
end

TaxodistBundle(; kwargs...) = TaxodistBundle(
    Dict{String,Any}(String(key) => value for (key, value) in kwargs),
)

function Base.getproperty(bundle::TaxodistBundle, name::Symbol)
    name === :data && return getfield(bundle, :data)
    data = getfield(bundle, :data)
    key = String(name)
    haskey(data, key) && return data[key]
    return getfield(bundle, name)
end

function Base.setproperty!(bundle::TaxodistBundle, name::Symbol, value)
    name === :data && return setfield!(bundle, :data, value)
    bundle.data[String(name)] = value
    return value
end

Base.getindex(bundle::TaxodistBundle, key::AbstractString) = bundle.data[String(key)]
Base.getindex(bundle::TaxodistBundle, key::Symbol) = bundle.data[String(key)]
Base.setindex!(bundle::TaxodistBundle, value, key::AbstractString) =
    setindex!(bundle.data, value, String(key))
Base.setindex!(bundle::TaxodistBundle, value, key::Symbol) =
    setindex!(bundle.data, value, String(key))
Base.haskey(bundle::TaxodistBundle, key) = haskey(bundle.data, String(key))
Base.keys(bundle::TaxodistBundle) = keys(bundle.data)
Base.delete!(bundle::TaxodistBundle, key) = delete!(bundle.data, String(key))

function Base.show(io::IO, ::MIME"text/plain", bundle::TaxodistBundle)
    validate_taxodist_bundle(bundle)
    counts = summary_counts(bundle.resolution)
    println(io, "TaxodistBundle schema $(bundle.schema_version)")
    println(io, "created: $(bundle.created_at)")
    println(io, "taxa: $(size(bundle.resolution, 1))")
    println(
        io,
        "resolved=$(counts.resolved), ambiguous=$(counts.ambiguous), " *
        "unresolved=$(counts.unresolved), retrieval_error=$(counts.retrieval_error)",
    )
end

"""
    taxo_bundle(taxa; ambiguity="warn", verbose=false, progress=true)

Combine an auditable resolution, its distance matrix, the shared metric
definition, and provenance in one portable object.
"""
function taxo_bundle(
    taxa;
    ambiguity="warn",
    verbose::Bool=false,
    progress::Bool=true,
)
    resolution = taxa isa TaxodistResolution ? taxa : taxo_resolve(
        taxa;
        ambiguity=ambiguity,
        verbose=verbose,
        progress=progress,
    )
    source_name = isempty(strip(resolution.source)) ? "The Taxonomicon" : resolution.source
    source_url = resolution.source_url
    if source_url === nothing && source_name == "The Taxonomicon"
        source_url = "http://taxonomicon.taxonomy.nl"
    end

    return TaxodistBundle(
        schema_version="1.0",
        created_at=_utc_timestamp(),
        source=Dict{String,Any}(
            "name" => source_name,
            "url" => source_url,
            "retrieved_at" => resolution.retrieved_at,
        ),
        software=Dict{String,Any}(
            "name" => "taxodist",
            "version" => string(TAXODIST_VERSION),
            "language" => "Julia",
        ),
        metric=Dict{String,Any}(
            "name" => "inverse_mrca_depth",
            "definition" => "0 for identical lineages; otherwise 1 / depth(MRCA)",
            "root_depth" => 1,
            "common_ancestor" => "continuous common lineage prefix",
        ),
        resolution=resolution,
        matrix=distance_matrix(resolution; progress=progress),
    )
end

_missing_value(value) = value === nothing || ismissing(value)

function _same_distance(left::Float64, right::Float64)
    isnan(left) && isnan(right) && return true
    left == right && return true
    isfinite(left) && isfinite(right) || return false
    return abs(left - right) <= sqrt(eps(Float64))
end

"""
    validate_taxodist_bundle(bundle)

Validate all scientific and structural invariants of a `TaxodistBundle`.
"""
function validate_taxodist_bundle(bundle)
    bundle isa TaxodistBundle || throw(ArgumentError(
        "bundle must be a TaxodistBundle object.",
    ))
    required = Set([
        "schema_version",
        "created_at",
        "source",
        "software",
        "metric",
        "resolution",
        "matrix",
    ])
    issubset(required, Set(keys(bundle))) || throw(ArgumentError(
        "Invalid TaxodistBundle: required fields are missing.",
    ))
    bundle.schema_version == "1.0" || throw(ArgumentError(
        "Unsupported taxodist bundle schema: $(repr(bundle.schema_version)).",
    ))

    resolution = bundle.resolution
    resolution isa TaxodistResolution || throw(ArgumentError(
        "Invalid bundle: resolution has the wrong type.",
    ))
    issubset(Set(RESOLUTION_COLUMNS), Set(names(resolution.data))) || throw(ArgumentError(
        "Invalid bundle: resolution fields are missing.",
    ))
    allowed = Set(["resolved", "ambiguous", "unresolved", "retrieval_error"])
    all(status -> status in allowed, resolution.status) || throw(ArgumentError(
        "Invalid bundle: unknown resolution status.",
    ))

    for (index, row) in enumerate(eachrow(resolution.data))
        candidates = row.candidates
        lineage = row.lineage
        count = row.n_candidates
        count isa Integer && count >= 0 || throw(ArgumentError(
            "Invalid bundle: candidate counts must be non-negative integers.",
        ))
        candidates isa DataFrame && issubset(Set(["id", "name"]), Set(names(candidates))) ||
            throw(ArgumentError(
                "Invalid bundle: malformed candidate table at row $(index).",
            ))
        nrow(candidates) == count || throw(ArgumentError(
            "Invalid bundle: candidate count mismatch at row $(index).",
        ))

        if lineage === nothing
            _missing_value(row.lineage_depth) || throw(ArgumentError(
                "Invalid bundle: lineage depth mismatch at row $(index).",
            ))
        else
            lineage isa AbstractVector && all(
                node -> node isa AbstractString && !isempty(strip(String(node))),
                lineage,
            ) || throw(ArgumentError(
                "Invalid bundle: malformed lineage at row $(index).",
            ))
            !_missing_value(row.lineage_depth) && length(lineage) == row.lineage_depth ||
                throw(ArgumentError(
                    "Invalid bundle: lineage depth mismatch at row $(index).",
                ))
        end

        if row.status in ("resolved", "ambiguous")
            (!_missing_value(row.id) && lineage !== nothing && nrow(candidates) > 0) ||
                throw(ArgumentError(
                    "Invalid bundle: incomplete resolved record at row $(index).",
                ))
            String(candidates.id[1]) == String(row.id) || throw(ArgumentError(
                "Invalid bundle: selected candidate mismatch at row $(index).",
            ))
        else
            (_missing_value(row.id) && lineage === nothing && _missing_value(row.resolved_name)) ||
                throw(ArgumentError(
                    "Invalid bundle: incomplete unresolved record at row $(index).",
                ))
            row.status == "unresolved" && nrow(candidates) != 0 && throw(ArgumentError(
                "Invalid bundle: unresolved record has candidates at row $(index).",
            ))
        end
    end

    matrix = bundle.matrix
    matrix isa TaxonomicDistanceMatrix || throw(ArgumentError(
        "Invalid bundle: matrix must be a TaxonomicDistanceMatrix.",
    ))
    matrix.taxa == String.(resolution.input) || throw(ArgumentError(
        "Invalid bundle: matrix labels do not match the resolution inputs.",
    ))
    expected = distance_matrix(resolution; progress=false)
    size(matrix) == size(expected) || throw(ArgumentError(
        "Invalid bundle: stored distances do not match the stored lineages.",
    ))
    all(
        _same_distance(matrix.values[index], expected.values[index])
        for index in eachindex(matrix.values)
    ) || throw(ArgumentError(
        "Invalid bundle: stored distances do not match the stored lineages.",
    ))
    return bundle
end

function distance_matrix(bundle::TaxodistBundle; kwargs...)
    validate_taxodist_bundle(bundle)
    return bundle.matrix
end

_analysis_matrix(bundle::TaxodistBundle; kwargs...) = distance_matrix(bundle)

function _portable_distance(value::Real)
    isnan(value) && return nothing
    isinf(value) && return value > 0 ? "Infinity" : "-Infinity"
    return Float64(value)
end

function _resolution_records(resolution::TaxodistResolution)
    records = Dict{String,Any}[]
    for row in eachrow(resolution.data)
        candidates = [
            Dict{String,Any}("id" => String(candidate.id), "name" => String(candidate.name))
            for candidate in eachrow(row.candidates)
        ]
        push!(records, Dict{String,Any}(
            "input" => String(row.input),
            "resolved_name" => _missing_value(row.resolved_name) ? nothing : String(row.resolved_name),
            "id" => _missing_value(row.id) ? nothing : String(row.id),
            "status" => String(row.status),
            "n_candidates" => Int(row.n_candidates),
            "lineage_depth" => _missing_value(row.lineage_depth) ? nothing : Int(row.lineage_depth),
            "lineage" => row.lineage === nothing ? nothing : String.(row.lineage),
            "candidates" => candidates,
        ))
    end
    return records
end

"""
    write_taxodist_bundle(bundle, file; pretty=true)

Write a validated bundle using the shared language-independent JSON schema
version 1.0. Returns the absolute output path.
"""
function write_taxodist_bundle(
    bundle::TaxodistBundle,
    file::AbstractString;
    pretty::Bool=true,
)
    validate_taxodist_bundle(bundle)
    matrix = bundle.matrix
    portable = Dict{String,Any}(
        "format" => "taxodist_bundle",
        "schema_version" => bundle.schema_version,
        "created_at" => bundle.created_at,
        "source" => bundle.source,
        "software" => bundle.software,
        "metric" => bundle.metric,
        "taxa" => _resolution_records(bundle.resolution),
        "matrix" => Dict{String,Any}(
            "labels" => copy(matrix.taxa),
            "values" => [
                [_portable_distance(matrix.values[row, column]) for column in axes(matrix.values, 2)]
                for row in axes(matrix.values, 1)
            ],
        ),
    )
    open(file, "w") do stream
        if pretty
            JSON3.pretty(stream, portable)
        else
            JSON3.write(stream, portable)
        end
        write(stream, '\n')
    end
    return abspath(file)
end

_bundle_dict(value) = Dict{String,Any}(String(key) => item for (key, item) in pairs(value))

function _required_bundle_keys(object, required, context)
    available = Set(String(key) for key in keys(object))
    issubset(Set(required), available) || throw(ArgumentError(
        "Invalid bundle JSON: required $(context) fields are missing.",
    ))
end

function _matrix_from_bundle_json(matrix_data)
    object = _bundle_dict(matrix_data)
    _required_bundle_keys(object, ["labels", "values"], "matrix")
    labels = String.(collect(object["labels"]))
    rows = collect(object["values"])
    length(rows) == length(labels) || throw(ArgumentError(
        "Invalid bundle JSON: matrix row count does not match labels.",
    ))
    values = Matrix{Float64}(undef, length(labels), length(labels))
    for row_index in eachindex(rows)
        row = collect(rows[row_index])
        length(row) == length(labels) || throw(ArgumentError(
            "Invalid bundle JSON: matrix must be square.",
        ))
        for column_index in eachindex(row)
            value = row[column_index]
            values[row_index, column_index] = if value === nothing
                NaN
            elseif value == "Infinity"
                Inf
            elseif value == "-Infinity"
                -Inf
            elseif value isa Real
                Float64(value)
            else
                throw(ArgumentError("Invalid bundle JSON: invalid matrix value."))
            end
        end
    end
    all(index -> values[index, index] == 0.0, eachindex(labels)) || throw(ArgumentError(
        "Invalid bundle JSON: matrix diagonal must contain zeros.",
    ))
    all(
        _same_distance(values[row, column], values[column, row])
        for row in axes(values, 1), column in axes(values, 2)
    ) || throw(ArgumentError(
        "Invalid bundle JSON: distance matrix must be symmetric.",
    ))
    return TaxonomicDistanceMatrix(values, labels)
end

function _resolution_from_bundle_json(records, source)
    rows = NamedTuple[]
    for raw_record in records
        record = _bundle_dict(raw_record)
        _required_bundle_keys(
            record,
            RESOLUTION_COLUMNS,
            "taxon",
        )
        candidate_table = _empty_candidates()
        for raw_candidate in record["candidates"]
            candidate = _bundle_dict(raw_candidate)
            _required_bundle_keys(candidate, ["id", "name"], "candidate")
            push!(candidate_table, (
                id=String(candidate["id"]),
                name=String(candidate["name"]),
            ))
        end
        lineage = record["lineage"] === nothing ? nothing : String.(collect(record["lineage"]))
        push!(rows, (
            input=String(record["input"]),
            resolved_name=record["resolved_name"] === nothing ? nothing : String(record["resolved_name"]),
            id=record["id"] === nothing ? nothing : String(record["id"]),
            status=String(record["status"]),
            n_candidates=Int(record["n_candidates"]),
            lineage_depth=record["lineage_depth"] === nothing ? nothing : Int(record["lineage_depth"]),
            lineage=lineage,
            candidates=candidate_table,
        ))
    end
    return _make_resolution(
        rows;
        source=String(source["name"]),
        source_url=source["url"] === nothing ? nothing : String(source["url"]),
        retrieved_at=source["retrieved_at"] === nothing ? nothing : String(source["retrieved_at"]),
    )
end

"""
    read_taxodist_bundle(file)

Read and validate a portable bundle written by R, Python, or Julia.
"""
function read_taxodist_bundle(file::AbstractString)
    isfile(file) || throw(ArgumentError("Bundle file not found: $(file)"))
    raw = try
        JSON3.read(read(file, String))
    catch error
        throw(ArgumentError("Could not parse bundle JSON: $(sprint(showerror, error))"))
    end
    raw isa JSON3.Object || throw(ArgumentError(
        "Invalid bundle JSON: expected an object.",
    ))
    object = _bundle_dict(raw)
    get(object, "format", nothing) == "taxodist_bundle" || throw(ArgumentError(
        "Invalid bundle JSON: unrecognized format.",
    ))
    schema_version = get(object, "schema_version", nothing)
    schema_version == "1.0" || throw(ArgumentError(
        "Unsupported taxodist bundle schema: $(repr(schema_version)).",
    ))
    _required_bundle_keys(
        object,
        ["created_at", "source", "software", "metric", "taxa", "matrix"],
        "top-level",
    )

    source = _bundle_dict(object["source"])
    _required_bundle_keys(source, ["name", "url", "retrieved_at"], "source")
    resolution = _resolution_from_bundle_json(object["taxa"], source)
    bundle = TaxodistBundle(
        schema_version=String(schema_version),
        created_at=object["created_at"] === nothing ? nothing : String(object["created_at"]),
        source=source,
        software=_bundle_dict(object["software"]),
        metric=_bundle_dict(object["metric"]),
        resolution=resolution,
        matrix=_matrix_from_bundle_json(object["matrix"]),
    )
    return validate_taxodist_bundle(bundle)
end