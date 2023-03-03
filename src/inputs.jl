
@enum GenVarInfo begin
    chr
    pos
    rsid
    chr_pos
    chr_colong_pos
    exposition_name
    exposition_isoform
    exposition_name_isoform
    β
    se
    odd_ratio
    # Add others ???
    other
end


struct Expositions
    path::String
    to_file::Vector{Any}
    to_file_regex::Regex
end

function _make_regex(x::Vector{Any})::Regex
    str = "";
    for i in x
        if (i isa String); str*=i; else; str*="(.)"; end;     # pourra etre augmenté pour prendre en compte les types d'infos et les conditions a respecter.
    end
    return Regex(str)
end

Expositions(path::String, to_file::Vector{Any}) = Expositions(path, to_file, _make_regex(to_file))


struct GWAS
    paths::Union{String, Expositions}
    columns::Vector{GenVarInfo} # regex ??
    separator::Union{String, Char}
end
