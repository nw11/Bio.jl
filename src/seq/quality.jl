
# A representation of positions-specific integer quality scores, as in FASTQ.


bitstype 16 QualityEncoding

function Base.convert(::Type{QualityEncoding}, nt::Uint16)
    return box(QualityEncoding, unbox(Uint16, nt))
end


function Base.convert(::Type{Uint16}, nt::QualityEncoding)
    return box(Uint16, unbox(QualityEncoding, nt))
end


function Base.|(a::QualityEncoding, b::QualityEncoding)
    return convert(QualityEncoding, convert(Uint16, a) | convert(Uint16, b))
end


function Base.&(a::QualityEncoding, b::QualityEncoding)
    return convert(QualityEncoding, convert(Uint16, a) & convert(Uint16, b))
end


const EMPTY_QUAL_ENCODING      = convert(QualityEncoding, uint16(0))
const SANGER_QUAL_ENCODING     = convert(QualityEncoding, uint16(0b00001))
const SOLEXA_QUAL_ENCODING     = convert(QualityEncoding, uint16(0b00010))
const ILLUMINA13_QUAL_ENCODING = convert(QualityEncoding, uint16(0b00100))
const ILLUMINA15_QUAL_ENCODING = convert(QualityEncoding, uint16(0b01000))
const ILLUMINA18_QUAL_ENCODING = convert(QualityEncoding, uint16(0b10000))
const ALL_QUAL_ENCODINGS =
    SANGER_QUAL_ENCODING | SOLEXA_QUAL_ENCODING | ILLUMINA13_QUAL_ENCODING |
    ILLUMINA15_QUAL_ENCODING | ILLUMINA18_QUAL_ENCODING


# Ranges and score of the first character in the range.
const qual_encoding_ranges = @compat Dict{QualityEncoding, (typeof(uint8(0):uint8(0)), Int8)}(
    SANGER_QUAL_ENCODING     => (uint8('!'):uint8('~'), int8(0)),
    SOLEXA_QUAL_ENCODING     => (uint8(';'):uint8('~'), int8(-5)),
    ILLUMINA13_QUAL_ENCODING => (uint8('@'):uint8('~'), int8(0)),
    ILLUMINA15_QUAL_ENCODING => (uint8('B'):uint8('~'), int8(3)),
    ILLUMINA18_QUAL_ENCODING => (uint8('!'):uint8('~'), int8(0)),
)

# Build an encoding lookup table
const compatible_qual_encoding = fill(EMPTY_QUAL_ENCODING, length('!':'~'))

for (encoding, (range, start_score)) in qual_encoding_ranges
    for c in range
        compatible_qual_encoding[c - uint8('!') + 1] |= encoding
    end
end

# When a quality string has multiple compatible encodings, we choose the first
# compatible alphabet in this list.
const preferred_quality_encodings = [
    ILLUMINA15_QUAL_ENCODING,
    ILLUMINA13_QUAL_ENCODING,
    SOLEXA_QUAL_ENCODING,
    SANGER_QUAL_ENCODING,
    ILLUMINA18_QUAL_ENCODING,
]

# Infer the encoding of a FASTQ quality string.
#
# Args:
#   data: sequence data in a string
#   start: first position to consider in data
#   stop: last position to consider in data
#   default: if there are multiple compatible encodings, default
#             to this one if it's compatible.
#
# Returns:
#   A quality encoding.
#
function infer_quality_encoding(data::Vector{Uint8}, start, stop, default)
    encodings = ALL_QUAL_ENCODINGS
    for i in start:stop
        c = data[i]
        if '!' <= c <= '~'
            encodings &= compatible_qual_encoding[convert(Char, c) - '!' + 1]
        else
            error("Character $(convert(Char, c)) is not compatible with any known quality encoding.")
        end
    end

    if count_ones(convert(Uint16, encodings)) == 0
        error("String is not compatible with any known sequence type.")
    elseif default != EMPTY_QUAL_ENCODING && (encodings & default) != EMPTY_ALPHABET
        return default
    else
        for encoding in preferred_quality_encodings
            if encoding & encodings != EMPTY_QUAL_ENCODING
                return encoding
            end
        end
    end
end


function decode_quality_string!(encoding::QualityEncoding, input::Vector{Uint8},
                                output::Vector{Int8}, start, stop)
    range, startqual = qual_encoding_ranges[encoding]
    for (i, j) in enumerate(start:stop)
        c = input[j]
        output[i] = startqual + (c - range.start)
    end
    return output
end


function decode_quality_string(encoding::QualityEncoding, input::Vector{Uint8},
                               start, stop)
    output = Array(Int8, stop - start + 1)
    return decode_quality_string!(encoding, input, output, start, stop)
end


