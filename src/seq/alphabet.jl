

# Assign flags for sequence types so we can maintain a set of compatible
# alphabets in an integer.

bitstype 16 Alphabet

function Base.convert(::Type{Alphabet}, nt::Uint16)
    return box(Alphabet, unbox(Uint16, nt))
end


function Base.convert(::Type{Uint16}, nt::Alphabet)
    return box(Uint16, unbox(Alphabet, nt))
end


function Base.|(a::Alphabet, b::Alphabet)
    return convert(Alphabet, convert(Uint16, a) | convert(Uint16, b))
end


function Base.&(a::Alphabet, b::Alphabet)
    return convert(Alphabet, convert(Uint16, a) & convert(Uint16, b))
end




const EMPTY_ALPHABET = convert(Alphabet, uint16(0))
const DNA_ALPHABET   = convert(Alphabet, uint16(0b0001))
const RNA_ALPHABET   = convert(Alphabet, uint16(0b0010))
const AA_ALPHABET    = convert(Alphabet, uint16(0b0100))

const ALL_ALPHABETS =
    DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET


const alphabet_type = [
    DNA_ALPHABET => DNASequence,
    RNA_ALPHABET => RNASequence,
    AA_ALPHABET  => AminoAcidSequence
]


# When a sequence has multiple compatible alphabets, we choose the first
# compatible alphabet in this list.
const preferred_sequence_alphabets = [
    DNA_ALPHABET, RNA_ALPHABET, AA_ALPHABET
]


# Lookup table mapping a character in 'A':'z' to an integer representing the set
# of alphabets that character is compatible with.
const compatible_alphabets = [
    # A                                          B
      DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET, EMPTY_ALPHABET,
    # C                                          D
      DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET, AA_ALPHABET,
    # E            F            G
      AA_ALPHABET, AA_ALPHABET, DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET,
    # H            I            J               K            L
      AA_ALPHABET, AA_ALPHABET, EMPTY_ALPHABET, AA_ALPHABET, AA_ALPHABET,
    # M            N                                          O
      AA_ALPHABET, DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET, EMPTY_ALPHABET,
    # P            Q            R            S
      AA_ALPHABET, AA_ALPHABET, AA_ALPHABET, AA_ALPHABET,
    # T                           U             V            W
      DNA_ALPHABET | AA_ALPHABET, RNA_ALPHABET, AA_ALPHABET, AA_ALPHABET,
    # X            Y            Z
      AA_ALPHABET, AA_ALPHABET, EMPTY_ALPHABET,

    # [               \               ]               ^
      EMPTY_ALPHABET, EMPTY_ALPHABET, EMPTY_ALPHABET, EMPTY_ALPHABET,
    #  _              `
      EMPTY_ALPHABET, EMPTY_ALPHABET,

    # a                                          b
      DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET, EMPTY_ALPHABET,
    # c                                          d
      DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET, AA_ALPHABET,
    # e            f            g
      AA_ALPHABET, AA_ALPHABET, DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET,
    # h            i            j               k            l
      AA_ALPHABET, AA_ALPHABET, EMPTY_ALPHABET, AA_ALPHABET, AA_ALPHABET,
    # m            n                                          o
      AA_ALPHABET, DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET, EMPTY_ALPHABET,
    # p            q            r            s
      AA_ALPHABET, AA_ALPHABET, AA_ALPHABET, AA_ALPHABET,
    # t                           u             v            w
      DNA_ALPHABET | AA_ALPHABET, RNA_ALPHABET, AA_ALPHABET, AA_ALPHABET,
    # x            y            z
      AA_ALPHABET, AA_ALPHABET, EMPTY_ALPHABET
]


# Infer the sequence type by inspecting a string.
#
# Args:
#   data: sequence data in a string
#   start: first position to consider in data
#   stop: last position to consider in data
#   default: if there are multiple compatible alphabets, default
#             to this one if it's compatible.
#
# Returns:
#   A type T to which the string data can be converted.
#
function infer_alphabet(data::Vector{Uint8}, start, stop, default)
    alphabets = ALL_ALPHABETS
    for i in start:stop
        c = data[i]
        if 'A' <= c <= 'z'
            alphabets &= compatible_alphabets[c - 'A' + 1]
        else
            error("Character $(c) is not compatible with any sequence type.")
        end
    end

    if count_ones(convert(Uint16, alphabets)) == 0
        error("String is not compatible with any known sequence type.")
    elseif alphabets & default != 0
        return default
    else
        for alphabet in preferred_sequence_alphabets
            if alphabet & alphabets != 0
                return alphabet
            end
        end
    end
end


