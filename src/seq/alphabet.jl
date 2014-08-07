

# Assign flags for sequence types so we can maintain a set of compatible
# alphabets in an integer.
const EMPTY_ALPHABET = uint16(0)
const DNA_ALPHABET   = uint16(0b0001)
const RNA_ALPHABET   = uint16(0b0010)
const AA_ALPHABET    = uint16(0b0100)

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
#
# Returns:
#   A type T to which the string data can be converted.
#
function infer_sequence_type(data::String)
    alphabets = ALL_ALPHABETS
    for c in data
        if 'A' <= c <= 'z'
            alphabets &= compatible_alphabets[c - 'A' + 1]
        else
            error("Character $(c) is not compatible with any sequence type.")
        end
    end

    if count_ones(alphabets) == 0
        error("String is not compatible with any known sequence type.")
    else
        for alphabet in preferred_sequence_alphabets
            if alphabet & alphabets != 0
                return alphabet_type[alphabet]
            end
        end
    end
end


