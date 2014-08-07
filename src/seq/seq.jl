
module Seq

using Base.Intrinsics

import Base: convert, getindex, show, length, start, next, done, copy, reverse,
             show, endof, ==

export Nucleotide, DNANucleotide, RNANucleotide,
       DNA_A, DNA_C, DNA_G, DNA_T, DNA_N,
       RNA_A, RNA_C, RNA_G, RNA_U, RNA_N,
       NucleotideSequence, DNASequence, RNASequence, @dna_str, @rna_str,
       complement, reverse_complement, mismatches, npositions, hasn, eachsubseq,
       eachkmer, nucleotide_count, Kmer, DNAKmer, RNAKmer, dnakmer, rnakmer,
       kmer, AminoAcid, AminoAcidSequence, @aa_str, translate, ncbi_trans_table,
       AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H, AA_I, AA_L,
       AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W, AA_Y, AA_V, AA_X


abstract Sequence


include("nucleotide.jl")
include("kmer.jl")
include("aminoacid.jl")
include("alphabet.jl")


# A sequence record is a named sequence, optionally with attached metadat.
immutable SeqRecord{S, T}
    name::String
    seq::S
    metadata::T
end


# Parsing of various file types
include("fasta.jl")

end # module Seq

