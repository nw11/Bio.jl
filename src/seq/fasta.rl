

# FASTA sequence types

type FASTAMetadata
    description::String

    function FASTAMetadata(description)
        return new(description)
    end

    function FASTAMetadata()
        return new("")
    end
end


typealias FASTASeqRecord{S}       SeqRecord{S, FASTAMetadata}
typealias FASTADNASeqRecord       DNASeqRecord{FASTAMetadata}
typealias FASTARNASeqRecord       RNASeqRecord{FASTAMetadata}
typealias FASTAAminoAcidSeqRecord AminoAcidSeqRecord{FASTAMetadata}


function Base.show(io::IO, seqrec::FASTASeqRecord)
    write(io, ">", seqrec.name, " ", seqrec.metadata.description, "\n")
    show(io, seqrec.seq)
end


module FASTAParserImpl

import Bio.Seq: FASTASeqRecord
import Bio.Ragel
using Switch
export FASTAParser


%%{
    machine fasta;

    action yield {
        yield = true;
        fbreak;
    }

    action identifier_start {
        Ragel.@pushmark!
    }

    action identifier_end {
        firstpos = Ragel.@popmark!
        input.namebuf = bytestring(Ragel.@spanfrom firstpos)
    }

    action description_start {
        Ragel.@pushmark!
    }

    action description_end {
        firstpos = Ragel.@popmark!
        input.descbuf = bytestring(Ragel.@spanfrom firstpos)
    }

    action letters_start {
        Ragel.@pushmark!
    }

    action letters_end {
        firstpos = Ragel.@popmark!
        append!(input.seqbuf, state.buffer, firstpos, p)
    }

    identifier  = (any - space)+ >identifier_start %identifier_end;
    description = [^\r\n]+       >description_start %description_end;
    newline     = [\n\r];
    hspace      = [ \t\v];
    letters     = alpha+         >letters_start %letters_end;
    sequence    = space* letters? (newline+ space* letters (hspace+ letters)*)*;
    fasta_entry = '>' identifier ( [ \t\v]+ description )? newline sequence space*;

    main := space* (fasta_entry %yield)*;
}%%


%% write data;


type FASTAParser
    state::Ragel.State
    seqbuf::Ragel.Buffer
    namebuf::String
    descbuf::String

    function FASTAParser(input::Union(IO, String, Vector{Uint8}))
        %% write init;

        return new(Ragel.State(cs, input), Ragel.Buffer{Uint8}(), "", "")
    end
end


function Ragel.ragelstate(parser::FASTAParser)
    return parser.state
end


function accept_state!{S}(input::FASTAParser, output::FASTASeqRecord{S})
    output.name = input.namebuf
    output.metadata.description = input.descbuf
    output.seq = S(input.seqbuf.data, 1, input.seqbuf.pos - 1)

    input.namebuf = ""
    input.descbuf = ""
    empty!(input.seqbuf)
end


Ragel.@generate_read_fuction("fasta", FASTAParser, FASTASeqRecord,
    begin
        @inbounds begin
            %% write exec;
        end
    end,
    begin
        accept_state!(input, output)
    end)

end # module FASTAParserImpl


using Bio.Seq.FASTAParserImpl


type FASTAIterator
    parser::FASTAParser

    # A type or function used to construct output sequence types
    default_alphabet::Alphabet
    isdone::Bool
    nextitem
end


function Base.read(input::IO, ::Type{FASTA}, alphabet::Alphabet=DNA_ALPHABET)
    return FASTAIterator(FASTAParser(input), alphabet, false, nothing)
end


function advance!(it::FASTAIterator)
    it.isdone = !FASTAParserImpl.advance!(it.parser)
    if !it.isdone
        alphabet = infer_alphabet(it.parser.seqbuf.data, 1, it.parser.seqbuf.pos - 1,
                                  it.default_alphabet)
        S = alphabet_type[alphabet]
        it.default_alphabet = alphabet
        it.nextitem =
            FASTASeqRecord{S}(it.parser.namebuf,
                              S(it.parser.seqbuf.data, 1, it.parser.seqbuf.pos - 1),
                              FASTAMetadata(it.parser.descbuf))
        empty!(it.parser.seqbuf)
    end
end


function start(it::FASTAIterator)
    advance!(it)
    return nothing
end


function next(it::FASTAIterator, state)
    item = it.nextitem
    advance!(it)
    return item, nothing
end


function done(it::FASTAIterator, state)
    return it.isdone
end

