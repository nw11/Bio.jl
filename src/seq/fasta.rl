

# FASTA sequence types

type FASTAMetadata
    description::String

    function FASTAMetadata()
        return new("")
    end
end


typealias FASTASeqRecord{S}       SeqRecord{S, FASTAMetadata}
typealias FASTADNASeqRecord       DNASeqRecord{FASTAMetadata}
typealias FASTARNASeqRecord       RNASeqRecord{FASTAMetadata}
typealias FASTAAminoAcidSeqRecord AminoAcidSeqRecord{FASTAMetadata}


module FASTAParserImpl

import Bio.Seq: FASTASeqRecord
import Bio.Ragel
using Switch
export FASTAParser


%%{
    machine fasta;

    action yield {
        fbreak;
    }

    action identifier_start {
        Ragel.@pushmark!
    }

    action identifier_end {
        firstpos = Ragel.@popmark!
        output.name = bytestring(Ragel.@spanfrom firstpos)
    }

    action description_start {
        Ragel.@pushmark!
    }

    action description_end {
        firstpos = Ragel.@popmark!
        output.metadata.description = bytestring(Ragel.@spanfrom firstpos)
    }

    action letter {
        write(input.seqbuf, Ragel.@char)
    }

    identifier  = (any - space)+ >identifier_start %identifier_end;
    description = [^\r\n]+       >description_start %description_end;
    newline     = [\n\r];
    letter      = (alpha >letter) space*;
    sequence    = letter*;
    fasta_entry = '>' identifier ( [ \t\v]+ description )? newline space* sequence;

    main := space* (fasta_entry %yield)*;
}%%


%% write data;


type FASTAParser
    state::Ragel.State
    seqbuf::IOBuffer

    function FASTAParser(filename::String)
        %% write init;

        return new(Ragel.State(cs, filename), IOBuffer())
    end
end


function Ragel.ragelstate(parser::FASTAParser)
    return parser.state
end


function accept_state!{S}(input::FASTAParser, output::FASTASeqRecord{S})
    seqstr = takebuf_string(input.seqbuf)
    output.seq = S(seqstr)
end


Ragel.@generate_read_fuction("fasta", FASTAParser, FASTASeqRecord,
    begin
        %% write exec;
    end,
    begin
        @show p
        accept_state!(input, output)
    end)

end # module FASTAParserImpl


using Bio.Seq.FASTAParserImpl

