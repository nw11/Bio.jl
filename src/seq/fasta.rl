

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
        yield = true;
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
    letters     = alpha+         >letters_start %letters_end;
    sequence    = (letters space+)*;
    fasta_entry = '>' identifier ( [ \t\v]+ description )? newline space* sequence space*;

    main := space* (fasta_entry %yield)*;
}%%


%% write data;


type FASTAParser
    state::Ragel.State
    seqbuf::Ragel.Buffer

    function FASTAParser(input::Union(IO, String))
        %% write init;

        return new(Ragel.State(cs, input), Ragel.Buffer())
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
        accept_state!(input, output)
    end)

end # module FASTAParserImpl


using Bio.Seq.FASTAParserImpl


