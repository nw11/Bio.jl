
include("../ragel.jl")

module FastaParserImpl

export FastaParser

import Ragel
using Switch

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
        println(bytestring(Ragel.@spanfrom firstpos))
    }

    action description {
        # TODO
    }

    action sequence {
        # TODO
    }

    identifier  = (any - space)+ >identifier_start %identifier_end;
    description = [^\r\n]+       %description;
    newline     = [\n\r];
    sequence    = space* alpha   @sequence;
    fasta_entry = '>' identifier ( [ \t\v]+ description )? newline sequence* > yield;

    main := (space* fasta_entry)* space* $/yield;
}%%


%% write data;


type FastaParser
    state::Ragel.State

    function FastaParser(filename::String)
        %% write init;

        return new(Ragel.State(cs, filename))
    end
end


function Ragel.ragelstate(parser::FastaParser)
    return parser.state
end


type DummyFastaType
end


Ragel.@generate_read_fuction("fasta", FastaParser, DummyFastaType,
    begin
        %% write exec;
    end,
    begin
    end)

end

