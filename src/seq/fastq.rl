

# FASTQ sequence types

type FASTQMetadata
    description::String
    quality::Vector{Int8}

    function FASTQMetadata(description, quality)
        return new(description, quality)
    end

    function FASTQ()
        return new("", Int8[])
    end
end



typealias FASTQSeqRecord DNASeqRecord{FASTQMetadata}


function Base.show(io::IO, seqrec::FASTQSeqRecord)
    write(io, "@", seqrec.name, " ", seqrec.metadata.description, "\n")
    for c in seqrec.seq
        show(io, c)
    end
    write(io, '\n')
    # print quality scores as a unicode bar chart
    for q in seqrec.metadata.quality
        if q <= 0
            write(io, '▁')
        elseif q <= 6
            write(io, '▂')
        elseif q <= 12
            write(io, '▃')
        elseif q <= 18
            write(io, '▄')
        elseif q <= 24
            write(io, '▅')
        elseif q <= 30
            write(io, '▆')
        elseif q <= 36
            write(io, '▇')
        else
            write(io, '█')
        end
    end
    write(io, '\n')
end


module FASTQParserImpl

import Bio.Seq: FASTQSeqRecord, QualityEncoding, ILLUMINA18_QUAL_ENCODING
import Bio.Ragel
using Switch
export FASTQParser


%%{
    machine fastq;

    action yield {
        println(STDERR, "yield")
        yield = true;
        fbreak;
    }

    action identifier_start {
        println(STDERR, "identifier_start")
        Ragel.@pushmark!
    }

    action identifier_end {
        println(STDERR, "identifier_end")
        firstpos = Ragel.@popmark!
        input.namebuf = bytestring(Ragel.@spanfrom firstpos)
    }

    action description_start {
        println(STDERR, "description_start")
        Ragel.@pushmark!
    }

    action description_end {
        println(STDERR, "description_end")
        firstpos = Ragel.@popmark!
        input.descbuf = bytestring(Ragel.@spanfrom firstpos)
    }

    action letters_start {
        println(STDERR, "letters_start")
        Ragel.@pushmark!
    }

    action letters_end {
        println(STDERR, "letters_end")
        firstpos = Ragel.@popmark!
        append!(input.seqbuf, state.buffer, firstpos, p)
    }

    action qletters_start {
        println(STDERR, "qletters_start")
        Ragel.@pushmark!
    }

    action qletters_end {
        println(STDERR, "qletters_end")
        firstpos = Ragel.@popmark!
        append!(input.qualbuf, state.buffer, firstpos, p)
    }

    action test_qlen {
        length(input.qualbuf) < length(input.seqbuf)
    }

    newline     = [\n\r];
    hspace      = [ \t\v];
    identifier  = (any - space)+ >identifier_start %identifier_end;
    description = [^\r\n]+       >description_start %description_end;

    letters     = alpha+         >letters_start %letters_end;
    sequence    = (newline+ hspace* letters (hspace+ letters)*)*;

    qletters    = ([!-h] when test_qlen)+         >qletters_start %qletters_end;
    quality     = (newline+ hspace* qletters (hspace+ letters)*)*;

    fastq_entry = '@' identifier ( [ \t\v]+ description )?
                  sequence
                  newline '+' hspace* # TODO: repeated identifiers
                  quality;

    main := space* (fastq_entry %yield)?;
}%%


%% write data;


type FASTQParser
    state::Ragel.State
    seqbuf::Ragel.Buffer
    qualbuf::Ragel.Buffer
    namebuf::String
    descbuf::String
    default_qual_encoding::QualityEncoding

    function FASTQParser(input::Union(IO, String, Vector{Uint8}),
                         default_qual_encoding=ILLUMINA18_QUAL_ENCODING)
        %% write init;

        return new(Ragel.State(cs, input), Ragel.Buffer{Uint8}(),
                   Ragel.Buffer{Uint8}(), "", "", default_qual_encoding)
    end
end


function Ragel.ragelstate(parser::FASTQParser)
    return parser.state
end


function accept_state!(input::FASTQParser, output::FASTQSeqRecord)
    output.name = input.namebuf
    output.metadata.description = input.descbuf
    output.seq = DNASequence(input.seqbuf.data, 1, input.seqbuf.pos - 1)

    encoding = infer_quality_encoding(input.qualbuf.data, 1,
                                      input.qualbuf.pos - 1,
                                      input.default_qual_encoding)
    input.default_qual_encoding = encoding
    output.metadata.quality = decode_quality_string(encoding, input.qualbuf.
                                                    1, input.qualbuf.pos - 1)

    input.namebuf = ""
    input.descbuf = ""
    empty!(input.seqbuf)
    empty!(input.qualbuf)
end


Ragel.@generate_read_fuction("fastq", FASTQParser, FASTQSeqRecord,
    begin
        @inbounds begin
            %% write exec;
        end
    end,
    begin
        accept_state!(input, output)
    end)

end # module FASTQParserImpl


using Bio.Seq.FASTQParserImpl


type FASTQIterator
    parser::FASTQParser

    # A type or function used to construct output sequence types
    default_qual_encoding::QualityEncoding
    isdone::Bool
    nextitem
end


function Base.read(input::IO, ::Type{FASTQ},
                   qual_encoding::QualityEncoding=ILLUMINA18_QUAL_ENCODING)
    return FASTQIterator(FASTQParser(input), qual_encoding, false, nothing)
end


function advance!(it::FASTQIterator)
    it.isdone = !FASTQParserImpl.advance!(it.parser)
    if !it.isdone
        encoding = infer_quality_encoding(it.parser.qualbuf.data, 1,
                                          it.parser.qualbuf.pos - 1,
                                          it.default_qual_encoding)
        it.default_qual_encoding = encoding

        qscores = decode_quality_string(encoding, it.parser.qualbuf.data,
                                        1, it.parser.qualbuf.pos - 1)

        it.nextitem =
            FASTQSeqRecord(it.parser.namebuf,
                           DNASequence(it.parser.seqbuf.data, 1, it.parser.seqbuf.pos - 1),
                           FASTQMetadata(it.parser.descbuf, qscores))
        empty!(it.parser.seqbuf)
        empty!(it.parser.qualbuf)
    end
end


function start(it::FASTQIterator)
    advance!(it)
    return nothing
end


function next(it::FASTQIterator, state)
    item = it.nextitem
    advance!(it)
    return item, nothing
end


function done(it::FASTQIterator, state)
    return it.isdone
end


