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


const fasta_start  = convert(Int , 7)
const fasta_first_final  = convert(Int , 7)
const fasta_error  = convert(Int , 0)
const fasta_en_main  = convert(Int , 7)
type FASTAParser
    state::Ragel.State
    seqbuf::Ragel.Buffer

    function FASTAParser(input::Union(IO, String, Vector{Uint8}))
        begin
cs = fasta_start;
	end
return new(Ragel.State(cs, input), Ragel.Buffer{Uint8}())
    end
end


function Ragel.ragelstate(parser::FASTAParser)
    return parser.state
end


function accept_state!{S}(input::FASTAParser, output::FASTASeqRecord{S})
    output.seq = S(input.seqbuf.data, 1, input.seqbuf.pos - 1)
    empty!(input.seqbuf)
end


Ragel.@generate_read_fuction("fasta", FASTAParser, FASTASeqRecord,
    begin
        @inbounds begin
            begin
if p == pe 
	@goto _test_eof

end
@switch cs  begin
    @case 7
@goto st_case_7
@case 0
@goto st_case_0
@case 1
@goto st_case_1
@case 2
@goto st_case_2
@case 3
@goto st_case_3
@case 4
@goto st_case_4
@case 8
@goto st_case_8
@case 5
@goto st_case_5
@case 6
@goto st_case_6

end
@goto st_out
@label st7
p+= 1;
	if p == pe 
	@goto _test_eof7

end
@label st_case_7
@switch ( data[1 + p ])  begin
    @case 32
begin
@goto st7
end
@case 62
begin
@goto st1
end

end
if 9 <= ( data[1 + p ]) && ( data[1 + p ]) <= 13 
	begin
@goto st7
end

end
begin
@goto st0
end
@label st_case_0
@label st0
cs = 0;
	@goto _out
@label ctr14
begin
	yield = true;
        begin
	p+= 1; cs = 1; @goto _out

end

    
end
@goto st1
@label st1
p+= 1;
	if p == pe 
	@goto _test_eof1

end
@label st_case_1
if ( data[1 + p ]) == 32 
	begin
@goto st0
end

end
if ( data[1 + p ]) < 14 
	begin
if 9 <= ( data[1 + p ]) 
	begin
@goto st0
end

end
end

elseif ( ( data[1 + p ]) > 31  )
	begin
if 33 <= ( data[1 + p ]) 
	begin
@goto ctr0
end

end
end

else
	begin
@goto ctr0
end

end
begin
@goto ctr0
end
@label ctr0
begin
	Ragel.@pushmark!
    
end
@goto st2
@label st2
p+= 1;
	if p == pe 
	@goto _test_eof2

end
@label st_case_2
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto ctr3
end
@case 10
begin
@goto ctr4
end
@case 11
begin
@goto ctr3
end
@case 12
begin
@goto st0
end
@case 13
begin
@goto ctr4
end
@case 32
begin
@goto ctr3
end

end
if ( data[1 + p ]) > 31 
	begin
if 33 <= ( data[1 + p ]) 
	begin
@goto st2
end

end
end

elseif ( ( data[1 + p ]) >= 14  )
	begin
@goto st2
end

end
begin
@goto st2
end
@label ctr3
begin
	firstpos = Ragel.@popmark!
        output.name = bytestring(Ragel.@spanfrom firstpos)
    
end
@goto st3
@label st3
p+= 1;
	if p == pe 
	@goto _test_eof3

end
@label st_case_3
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto ctr6
end
@case 11
begin
@goto ctr6
end
@case 12
begin
@goto ctr5
end
@case 32
begin
@goto ctr6
end

end
if ( data[1 + p ]) < 14 
	begin
if 10 <= ( data[1 + p ]) 
	begin
@goto st0
end

end
end

elseif ( ( data[1 + p ]) > 31  )
	begin
if 33 <= ( data[1 + p ]) 
	begin
@goto ctr5
end

end
end

else
	begin
@goto ctr5
end

end
begin
@goto ctr5
end
@label ctr5
begin
	Ragel.@pushmark!
    
end
@goto st4
@label st4
p+= 1;
	if p == pe 
	@goto _test_eof4

end
@label st_case_4
@switch ( data[1 + p ])  begin
    @case 10
begin
@goto ctr8
end
@case 13
begin
@goto ctr8
end

end
if ( data[1 + p ]) > 12 
	begin
if 14 <= ( data[1 + p ]) 
	begin
@goto st4
end

end
end

elseif ( ( data[1 + p ]) >= 11  )
	begin
@goto st4
end

end
begin
@goto st4
end
@label ctr4
begin
	firstpos = Ragel.@popmark!
        output.name = bytestring(Ragel.@spanfrom firstpos)
    
end
@goto st8
@label ctr8
begin
	firstpos = Ragel.@popmark!
        output.metadata.description = bytestring(Ragel.@spanfrom firstpos)
    
end
@goto st8
@label ctr9
begin
	firstpos = Ragel.@popmark!
        append!(input.seqbuf, state.buffer, firstpos, p)
    
end
@goto st8
@label st8
p+= 1;
	if p == pe 
	@goto _test_eof8

end
@label st_case_8
@switch ( data[1 + p ])  begin
    @case 32
begin
@goto st8
end
@case 62
begin
@goto ctr14
end

end
if ( data[1 + p ]) < 65 
	begin
if 9 <= ( data[1 + p ]) && ( data[1 + p ]) <= 13 
	begin
@goto st8
end

end
end

elseif ( ( data[1 + p ]) > 90  )
	begin
if 97 <= ( data[1 + p ]) && ( data[1 + p ]) <= 122 
	begin
@goto ctr15
end

end
end

else
	begin
@goto ctr15
end

end
begin
@goto st0
end
@label ctr15
begin
	Ragel.@pushmark!
    
end
@goto st5
@label st5
p+= 1;
	if p == pe 
	@goto _test_eof5

end
@label st_case_5
if ( data[1 + p ]) == 32 
	begin
@goto ctr9
end

end
if ( data[1 + p ]) < 65 
	begin
if 9 <= ( data[1 + p ]) && ( data[1 + p ]) <= 13 
	begin
@goto ctr9
end

end
end

elseif ( ( data[1 + p ]) > 90  )
	begin
if 97 <= ( data[1 + p ]) && ( data[1 + p ]) <= 122 
	begin
@goto st5
end

end
end

else
	begin
@goto st5
end

end
begin
@goto st0
end
@label ctr6
begin
	Ragel.@pushmark!
    
end
@goto st6
@label st6
p+= 1;
	if p == pe 
	@goto _test_eof6

end
@label st_case_6
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto ctr6
end
@case 10
begin
@goto ctr8
end
@case 11
begin
@goto ctr6
end
@case 12
begin
@goto ctr5
end
@case 13
begin
@goto ctr8
end
@case 32
begin
@goto ctr6
end

end
if ( data[1 + p ]) > 31 
	begin
if 33 <= ( data[1 + p ]) 
	begin
@goto ctr5
end

end
end

elseif ( ( data[1 + p ]) >= 14  )
	begin
@goto ctr5
end

end
begin
@goto ctr5
end
@label st_out
@label _test_eof7
cs = 7; @goto _test_eof
@label _test_eof1
cs = 1; @goto _test_eof
@label _test_eof2
cs = 2; @goto _test_eof
@label _test_eof3
cs = 3; @goto _test_eof
@label _test_eof4
cs = 4; @goto _test_eof
@label _test_eof8
cs = 8; @goto _test_eof
@label _test_eof5
cs = 5; @goto _test_eof
@label _test_eof6
cs = 6; @goto _test_eof
@label _test_eof
begin
end
if p == eof 
	begin
@switch cs  begin
    @case 8
begin
	yield = true;
        begin
	p+= 1; cs = 0; @goto _out

end

    
end

	break;
	
end
end

end
@label _out
begin
end
end
end
    end,
    begin
        accept_state!(input, output)
    end)

end # module FASTAParserImpl


using Bio.Seq.FASTAParserImpl


