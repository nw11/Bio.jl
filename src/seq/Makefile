
 # Run ragel to regenerate parsers after updating the .rl source file

parsers=fasta.jl fasta-manual.jl

all: ${parsers}

clean:
	rm -f ${parsers}

%.jl: %.rl
	ragel -U -G2 $<

