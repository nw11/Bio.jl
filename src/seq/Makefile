
 # Run ragel to regenerate parsers after updating the .rl source file

parsers=fasta.jl

all: fasta.jl

clean:
	rm -f ${parsers}

%.jl: %.rl
	ragel -U -G2 $<

