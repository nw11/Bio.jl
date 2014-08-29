
# Update the BioFmtSpecimen repository if it doesn't exist
if !isdir("BioFmtSpecimens")
    run(`git clone https://github.com/BioJulia/BioFmtSpecimens.git`)
end

include("align/test_align.jl")
include("phylo/test_phylo.jl")
include("ranges/test_ranges.jl")
include("seq/test_seq.jl")
include("services/test_services.jl")
include("tools/test_tools.jl")
