
module Ragel

using Switch

# Nuber of bytes to read at a time
const RAGEL_PARSER_INITIAL_BUF_SIZE = 100000


# A type keeping track of a ragel-based parser's state.
type State
    # Input source. If nothing, it's assumed that buffer contains
    # all the input.
    input::Union(IO, Nothing)

    # Buffered data.
    buffer::Vector{Uint8}

    # Indexes into buffer that will be updated and stay valid upon refill.
    marks::Vector{Uint}

    # True when parsing has completed
    finished::Bool

    # Internal ragel state:
    p::Int  # index into the input stream (0-based)
    pe::Int # index after the last symbol in the input stream (0-based)
    cs::Int # current DFA stae

    function State(cs, data::Vector{Uint8})
        return new(nothing, data, Uint[], false, 0, length(data), cs)
    end

    function State(cs, input::IO)
        return new(input, Array(Uint8, RAGEL_PARSER_INITIAL_BUF_SIZE),
                   Uint[], false, 0, 0, cs)
    end

    function State(cs, filename::String, memory_map=true)
        if memory_map
            data = mmap_array(Uint8, (filesize(filename),), open(filename))
            return new(nothing, data, Uint[], false, 0, length(data), cs)
        else
            return new(input, Array(Uint8, RAGEL_PARSER_INITIAL_BUF_SIZE),
                       Uint[], false, 0, 0, cs)
        end
    end
end


# Get a State object from a parser. Parser implementations may want
# to define a more specific method.
function ragelstate(x)
    return x.state
end


# Macros for push and popping marks from within a ragel parser
macro pushmark!()
    quote
        push!($(esc(:state)).marks, 1 + $(esc(:p)))
    end
end

macro popmark!()
    quote
        pop!($(esc(:state)).marks)
    end
end

macro position()
    quote
        1 + $(esc(:p))
    end
end

macro spanfrom(firstpos)
    quote
        $(esc(:state)).buffer[$(esc(firstpos)):$(esc(:p))]
    end
end

# return the current character
macro char()
    quote
        $(esc(:state)).buffer[$(esc(:p))+1]
    end
end


# Refill a buffer, keeping some portion of it.
#
# The range buffer[parser.marks[1]:parser.pe] will be shifted to the beginning of
# the buffer, then the rest of the buffer will be filled by reading bytes from
# the input.
#
# This is useful if a parser is the middle of matching something when the buffer
# runs out
#
# Args:
#   parser: Current ragel parser state.
#
# Modifies:
#   parser.buffer is refilled, and indexes are updated to the correct positions
#   in the refilled buffer
#
# Returns:
#   Number of bytes read
#
function fillbuffer!(parser::State)
    if parser.input === nothing
        return 0
    end

    buflen = length(parser.buffer)
    keeplen = 0
    first_mark = 0
    if !isempty(parser.marks) > 0
        first_mark = parser.marsk[1]
        keeplen = parser.pe - first_mark + 1
        if keeplen == buflen
            buflen = 2 * buflen
            resize!(parser.buffer, buflen)
        end
        copy!(parser.buffer, 1, parser.buffer, first_mark, keeplen)
    end

    i = keeplen
    while i < buflen && !eof(parser.input)
        i += 1
        parser.buffer[i] = read(parser.input, Uint8)
    end

    parser.p = 0
    parser.pe = i
    for i in 1:length(parser.marks)
        parser.marks[i] -= first_mark - 1
    end

    return i - keeplen
end


# Define a read! function wrapping ragel-generated parser.
#
# This macro handles some the dirty work of maintaining state, refilling
# buffers, etc.
#
# Args:
#
macro generate_read_fuction(machine_name, input_type, output_type, ragel_body, accept_body)
    start_state = esc(symbol(string(machine_name, "_start")))
    accept_state = esc(symbol(string(machine_name, "_first_final")))
    error_state = esc(symbol(string(machine_name, "_error")))

    # ragel needs these specific variable names so we have to escape them
    # throughout
    p = esc(:p)
    pe = esc(:pe)
    cs = esc(:cs)
    data = esc(:data)
    state = esc(:state)

    quote
        function Base.read!(input::$(esc(input_type)), output::$(esc(output_type)))
            # TODO: is there a more idiomatic way to do this?
            local $(esc(:input)) = input
            local $(esc(:output)) = output

            $(state) = ragelstate(input)
            if $(state).finished
                return false
            end

            $(p) = $(state).p
            $(pe) = $(state).pe
            $(cs) = $(state).cs
            $(data) = $(state).buffer

            # run the parser until all input is consumed or a match is found
            while true
                if $(p) == $(pe)
                    if fillbuffer!($(state)) == 0
                        break
                    end

                    $(p) = $(state).p
                    $(pe) = $(state).pe
                end

                $(esc(ragel_body))

                @show $(cs)

                if $(cs) == $(error_state)
                    # TODO: better error messages would be nice. E.g. keeping
                    # track of the line number at the very least.
                    error($("Error parsing $(machine_name)"))
                #elseif $(p) != $(pe) && $(cs) >= $(accept_state)
                else
                    break
                end
            end
            $(state).p = $(p)
            $(state).pe = $(pe)
            $(state).cs = $(cs)

            # TODO: How do we actually check if anything was read???

            if $(p) == $(pe)
                $(state).finished = true
            end
            $(esc(accept_body))
            return true
        end
    end
end


end # module Ragel

