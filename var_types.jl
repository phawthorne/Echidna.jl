#
# Creates types MOGA_Real, MOGA_Binary, and MOGA_Integer which can be used
# to create decision variable lists. MOGA_Real and MOGA_Integer have min
# and max values as part of the struct.
#
# Extends Base.Random.rand so that when given an instance of one of these types,
# returns a valid random value (i.e. within the specified bounds)


import Base.Random.rand

abstract type MOGA_Type end

struct MOGA_Real <: MOGA_Type
    min_value::Float64
    max_value::Float64
end

struct MOGA_Binary <: MOGA_Type
end

struct MOGA_Integer <: MOGA_Type
    min_value::Int64
    max_value::Int64
end

function rand(t::MOGA_Real)
    return (t.max_value - t.min_value) * rand() + t.min_value
end

function rand(t::MOGA_Binary)
    return rand() > 0.5
end

function rand(t::MOGA_Integer)
    return Int64(floor((1 + t.max_value - t.min_value) * rand()) + t.min_value)
end
