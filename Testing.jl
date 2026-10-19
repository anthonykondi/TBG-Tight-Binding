using AutoHashEquals
using StaticArrays

struct my_struct
    x::Float64
    y::SVector{2, Int}
end

a = my_struct(32, (2, 3))
c = my_struct(32, (2, 3))

b = Dict(a => "yes")

println(b[c])