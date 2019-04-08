module StopIon
using Unitful

export Particle, hello, domath

struct Particle
    z::Float64
    m::Unitful.Mass
end

hello(who::String) = "Hello, $who"
domath(x::Number) = (x + 5)

end
