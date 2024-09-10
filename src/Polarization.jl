mutable struct Photon
    E::Float64
    Q::Float64
    U::Float64
    V::Float64
    n::Array{Float64,1}{3}
    k::Array{Float64,1}{3}
end