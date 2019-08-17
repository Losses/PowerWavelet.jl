function SurrogateData(
    x::Union{Vector{Float64}, Vector{Int64}};
    method::String = "white.noise", params::Dict{String, Real} = Dict(
      "p"      => 1, "q"    => 1,     "include_mean" => true,
      "sd_fac" => 1, "trim" => false, "trim_prop"    => 0.01))

    method == "white.noise" && return randn(length(x))
    method == "shuffle" && return shuffle(length(x))

    throw(ArgumentError(method, "surrgate method not implemented or not exists"))
end
