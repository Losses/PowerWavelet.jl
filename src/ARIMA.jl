using ARIMA

DEFAULT_ARIMA_PARAMS = Dict("p" => 1, "q" => 1, "include_mean" => true, "sd_fac" => 1, "trim" => false, "trim_prop" => 0.01)
function ARIMA(x::Union{Vector{Int64}, Vector{Float64}}, params::Dict{String, Any} = DEFAULT_ARIMA_PARAMS)
    isa(params["p"], Number) || throw(TypeError(params["p"], "`p` should be a number"))
    isa(params["q"], Number) || throw(TypeError(params["q"], "`q` should be a number"))
    isa(params["include_mean"], Bool) || throw(TypeError(params["include_mean"], "`include_mean` should be a boolean"))
    isa(params["sd_fac"], Number) || throw(TypeError(params["sd_fac"], "`sd_fac` should be a number"))
    isa(params["trim"], Bool) || throw(TypeError(params["trim"], "`trim` should be a boolean"))
    isa(params["trim_prop"], Number) || throw(TypeError(params["trim_prop"], "`trim_prop` should be a number"))

    n = length(x)

    my_arima = ()
end
