using Loess

include("./wt.jl")
include("./ridge.jl")

function analyze_wavelet(
    x:: Union{Vector{Float64}, Vector{Int64}};
    loess_span::Union{Float64, Int64} = 0.75,
    dt::Union{Float64, Int64} = 1, dj::Union{Float64, Int64} = 1/20,
    lowerPeriod::Union{Float64, Int64} = 2*dt,
    upperPeriod::Union{Float64, Int64} = floor(length(x)*dt/3),
    make_pval::Bool = true, method::String = "white.noise",
    n_sim::Int64 = 100,
    params::Dict{String, Real} = Dict(
      "p"      => 1, "q"    => 1,     "include_mean" => true,
      "sd_fac" => 1, "trim" => false, "trim_prop"    => 0.01))

    typeof(x) == Vector{Int64} && (x = convert(Vector{Float64}, x))

    ###################################################################################################
    ## Smooth the data (if requested)
    ###################################################################################################

    if (loess_span != 0)
        day_index = convert(Vector{Float64}, [1:1:length(x); ])
        my_loess_x = loess(day_index, x, span = loess_span)
        x_trend = predict(my_loess_x, day_index)

        x = x .- x_trend
    end

    ###################################################################################################
    ## Start the analysis of wavelets
    ###################################################################################################

    my_wt = wt(
        x,
        start = 1,
        dt = dt, dj = dj,
        lowerPeriod = lowerPeriod, upperPeriod = upperPeriod,
        make_pval = make_pval, method = method,
        params = params,
        n_sim = n_sim, save_sim = false)

    ##################################################################################################
    ## Compute the power ridge
    ##################################################################################################

    Ridge = ridge(my_wt.Power)

    return Wt(
        Wave        = my_wt.Wave,
        Phase       = my_wt.Phase,       Ampl           = my_wt.Ampl,
        Period      = my_wt.Period,      Scale          = my_wt.Scale,
        Power       = my_wt.Power,
        nc          = my_wt.nc,          nr             = my_wt.nr,
        dt          = my_wt.dt,          dj             = my_wt.dj,
        lowerPeriod = my_wt.lowerPeriod,
        upperPeriod = my_wt.upperPeriod,
        Power_avg   = my_wt.Power_avg,
        Power_pval  = my_wt.Power_pval,  Power_avg_pval = my_wt.Power_avg_pval,
        COI         = my_wt.COI       ,  series_sim     = my_wt.series_sim,
        Ridge       = Ridge)
end
