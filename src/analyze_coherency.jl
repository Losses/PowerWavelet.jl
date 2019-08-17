include("./wc.jl")
inculde("./ridge.jl")

function analyze_coherency(
    x::Union{Vector{Float64}, Vector{Int64}},
    y::Union{Vector{Float64}, Vector{Int64}};
    loess_span::Union{Float64, Int64} = 0.75,
    start::Int = 1,
    window_type_t::String = "bar",
    window_type_s::String = "bar",
    window_size_t::Union{Float64, Int64} = 5,
    window_size_s::Union{Float64, Int64} = 1/4,
    make_pval::Bool = true, method::String = "white.noise",
    params::Dict{String, Real} = Dict(
      "p"      => 1, "q"    => 1,     "include_mean" => true,
      "sd_fac" => 1, "trim" => false, "trim_prop"    => 0.01),
    n_sim::Int64 = 100, save_sim::Bool = false)

    typeof(x) == Vector{Int64} && (x = convert(Vector{Float64}, x))
    typeof(y) == Vector{Int64} && (y = convert(Vector{Float64}, y))

    ###################################################################################################
    ## Smooth the data (if requested)
    ###################################################################################################

    if (loess_span != 0)
        day_index = convert(Vector{Float64}, [1:1:length(x); ])
        my_loess_x = loess(day_index, x, span = loess_span)
        x_trend = predict(my_loess_x, day_index)
        my_loess_y = loess(day_index, y, span = loess_span)
        y_trend = predict(my_loess_y, day_index)

        x = x .- x_trend
        y = y .- y_trend
    end

    ###################################################################################################
    ## Start the analysis of wavelet coherency
    ###################################################################################################

    my_wc = wc(
        x, y,
        start = 1,
        dt = dt, dj = dj,
        lowerPeriod = lowerPeriod, upperPeriod = upperPeriod,
        window_type_t = window_type_t,
        window_type_s = window_type_s,
        window_size_t = window_size_t,
        window_size_s = window_size_s,
        make_pval = make_pval, method = method,
        params = params,
        n_sim = n_sim, save_sim = false)

    ##################################################################################################
    ## Compute the ridges
    ##################################################################################################

    Ridge_xy = ridge(my_wc.Power)
    Ridge_co = ridge(my_wc.Coherence)

    Ridge_x = ridge(my_wc.Wt_x.Power)
    Ridge_y = ridge(my_wc.Wt_y.Power)

    
end
