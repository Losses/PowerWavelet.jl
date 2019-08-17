include("./WaveletCoherency.jl")
include("./wt.jl")
include("./types.jl")

function wc(
    x::Wc;
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

    #check here!

    dt = x.Wt_x.dt
    dj = x.Wt_x.dj
    nc = x.nc
    nr = x.nr

    lowerPeriod = x.Wt_x.lowerPeriod
    upperPeriod = x.Wt_x.upperPeriod

    Wave_xy = x.Wave
    sWave_xy = x.sWave

    Angle = angle.(Wave_xy)
    sAngle = angle.(sWave_xy)

    Power_xy = x.Power
    Power_xy_avg = Statistics.mean(Power_xy, dims = 2)[:, 1]

    Coherency = x.Coherency
    Coherence = x.Coherence
    Coherence_avg = Statistics.mean(Coherence, dims = 2)[:, 1]

    # individual wavelet results

    Wave_x = x.Wt_x.Wave
    Wave_y = x.Wt_y.Wave
    Phase_x = x.Wt_x.Phase
    Phase_y = x.Wt_y.Phase
    Ampl_x = x.Wt_x.Ampl
    Ampl_y = x.Wt_y.Ampl

    Power_x = x.Wt_x.Power
    Power_y = x.Wt_y.Power
    Power_x_avg = x.Wt_x.Power_avg
    Power_y_avg = x.Wt_y.Power_avg

    sPower_x = x.Wt_x.sPower
    sPower_y = x.Wt_y.sPower

    # parameters

    Period = x.Period
    Scale  = x.Scale
    nr  = x.nr
    nc  = x.nc

    ###############################################################################
    ## Compute p values for significance check
    ###############################################################################

    Coherence_pval = nothing
    Power_xy_pval = nothing
    Power_x_pval = nothing
    Power_y_pval = nothing

    Coherence_avg_pval = nothing
    Power_xy_avg_pval = nothing
    Power_x_avg_pval = nothing
    Power_y_avg_pval = nothing

    series_sim_x = nothing
    series_sim_y = nothing

    if (make_pval)
        Coherence_pval = zeros(Float64, nr, nc)
        Power_xy_pval = zeros(Float64, nr, nc)
        Power_x_pval = zeros(Float64, nr, nc)
        Power_y_pval = zeros(Float64, nr, nc)

        Coherence_avg_pval = repeat([0], nr)
        Power_xy_avg_pval = repeat([0], nr)
        Power_x_avg_pval = repeat([0], nr)
        Power_y_avg_pval = repeat([0], nr)

        if (save_sim)
            series_sim_x = reshape(repeat([0.0], nc * n_sim), (n_sim, :))
            series_sim_y = reshape(repeat([0.0], nc * n_sim), (n_sim, :))
        end

        for ind_sim = 1:n_sim
            x_sim = SurrogateData(x, method = method, params = params)
            y_sim = SurrogateData(x, method = method, params = params)

            if (save_sim)
                series_sim_x[ind_sim, :] = x_sim
                series_sim_y[ind_sim, :] = y_sim
            end

            Wc_sim = WaveletCoherency(
                x_sim, y_sim,
                lowerPeriod = lowerPeriod,
                upperPeriod = upperPeriod,
                window_size_s = window_size_s,
                window_size_t = window_size_t,
                window_type_s = window_type_s,
                window_type_t = window_type_t)

            Coherence_sim = Wc_sim.Coherence
            Power_sim_xy = Wc_sim.Power
            Power_sim_x = Wc_sim.Wt_x.Power
            Power_sim_y = Wc_sim.Wt_y.Power

            Coherence_avg_sim = Statistics.mean(Coherence_sim, dims = 2)[:, 1]
            Power_xy_avg_sim = Statistics.mean(Power_avg_xy_sim, dims = 2)[:, 1]
            Power_x_avg_sim = Statistics.mean(Power_avg_x_sim, dims = 2)[:, 1]
            Power_y_avg_sim = Statistics.mean(Power_avg_y_sim, dims = 2)[:, 1]

            Coherence_selector = abs(Coherence_sim) .>= abs(Coherence)
            Power_xy_selector = abs(Power_sim_xy) .>= abs(Power_xy)
            Power_x_selector = abs(Power_sim_x) .>= abs(Power_x)
            Power_y_selector = abs(Power_sim_y) .>= abs(Power_y)

            Coherence_avg_selector = abs(Coherence_sim_avg) .>= abs(Coherence_avg)
            Power_xy_avg_selector = abs(Power_sim_xy_avg) .>= abs(Power_xy_avg)
            Power_x_avg_selector = abs(Power_sim_x_avg) .>= abs(Power_x_avg)
            Power_y_avg_selector = abs(Power_sim_y_avg) .>= abs(Power_y_avg)

            Coherence_pval[Coherence_selector] = Coherence_pval[Coherence_selector] .+ 1
            Power_xy_pval[Power_xy_selector] = Power_xy_pval[Power_xy_selector] .+ 1
            Power_x_pval[Power_x_selector] = Power_x_pval[Power_x_selector] .+ 1
            Power_y_pval[Power_y_selector] = Power_y_pval[Power_y_selector] .+ 1

            Coherence_avg_pval[Coherence_avg_selector] = Coherence_avg_pval[Coherence_avg_selector] .+ 1
            Power_xy_avg_pval[Power_xy_avg_selector] = Coherence_avg_pval[Coherence_avg_selector] .+ 1
            Power_x_avg_pval[Power_x_avg_selector] = Coherence_avg_pval[Coherence_avg_selector] .+ 1
            Power_y_avg_pval[Power_y_avg_selector] = Coherence_avg_pval[Coherence_avg_selector] .+ 1
        end

        #p-values

        Coherence_pval = Coherence_pval ./ n_sim
        Power_xy_pval = Power_xy_pval ./ n_sim
        Power_x_pval = Power_x_pval ./ n_sim
        Power_y_pval = Power_y_pval ./ n_sim

        Coherence_avg_pval = Coherence_avg_pval ./ n_sim
        Power_xy_avg_pval = Power_xy_avg_pval ./ n_sim
        Power_x_avg_pval = Power_x_avg_pval ./ n_sim
        Power_y_avg_pval = Power_y_avg_pval ./ n_sim
    end

    Coi = COI(start = start, dt = dt, nc = nc, nr = nr, Period = Period)

    return blend_wc_with_p(
        x,
        blend_wt_with_p(
            x.Wt_x,
            Power_x_avg, Power_x_pval,
            Power_x_avg_pval,
            Coi, series_sim_x),
        blend_wt_with_p(
            x.Wt_y,
            Power_y_avg, Power_y_pval,
            Power_y_avg_pval,
            Coi, series_sim_y),
        Power_xy_avg,
        Power_xy_pval,
        Power_xy_avg_pval,
        Coi,
        series_sim_x,
        series_sim_y)
end


function wc(
    x::Union{Vector{Float64}, Vector{Int64}},
    y::Union{Vector{Float64}, Vector{Int64}};
    start::Int = 1,
    dt::Union{Float64, Int64} = 1, dj::Union{Float64, Int64} = 1/20,
    lowerPeriod::Union{Float64, Int64} = 2*dt,
    upperPeriod::Union{Float64, Int64} = floor(length(x)*dt/3),
    window_type_t::String = "bar",
    window_type_s::String = "bar",
    window_size_t::Union{Float64, Int64} = 5,
    window_size_s::Union{Float64, Int64} = 1/4,
    make_pval::Bool = true, method::String = "white.noise",
    params::Dict{String, Real} = Dict(
      "p"      => 1, "q"    => 1,     "include_mean" => true,
      "sd_fac" => 1, "trim" => false, "trim_prop"    => 0.01),
    n_sim::Int64 = 100, save_sim::Bool = false)
    Wt_x = wt(
        x,
        start = start,
        dt = dt,
        dj = dj,
        lowerPeriod = lowerPeriod,
        upperPeriod = upperPeriod,
        make_pval = false)

    Wt_y = wt(
        y,
        start = start,
        dt = dt,
        dj = dj,
        lowerPeriod = lowerPeriod,
        upperPeriod = upperPeriod,
        make_pval = false)


    Wc_xy = WaveletCoherency(
        Wt_x, Wt_y,
        window_type_t = window_type_t,
        window_type_s = window_type_s,
        window_size_t = window_size_t,
        window_size_s = window_size_s)

    return wc(
        Wc_xy,
        start = start,
        make_pval = make_pval,
        method = method,
        params = params,
        n_sim = n_sim,
        save_sim = save_sim)
end
