include("./WaveletTransform.jl")
include("./SurrogateData.jl")
include("./COI.jl")
include("./types.jl")

function wt(
    x:: Union{Vector{Float64}, Vector{Int64}};
    start::Int = 1,
    dt::Union{Float64, Int64} = 1, dj::Union{Float64, Int64} = 1/20,
    lowerPeriod::Union{Float64, Int64} = 2*dt,
    upperPeriod::Union{Float64, Int64} = floor(length(x)*dt/3),
    make_pval::Bool = true, method::String = "white.noise",
    params::Dict{String, Real} = Dict(
      "p"      => 1, "q"    => 1,     "include_mean" => true,
      "sd_fac" => 1, "trim" => false, "trim_prop"    => 0.01),
    n_sim::Int64 = 100, save_sim::Bool = false)
    ###############################################################################
    ## Call function WaveletTransform
    ## Retrieve the wavelet transform, power, phases, amplitudes
    ###############################################################################
    # wavelet transform
    WT = WaveletTransform(x, dt = dt, dj = dj, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod)

    Wave = WT.Wave
    Phase = WT.Phase
    Ampl = WT.Ampl

    Power = WT.Power
    Power_avg = Statistics.mean(Power, dims = 2)[:, 1]

    Period = WT.Period
    Scale = WT.Scale
    nr = WT.nr
    nc = WT.nc


    ###############################################################################
    ## Compute p values for significance check
    ###############################################################################

    Power_pval = nothing
    Power_avg_pval = nothing
    series_sim = nothing

    if (make_pval)
        Power_pval = zeros(Float64, nr, nc)
        Power_avg_pval = repeat([0], nr)

        save_sim && (
            series_sim = reshape(repeat([0.0], nc * n_sim), (n_sim, :))
        )

        for ind_sim = 1:n_sim
            x_sim = SurrogateData(x, method = method, params = params)

            save_sim && (series_sim[ind_sim, :] = x_sim)

            WT_sim = WaveletTransform(x_sim, dt = dt, dj = dj, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod)

            Power_sim = WT_sim.Power
            Power_avg_sim = Statistics.mean(Power_sim, dims = 2)

            Power_pval[Power_sim .>= Power] .= Power_pval[Power_sim .>= Power] .+ 1

            _selector = (Power_avg_sim .>= Power_avg)[:, 1]
            Power_avg_pval[_selector] .= Power_avg_pval[_selector] .+ 1
        end

        #p-values
        Power_pval = Power_pval ./ n_sim
        Power_avg_pval = Power_avg_pval ./ n_sim
    end

    ###############################################################################
    ## Compute the cone of influence COI
    ###############################################################################

    Coi = COI(start = start, dt = dt, nc = nc, nr = nr, Period = Period)

    ###############################################################################
    ## Prepare the output
    ###############################################################################

    return blend_wt_with_p(
        WT,
        Power_avg, Power_pval,
        Power_avg_pval,
        Coi, series_sim)
end
