using Statistics, FFTW

include("./WaveletTransform.jl")
include("./smooth2D.jl")
include("./windowFn.jl")
include("./types.jl")

function WaveletCoherency(
    x::Wt, y::Wt;
    window_type_t::String = "bar",
    window_type_s::String = "bar",
    window_size_t::Union{Float64, Int64} = 5,
    window_size_s::Union{Float64, Int64} = 1/4)

    (x.lowerPeriod != y.lowerPeriod) && throw(ArgumentError("the lower period of x and y should be identical"))
    (x.upperPeriod != y.upperPeriod) && throw(ArgumentError("the upper period of x and y should be identical"))
    (x.dt != y.dt) && throw(ArgumentError("the dt of x and y should be identical"))
    (x.dj != y.dj) && throw(ArgumentError("the dj of x and y should be identical"))

    Period = x.Period
    Scale  = x.Scale
    nc = x.nc
    nr = x.nr
    dt = x.dt
    dj = x.dj

    #############################################################################
    ## compute cross-wavelet transform and cross-wavelet power
    #############################################################################

    Wave_xy = (x.Wave .* conj(y.Wave)) ./ reshape(repeat(Scale, nc), (nr, :))
    Power_xy = abs.(Wave_xy)

    # odd window sizes, given in terms of dj, dt resolution
    window_size_s = 2*floor(Int64, window_size_s/(2*dj))+1
    window_size_t = 2*floor(Int64, window_size_t/(2*dt))+1

    # window 2D matrix
    window2D = windowFn(window_type_s, window_size_s) * windowFn(window_type_t, window_size_t)'

    # smooth individual wavelet powers and smooth the cross-wavelet transform
    sPower_x = real(smooth2D(x.Power, window2D))
    sPower_y = real(smooth2D(y.Power, window2D))
    sWave_xy = smooth2D(Wave_xy, window2D)

    sPower_x_y = sPower_x .* sPower_y
    coh_replace = (sPower_x_y) .== 0
    coh_replace_len = sum(coh_replace)

    Coherency = sWave_xy ./ sqrt.(sPower_x_y)
    Coherence = abs.(sWave_xy).^2 ./ (sPower_x_y)

    Coherency[coh_replace] = repeat([0.0 + 0im], coh_replace_len)
    Coherence[coh_replace] = zeros(Float64, coh_replace_len)

    return Wc(  Wave      = Wave_xy,     sWave     = sWave_xy,
                Coherence = Coherence,   Coherency = Coherency,
                Period    = Period,      Scale     = Scale,
                Power     = Power_xy,
                nc        = nc,          nr        = nr,
                Wt_x      = blend_wt_with_sPower(x, sPower_x),
                Wt_y      = blend_wt_with_sPower(y, sPower_y))
end

function WaveletCoherency(
    x::Union{Vector{Float64}, Vector{Int64}},
    y::Union{Vector{Float64}, Vector{Int64}};
    dt::Union{Float64, Int64} = 1, dj::Union{Float64, Int64} = 1/20,
    lowerPeriod::Union{Float64, Int64} = 2*dt,
    upperPeriod::Union{Float64, Int64} = floor(length(x)*dt/3),
    window_type_t::String = "bar",
    window_type_s::String = "bar",
    window_size_t::Union{Float64, Int64} = 5,
    window_size_s::Union{Float64, Int64} = 1/4)

    Wt_x = WaveletTransform(
        x,
        dt = dt, dj = dj,
        lowerPeriod = lowerPeriod,
        upperPeriod = upperPeriod)

    Wt_y = WaveletTransform(
        y,
        dt = dt, dj = dj,
        lowerPeriod = lowerPeriod,
        upperPeriod = upperPeriod)

    return WaveletCoherency(
        Wt_x, Wt_y,
        window_type_t = window_type_t,
        window_type_s = window_type_s,
        window_size_t = window_size_t,
        window_size_s = window_size_s)
end
