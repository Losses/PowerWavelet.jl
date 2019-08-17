using Statistics, FFTW

include("./types.jl")

function WaveletTransform(
    x::Union{Vector{Float64}, Vector{Int64}};
    dt::Union{Float64, Int64} = 1, dj::Union{Float64, Int64} = 1/20,
    lowerPeriod::Union{Float64, Int64} = 2*dt,
    upperPeriod::Union{Float64, Int64} = floor(length(x)*dt/3))

    ###############################################################################
    ## Provide parameters (which could be useful for other transforms as well)
    ###############################################################################

    # Original length and length of zero padding:
    series_length = length(x)
    pot2 = trunc(log2(series_length) + 0.5)
    pad_length = trunc(Int, 2 ^ (pot2 + 1) - series_length)

    # Define central angular frequency omega0 and fourier factor:
    omega0 = 6
    fourier_factor = (2 * pi) / omega0

    # Compute scales and periods:
    min_scale = lowerPeriod / fourier_factor             # Convert lowerPeriod to minimum scale
    max_scale = upperPeriod / fourier_factor             # Convert upperPeriod to maximum scale
    J = floor(Int, log2(max_scale / min_scale) / dj)     # Index of maximum scale -1

    scales = min_scale .* 2 .^ ([0:(J/abs(J)):J;] .* dj)          # sequence of scales
    scale_length = length(scales)                                 # J + 1
    periods = fourier_factor * scales                             # sequence of periods

    # Computation of the angular frequencies
    N = series_length + pad_length
    _omega_k = [1:1:floor(N/2);] .* ((2*pi)/(N*dt))
    _omega_k_3 = -1 * _omega_k[reverse([1:1:floor(Int, (N-1) / 2);])]

    omega_k = [0; _omega_k; _omega_k_3]

    ###############################################################################
    ## Define the Morlet wavelet transform function
    ###############################################################################

    function morletWaveletTransform(x::Union{Vector{Float64}, Vector{Int64}})

        # Standardize x and pad with zeros
        x = (x .- Statistics.mean(x)) ./ Statistics.std(x)
        xpad = [x; repeat([0.0], pad_length)]

        # Compute Fast Fourier Transform of xpad
        fft_xpad = fft(xpad)
        wave = repeat([0.0], scale_length, N)
        wave = wave .+ 1im .* wave

        # Computation for each scale...
        # ... simultaneously for all time instances
        for ind_scale = 1:scale_length
            my_scale = scales[ind_scale]

            norm_factor = pi^(1/4) * sqrt(2 * my_scale / dt)
            expnt = -((my_scale .* omega_k .- omega0).^2 ./ 2) .* (omega_k .> 0)
            daughter = (norm_factor .* exp.(expnt)) .* (omega_k .> 0)

            wave[ind_scale, :] = bfft(fft_xpad .* daughter) ./ N
        end

        wave = wave[:, 1:series_length]
        return wave
    end

    ###############################################################################
    ## Compute the wavelet transform, power, phases, amplitudes
    ###############################################################################

    Wave = morletWaveletTransform(x)

    # Compute wavelet power
    Power = abs.(Wave) .^2 ./ reshape(repeat(scales, series_length), (length(scales), series_length))

    # Phase
    Phase = angle.(Wave)

    # Amplitude
    Ampl  = abs.(Wave) ./ reshape(repeat(scales, series_length), (length(scales), series_length))

    return Wt(
        Wave        = Wave,
        Phase       = Phase,         Ampl  = Ampl,
        Period      = periods,       Scale = scales,
        Power       = Power,
        nc          = series_length, nr    = scale_length,
        dt          = dt,            dj    = dj,
        lowerPeriod = lowerPeriod,
        upperPeriod = upperPeriod)
end
