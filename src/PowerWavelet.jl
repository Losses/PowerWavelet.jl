module PowerWavelet

include("WaveletTransform.jl")
include("WaveletCoherency.jl")
include("wt.jl")
include("wc.jl")
include("analyze_wavelet.jl")
include("analyze_coherence.jl")
include("rnd_signal.jl")

export WaveletTransform, WaveletCoherency,
    Wt, Wc,
    analyze_wavelet, analyze_coherence,
    rnd_signal

end # module
