function windowFn(type::String = "bar", n::Int64 = 5)
"""
| type::str | Description                                       |
|:---------:|:--------------------------------------------------|
| none      | No smoothing                                      |
| bar       | Bartlett  window   (triangular window with L=n-1) |
| tri       | Triangular window  (L=n)                          |
| box       | Rectangular  window  (Boxcar or Dirichlet window) |
| han       | Hanning  window                                   |
| ham       | Hamming  window                                   |
| bla       | Blackmann  window                                 |
"""

############################################################################
## Windows to use for smoothing cross-wavelet transform and individual wavelet power
## Inspired by:
## Luis Aguiar-Conraria and Maria Joana Soares, "GWPackage"
############################################################################

  if (type == "none")
      window = 1
  elseif (type == "bar")
      (n<=2) && throw(ArgumentError(n, "Bartlett window requires minimum size 3"))

      window = 1 .- abs.([0:1:(n-1);] .- (n-1)/2) ./ ((n-1)/2)
  elseif (type == "tri")
      (n<=1) && throw(ArgumentError(n, "Triangular (non-Bartlett) window requires minimum size 2"))

      window = 1 .- abs.([0:1:(n-1);] .- (n-1)/2) ./ (n/2)
  elseif (type == "box")
      window = rep([1], n)
  elseif (type == "han")
      (n<=2) && throw(ArgumentError(n, "Hanning window requires minimum size 3"))

      window = 0.5 .- 0.5 .* cos.(2*pi .* [0:1:(n-1);] ./ (n-1))
  elseif (type == "ham")
      (n<=1) && throw(ArgumentError(n, "Hamming window requires minimum size 2"))

      window = 0.53836 .- (1-0.53836) .* cos.(2*pi .* [0:1:(n-1);] ./ (n-1))
  elseif (type == "bla")
      (n<=1) && throw(ArgumentError(n, "Blackman window requires minimum size 2"))

      window = 7938/18608 .- (9240/18608) .* cos.(2*pi .* [0:1:(n-1);]/(n-1)) .+ (1430/18608) .* cos.(4*pi .* [0:1:(n-1);] ./ (n-1))
  else
    throw(ArgumentError(type, "invalid window type"))
  end

  return window ./ sum(window)
end
