function rnd_signal(
    len::Int64 = 100,
    power_noise::Float64 = 0.3,
    freq_noise::Float64 = 0.3,
    freq_times::Int64 = 3)

    x = [1:1:len;]
    y = sin.(x)

    for i = 1:freq_times
        freq = rand(Float64) * freq_noise
        y = y + sin.(x * freq)
    end

    y .+ (rand(Float64, len) .* power_noise)
end
