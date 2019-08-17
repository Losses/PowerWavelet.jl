function COI(
    ;start::Int64, dt::Union{Float64, Int64},
    nc::Int64, nr::Int64,
    Period::Union{Vector{Float64}, Vector{Int64}})

    axis1 = [range(start, step = dt, length = nc); ]
    axis2 = log2.(Period)

    # Define central angular frequency omega0 and fourier factor:
    omega0 = 6.0
    # fourier.factor = (4*pi)/(omega0 + sqrt(2+omega0^2))
    fourier_factor = (2*pi)/omega0

    coi = fourier_factor .* sqrt(2) .* dt .* [1e-5; 1:1:((nc+1)/2-1); reverse(1:1:(nc/2-1)); 1e-5]
    coi_x = [repeat([axis1[1]], 2) .- (dt.*0.5); axis1; repeat([axis1[nc]], 2) .+ (dt.*0.5)]
    logyint = axis2[2] - axis2[1]
    yl = [log2(Period[nr]) + 0.5*logyint; log2(Period[1]) - 0.5*logyint]
    yr = reverse(yl)

    coi_y = [yl; log2.(coi); yr]

    return Coi(coi_1 = coi_x, coi_2 = coi_y, axis_1 = axis1, axis_2 = axis2)
end

## Major part of this code: Tian, H. and Cazelles, B., \code{WaveletCo}
