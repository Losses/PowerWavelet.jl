function ridge(
    wavelet_spectrum::Matrix{Float64},
    band::Int64 = 5, scale_factor::Float64 = 0.1)

    min_level = scale_factor * maximum(wavelet_spectrum)
    nrows = size(wavelet_spectrum)[1]
    ncols = size(wavelet_spectrum)[2]

    Ridge = reshape(zeros(nrows*ncols), nrows, :)

    for col = 1:ncols
        column_vec = wavelet_spectrum[:, col]

        ind = [1:1:nrows;]

        band_max_vec = column_vec

        for i =1:band
            lower_ind  = ind .- i
            lower_ind[lower_ind .< 1] = ones(sum(lower_ind .< 1))
            upper_ind  = ind .+ i
            upper_ind[upper_ind .> nrows] = repeat([nrows], sum(upper_ind .> nrows))

            band_max_vec = max.(band_max_vec, column_vec[lower_ind], column_vec[upper_ind])
        end

        my_ridge_column = zeros(nrows)
        my_ridge_column_selector =  band_max_vec .== column_vec
        my_ridge_column[my_ridge_column_selector] = repeat([1], sum(my_ridge_column_selector))

        Ridge[:, col] = my_ridge_column
    end

    return Ridge .* (wavelet_spectrum .> min_level)
end
