using FFTW;

function smooth2D(mat::Union{Matrix{Float64}, Matrix{ComplexF64}}, window2D::Matrix{Float64})
    mat_rows = size(mat)[1]
    mat_cols = size(mat)[2]
    win_rows = size(window2D)[1]
    win_cols = size(window2D)[2]

    mat_pad = zeros(typeof(mat[1]), mat_rows + win_rows - 1, mat_cols + win_cols - 1)
    win_pad = zeros(Float64, mat_rows + win_rows - 1, mat_cols + win_cols - 1)

    mat_pad[1:mat_rows, 1:mat_cols] = mat
    win_pad[1:win_rows, 1:win_cols] = window2D

    smooth_mat = bfft(fft(mat_pad) .* fft(win_pad)) ./ length(mat_pad)

    return smooth_mat[floor(Int64, win_rows/2) .+ (1:mat_rows),
                      floor(Int64, win_cols/2) .+ (1:mat_cols)]
end
