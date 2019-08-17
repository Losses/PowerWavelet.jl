using Parameters

@with_kw struct Coi
    coi_1::Vector{Float64}
    coi_2::Vector{Float64}
    axis_1::Union{Vector{Float64}, Vector{Int64}}
    axis_2::Union{Vector{Float64}, Vector{Int64}}
end

@with_kw struct Wt
    Wave::Matrix{Complex}
    Ampl::Matrix{Float64}
    Period::Union{Vector{Float64}, Vector{Int64}}
    Scale::Union{Vector{Float64}, Vector{Int64}}
    Power::Matrix{Float64}
    Phase::Matrix{Float64}
    nc::Int64
    nr::Int64
    dt::Union{Float64, Int64}
    dj::Union{Float64, Int64}
    lowerPeriod::Union{Float64, Int64}
    upperPeriod::Union{Float64, Int64}
    Power_avg::Union{Vector{Float64}, Nothing} = nothing
    Power_pval::Union{Matrix{Float64}, Nothing} = nothing
    Power_avg_pval::Union{Vector{Float64}, Nothing} = nothing
    COI::Union{Coi, Nothing} = nothing
    series_sim::Union{Matrix{Float64}, Nothing} = nothing
    Ridge::Union{Matrix{Float64}, Nothing} = nothing
    sPower::Union{Matrix{Float64}, Nothing} = nothing
end

@with_kw struct Wc
    Wave::Matrix{Complex}
    sWave::Matrix{Complex}
    Coherence::Matrix{Float64}
    Coherency::Matrix{Complex}
    Period::Union{Vector{Float64}, Vector{Int64}}
    Scale::Union{Vector{Float64}, Vector{Int64}}
    Power::Matrix{Float64}
    nc::Int64
    nr::Int64
    Wt_x::Wt
    Wt_y::Wt
    Power_avg::Union{Vector{Float64}, Nothing} = nothing
    Power_pval::Union{Matrix{Float64}, Nothing} = nothing
    Power_avg_pval::Union{Vector{Float64}, Nothing} = nothing
    COI::Union{Coi, Nothing} = nothing
    series_sim_x::Union{Matrix{Float64}, Nothing} = nothing
    series_sim_y::Union{Matrix{Float64}, Nothing} = nothing
    Ridge::Union{Matrix{Float64}, Nothing} = nothing
end

function blend_wt_with_p(
    x::Wt,
    Power_avg::Vector{Float64},
    Power_pval::Union{Matrix{Float64}, Nothing},
    Power_avg_pval::Union{Matrix{Float64}, Nothing},
    COI::Coi,
    series_sim::Union{Matrix{Float64}, Nothing})

    Wt(
        Wave           = x.Wave,        Ampl          = x.Ampl,
        Period         = x.Period,      Scale         = x.Scale,
        Power          = x.Power,       Phase         = x.Phase,
        nc          = x.nc,         nr             = x.nr,
        dt          = x.dt,         dj             = x.dj,
        lowerPeriod = x.lowerPeriod,
        upperPeriod = x.upperPeriod,
        Power_avg   = Power_avg,
        Power_pval  = Power_pval,   Power_avg_pval = Power_avg_pval,
        COI         = COI       ,   series_sim     = series_sim,
        sPower      = x.sPower)
end

function blend_wt_with_sPower(
    x::Wt, sPower::Matrix{Float64})
    Wt(
        Wave           = x.Wave,         Ampl          = x.Ampl,
        Period         = x.Period,       Scale         = x.Scale,
        Power          = x.Power,        Phase         = x.Phase,
        nc             = x.nc,           nr            = x.nr,
        dt             = x.dt,           dj            = x.dj,
        lowerPeriod    = x.lowerPeriod,
        upperPeriod    = x.upperPeriod,
        Power_avg      = x.Power_avg,
        Power_pval     = x.Power_pval,
        Power_avg_pval = x.Power_avg_pval,
        COI            = x.COI,          series_sim    = x.series_sim,
        Ridge          = x.Ridge,
        sPower         = sPower)
end

function blend_wc_with_p(
    x::Wc,
    Wt_x::Wt,
    Wt_y::Wt,
    Power_avg::Vector{Float64},
    Power_pval::Union{Matrix{Float64}, Nothing},
    Power_avg_pval::Union{Matrix{Float64}, Nothing},
    COI::Coi,
    series_sim_x::Union{Matrix{Float64}, Nothing},
    series_sim_y::Union{Matrix{Float64}, Nothing})

    Wc(
        Wave         = x.Wave,       sWave          = x.sWave,
        Coherence    = x.Coherence,  Coherency      = x.Coherency,
        Period       = x.Period,     Scale          = x.Scale,
        Power        = x.Power,
        nc           = x.nc,         nr             = x.nr,
        Wt_x         = Wt_x,         Wt_y           = Wt_y,
        Power_avg    = Power_avg,
        Power_pval   = Power_pval,   Power_avg_pval = Power_avg_pval,
        COI          = COI,
        series_sim_x = series_sim_x, series_sim_y   = series_sim_y)
end
