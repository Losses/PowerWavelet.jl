using Random, RCall, Test
@rimport WaveletComp as rWaveletComp;

function rnd(x, digits = 5)
    round.(x, digits = digits)
end

function check(x, y, digits = 5)
    sum(abs.(rnd(x, digits) .- rnd(y, digits))) == 0
end

Random.seed!(2333);

x = rnd_signal(1000);
y = rnd_signal(1000);

r_wc = rcopy(rWaveletComp.wc(x, y, var"make.pval" = false));
j_wc = wc(x, y, make_pval = false);

@testset "PowerWavelet.jl" begin
    ############################################################
    # Test Wave
    ############################################################

    @test check(r_wc[:Wave_x],  j_wc.Wt_x.Wave)
    @test check(r_wc[:Wave_y],  j_wc.Wt_y.Wave)
    @test check(r_wc[:Wave_xy], j_wc.Wave)

    @test check(r_wc[:Power_x],  j_wc.Wt_x.Power)
    @test check(r_wc[:Power_y],  j_wc.Wt_y.Power)
    @test check(r_wc[:Power_xy], j_wc.Power)

    @test check(r_wc[:sPower_x],  j_wc.Wt_x.sPower)
    @test check(r_wc[:sPower_y],  j_wc.Wt_y.sPower)
    @test check(r_wc[:Power_xy], j_wc.Power)

    @test check(r_wc[:Power_x_avg],  j_wc.Wt_x.Power_avg)
    @test check(r_wc[:Power_y_avg],  j_wc.Wt_y.Power_avg)
    @test check(r_wc[:Power_xy_avg], j_wc.Power_avg)

    @test check(r_wc[:Coherence], j_wc.Coherence)
    @test check(r_wc[:Coherency], j_wc.Coherency)

    ############################################################
    # Test Period and Scale related
    ############################################################

    @test check(r_wc[:Period], j_wc.Period)
    @test check(r_wc[:Scale], j_wc.Scale)

    ############################################################
    # Something else
    ############################################################

    @test check(r_wc[:nc], j_wc.nc)
    @test check(r_wc[:nr], j_wc.nr)

    @test check(r_wc[:axis_1], j_wc.COI.axis_1)
    @test check(r_wc[:axis_2], j_wc.COI.axis_2)
    @test check(r_wc[:coi_1], j_wc.COI.coi_1)
    @test check(r_wc[:coi_2], j_wc.COI.coi_2)
end
