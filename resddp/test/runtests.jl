using ReSDDP
using Test
using TestReports
using JLD2
using FileIO

@testset "ReSDDP.jl" begin
    @testset "4-area convergence" begin
        include("data/convergence_test/params.jl") # Load parameters
        model_file = File{format"JLD2"}("data/convergence_test/4area_deterministic.jld2")
        data = JLD2.load(model_file)
        model = data["model"]
        inflow_model = data["inflow_model"]
        feas_spaces = data["feas_spaces"]

        strategy = init_strategy(model, parameters)
        init_val = init_system(model, parameters)

        (sum_forward, sum_backward) = train!(strategy, init_val, model, inflow_model, feas_spaces, parameters)
        ϵ = 1e-9
        for i in eachindex(sum_backward)
            if i > 1
                # Check that sumback is increasing
                @test sum_backward[i] + ϵ > sum_backward[i-1]
            end
            # Check that SumForw is always greater than SumBack
            @test sum_forward[i] + ϵ > sum_backward[i]
        end
        # Check that case has converged
        @test sum_forward[end] ≈ sum_backward[end] atol=ϵ

    end
end
