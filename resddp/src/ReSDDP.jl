module ReSDDP
    using JuMP
    using MathOptInterface
    using HiGHS
    using Clp
    using HDF5
    
    include("./structures.jl")
    include("./reserves.jl") #Added
    include("./reademps.jl")
    include("./aggregate.jl")
    include("./readinput.jl")
    include("./sampletree.jl")
    include("./solvefwd.jl")
    include("./solvebwd.jl")
    include("./stageprob.jl")
    include("./stageprob_full.jl")
    include("./stageprob_det.jl")
    include("./costprob.jl")
    include("./hprob.jl")

    include("./load.jl")
    include("./init.jl")
    include("./feasibility.jl")
    include("./train.jl")
    include("./simulate.jl")
    include("./save.jl")
    include("./print.jl")

    function print(model::Model, parameters::Parameters,LFeasibility,LStrategy)
        println("  ")
        println("Number of hydro systems           :     ",model.NHSys)
        println("Number of H2 areas                :     ",model.H2Data.NArea)
        println("Number of time steps              :     ",parameters.Time.NK)
        if LStrategy
            println("Number of SDDP forward scenarios  :     ",parameters.Control.NScen)
            println("Number of SDDP stages             :     ",parameters.Control.NStage)
            println("Residual function approximation   :     ",parameters.Control.LCostApprox)
            if parameters.Control.LCostApprox
                if parameters.Control.LCostApproxNewCuts
                    println("..new cost cuts per scenario      :     ",parameters.Control.CCMaxIter)
                else
                    println("..using old cost cuts")
                end
            end
            println("VRES uncertainty                  :     ",parameters.Control.LWindStoch)
            if parameters.Control.LWindStoch
                println("..number of scenarios             :     ",parameters.Control.NWindScen)
            end
        end
        if LFeasibility
            println("Feasibility cuts in optimization  :     ",parameters.Control.LFeasSpace)
        end
        println("Number of scenarios in simulation :     ",parameters.Control.NScenSim)
        println("Number of stages in simulation    :     ",parameters.Control.NStageSim)
        if parameters.Control.LIgnoreCrossCorr
            println("Inflow cross correlation neglected")
        end
        println(" ")
    end

    # Export functions for loading, training and simulating model
    export load, load_inflow, load_feas_cut, feasibility, train!, simulate_detailed, simulate_aggregated, init_strategy, extend_strategy, init_system

    # Export structures for model, strategy, scenario, results and parameters
    export Model, Strategy, Scenario, Result, Parameters
    
    # Export functions for printing results
    export print_results, print_results_h5, print_detailed_results, print_detailed_results_h5, print_dims, print_strategy, print_feas
end
