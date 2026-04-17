using Revise # Sikrer at oppdateringer av kildekode kompileres automatisk
using Pkg
Pkg.activate(@__DIR__) # Aktiverer lokalt pakkemiljø
# Pkg.resolve()
Pkg.instantiate()

using ReSDDP
include("params.jl") # Import parameters

import JuMP
import CPLEX
optimizer = JuMP.optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPX_PARAM_SCRIND" => 0,
    "CPX_PARAM_PREIND" => 1,
    "CPX_PARAM_LPMETHOD" => 2 #1=primal, 2=dual, 3=network, 4=barrier
)

println("Threads available: ",Threads.nthreads())

datapath = "C:\\datasets\\RES100_tutorial\\4area_test\\"
#datapath = "C:\\Users\\arildh\\data\\res100\\norge30_test\\"
#datapath = "C:\\Users\\arildh\\data\\res100\\HydroCen_LowEmission_V10\\"

model = load(datapath, parameters) 
inflow_model = load_inflow(datapath, model, parameters)

using JLD2 
using FileIO 

#Compute feasibility cuts 
println("Compute feasibility cuts..")
feas_spaces = feasibility(model, inflow_model, parameters, datapath; optimizer=optimizer)

# Save feasibility cuts to file
file = File(format"JLD2", joinpath(datapath, "feas_spaces.jld2"))
save(file, "feas_spaces", feas_spaces)

# Load feasibility cuts from file
data = JLD2.load(file) 
feas_spaces = data["feas_spaces"]

ReSDDP.print(model, parameters, true, true)

strategy = init_strategy(model, parameters)
init_val = init_system(model, parameters)

#Compute strategy by SDDP
println("Start strategy computation..")

train!(strategy, init_val, model, inflow_model, feas_spaces, parameters; optimizer=optimizer)
# using Serialization
# serialize(joinpath(@__DIR__, "strategy.jls"), strategy) # Save cuts to file
# strategy = deserialize(joinpath(@__DIR__, "strategy.jls")) # Load cuts from file

# Save strategy to file
file = File(format"JLD2", joinpath(datapath, "strategy.jld2"))
save(file, "strategy", strategy)

# Load strategy from file
data = JLD2.load(file) 
strategy = data["strategy"]

# Simulate aggregated
println("Start simulation ..")
results_agg = simulate_aggregated(model, inflow_model, parameters, strategy, feas_spaces, init_val; optimizer = optimizer)

# Print results to ASCII files 
println("Write results ..")
print_results(datapath,results_agg,model,parameters)
print_results_h5(datapath,results_agg,model,parameters)

print_dims(datapath,model.NHSys,model.H2Data.NArea,parameters.Control.NStage,parameters.Control.NScenSim,strategy.NCut,parameters.Control.MaxIter,parameters.Control.CCMaxIter)
print_strategy(datapath,strategy,parameters.Control.LCostApprox)
print_feas(datapath,feas_spaces[1],model.NHSys)

println("Start detailed simulation ..")
results_det = simulate_detailed(model, inflow_model, parameters, strategy; optimizer = optimizer)

println("Write detailed results ..")
print_detailed_results(datapath,results_det,model,parameters)
print_detailed_results_h5(datapath,results_det,model,parameters)
println("Program finished.")