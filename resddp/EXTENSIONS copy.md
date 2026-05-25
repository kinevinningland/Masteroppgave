# ReSDDP Extensions — Operating Reserves

This document describes all extensions added to the ReSDDP model as part of a master's thesis.
All additions are marked with `ADDED` and modifications to existing code are marked with `CHANGED`.

## Overview

The main extension introduces **operating reserve constraints** into both aggregated and detailed hydropower simulation.

When enabled `LOperatingReserves = true`, the model requires that each price zone holds
sufficient upward and downward reserve capacity at every time step, supplied by hydropower,
wind power, and (optionally) market steps. The reserve requirements follow the ENTSO-E
methodology.

When `false` (the default), the model runs exactly as before — no reserve variables or
constraints are added.


## Enabling the Extension

Set `LOperatingReserves = true` in [`examples/params.jl`](examples/params.jl) 

## Input Data Files 

The following files can be present in the dataset folder when `LOperatingReserves = true`:

| File | Description |
| `ORData_zone_reqs.csv` | One row per price zone: `zone_name, RI_up, RI_down, NI_up, NI_down, NI_up_OWP, NI_down_OWP` |
| `ORData_area_zones.csv` | One row per area: `area_index, zone_name` |
| `ORData_mark.csv` | Reserve market eligibility: `power_cat, zone_name, pos (true/false), neg (true/false)` |

`ORData_area_zones.csv` - Maps each area index to a price zone. Areas within the same zone jointly satisfy the reserve requirement. Areas not listed in the file are excluded from reserve accounting.

`ORData_zone_reqs.csv` - Specifies the reserve requirement parameters for each price zone: the reference incident, the nominal imbalance factor for upward and downward reserves for both onshore and offshore wind.

If `ORData_mark.csv` - Optional file that allows market steps to contribute to reserves. Each row specifies a power category and zone, along with boolean flags indicating whether positive and negative market steps of that category are eligible for upward and downward reserve provision.


## Data Flow

```
params.jl
  └─ LOperatingReserves = true
       │
       ▼
load.jl → ReadOperatingReserves()   [in readinput.jl]
       │   reads: ORData_zone_reqs.csv, ORData_area_zones.csv,
       │           ORData_mark.csv, model.h5
       │   → ORData stored on Model
       ▼
solvefwd.jl / solvebwd.jl   (SDDP training)
  └─ StageProbFull.Build(..., LOperatingReserves=false, ORData)
       └─ Reserves NOT active during training
       │
       ▼
simulate_detailed / simulate_aggregated
  └─ StageProbDet/Full.Build(..., LOperatingReserves=true, ORData)
       └─ Adds reserve variables + constraints to JuMP model
       └─ Updates wp_avail RHS from realized wind each scenario/stage
       │
       ▼
save! / save_detailed!
  └─ Extracts reserve quantities, duals, water values, H2 results into Result tables
       │
       ▼
print_results_h5 / print_detailed_results_h5
  └─ Writes XNordic_Reserve_req group to HDF5 output
  └─ Writes WaterValue, ObjectiveValue, H2 results to HDF5 output
```

## Issue fix - modifications to reademps.jl

An issue was identified in how hydropower areas are constructed in the detailed simulation. The dataset contains 24 hydropower areas. Rather than creating hydro variables for the correct price area IDs, the code created variables for the first 24 (lenght of the hydropower areas) price areas by index.

To fix the area indexing issue, the detailed simulation was extended to create hydropower variables for all 78 price areas rather than only the 
first 24. This ensures that all areas with hydropower are correctly included. This was done by insted of iterating throught the size of `NHSys` for hydropower variables in `Stageprob_det`, they are iterating trough `NArea`. 

In the original framework, characteristics for hydropower, batteries, and hydrogen are all stored in the struct `AHData`, which is used in 
the detailed simulation when hydropower variables are created. This struct holds detailed module-level characteristics for each area. In the aggregated simulation, `AHData` is used to construct `HSys`, which stores aggregated hydropower characteristics only, excluding hydrogen and batteries. Hydrogen characteristics are separately extracted from `AHData` and stored in `H2Data` for the areas designated as hydrogen areas in the dataset.

However, creating hydropower variables for hydrogen areas is 
undesirable, as hydrogen storage is handled through separate variables which is added though this thesis work. To resolve this, the way data is read into the optimization problem was modified. Rather than reading hydrogen characteristics into `AHData`, dummy objects are now assigned to hydrogen areas in that struct:

```julia
if AreaName[iArea] in MyKeys && !contains(AreaName[iArea], "H2")
```

The construction of `H2Data` in `ReadH2`was also changed to read hydrogen characteristics directly from the input dataset instead of extracting them from `AHData`:

Since hydrogen areas are now assigned dummy objects in `AHData`, the original approach of reading `H2MaxDis` and `H2MaxRes` from `AHData[iArea].MData[1]` is no longer valid. The updated `ReadH2` therefore opens the HDF5 input file directly and reads module-level characteristics from `hydro_data/<AreaName>/Module_data`. The `AHData` argument was removed from the function signature entirely, as the function no longer depends on the pre-processed hydropower struct.

This ensures that hydrogen storage is represented only through the dedicated hydrogen variables and not incorrectly duplicated through the hydropower variables.

A similar adjustment was required in `ReadDynmod`, which reads time-dependent module data from the binary `DYNMODELL.OM` files. In the original implementation, all areas used `AHData[iArea].MData[iMod].ModId` to look up the correct module index in the file via findfirst, ensuring that data was mapped to the right module regardless of storage order. Since hydrogen areas are no longer stored in `AHData`, this lookup would fail for those areas. The updated function therefore branches on whether the current area is an H2 area: 

* for hydrogen areas, module data is copied directly in file order, relying on the assumption that the module ordering in DYNMODELL.OM is consistent with that in model.h5
* for all other areas, the original index-lookup logic is preserved unchanged

In simulate_detailed in `simulate.jl`, multiple loop bounds were changed from `NHSys` to `NArea`, marked by #Changed from NHSys to NArea. This was done to fix the indexing issue mentioned. Another change to fix this problem was adding this line:
```julia
iArea = model.HSys[iSys].AreaNo 
```
to loops that iterate over hydropower systems (`iSys = 1:model.NHSys`) and access `AHData` when updating the `cut`constraint. Previously, `AHData` was indexed directly by `iSys`, implicitly assuming that the hydropower system index corresponds to the area index. This assumption no longer holds after the area indexing fix, since `HSys` contains only hydropower systems and their indices do not necessarily match the area indices in `AHData`.


## File-by-File Summary


### `src/structures.jl`
New struct added: `ReserveZoneReq` — reserve requirements for a price zone:

| Field | Type | Description 
| `zone_name` | `String` | Zone identifier
| `RI_up` | `Float64` | Reference incident upward [GW] 
| `RI_down` | `Float64` | Reference incident downward [GW] 
| `NI_up` | `Float64` | Normal imbalance factor for onshore wind, upward 
| `NI_down` | `Float64` | Normal imbalance factor for onshore wind, downward 
| `NI_up_OWP` | `Float64` | Normal imbalance factor for offshore wind, upward 
| `NI_down_OWP` | `Float64` | Normal imbalance factor for offshore wind, downward 
| `MaxLoad` | `Float64` | Maximum load in zone (computed from `AMData`) 
| `owp_areas_in_zone` | `Vector{Int}` | Indices of offshore wind areas in this zone 

New struct added: `OperatingReserves` — top-level reserve data object stored in `Model`:

| Field | Type | Description 
| `NZ` | `Int` | Total number of price zones 
| `price_zones` | `Vector{String}` | Zone name strings 
| `zone_reqs` | `Vector{ReserveZoneReq}` | Per-zone requirements 
| `area_to_zone` | `Vector{Int}` | Maps each model area index to a zone index
| `areas_in_zone` | `Vector{Vector{Int}}` | Maps each zone index to its list of model area indices 
| `hydrosys_to_area` | `Vector{Int}` | Maps each aggregated hydro system index to its model area 
| `LMarkReserves` | `Bool` | Flag enabling market-step reserve contributions 
| `a` | `Int` | Empirical ENTSO-E constant for load-dependent reserve term (= 10) 
| `b` | `Int` | Empirical ENTSO-E constant for load-dependent reserve term (= 150) 
| `pos_by_area` | `Dict{Int, Set{Int}}` | Per-area set of positive-capacity market step indices eligible for reserves 
| `neg_by_area` | `Dict{Int, Set{Int}}` | Per-area set of negative-capacity market step indices eligible for reserves 

Added field in existing struct `Control`:
* LOperatingReserves::Bool

Added field in existing struct `Model`:
* ORData::OperatingReserves

Added field in existing struct `Result` (all below `H2DisTable`):

| Field | Dimensions | Description 
| `CapZoneUpTable` | `NZ × NScen × NStage × NK` | Total upward reserve provided per zone 
| `CapZoneDownTable` | `NZ × NScen × NStage × NK` | Total downward reserve provided per zone 
| `HydroCapDownTable` | `NHSys × NScen × NStage × NK` | Hydro downward reserve
| `HydroCapUpTable` | `NHSys × NScen × NStage × NK` | Hydro upward reserve  
| `WindCapDownTable` | `NArea × NScen × NStage × NK` | Wind downward reserve 
| `CapDualUpTable` | `NZ × NScen × NStage × NK` | Shadow price on upward reserve constraint 
| `CapDualDownTable` | `NZ × NScen × NStage × NK` | Shadow price on downward reserve constraint 
| `ObjTable` | `NScen × NStage` | Stage-wise cost (objective minus future-cost `alpha`) 
| `WaterValueTable` | `NHSys × NScen × NStage` | Shadow price on reservoir state constraint 
| `MarkCapUpTablePos` | `NArea × NScen × NStage × NK` | Positive-step market upward reserve 
| `MarkCapDownTablePos` | `NArea × NScen × NStage × NK` | Positive-step market downward reserve 
| `MarkCapUpTableNeg` | `NArea × NScen × NStage × NK` | Negative-step market upward reserve 
| `MarkCapDownTableNeg` | `NArea × NScen × NStage × NK` | Negative-step market downward reserve 

Same fields added in existing struct `DetailedResult`, plus two additional:

| Field | Dimensions | Description 
| `SlackUpTable` | `NZ × NScen × NStage × NK` | Slack on upward reserve constraint 
| `SlackDownTable` | `NZ × NScen × NStage × NK` | Slack on downward reserve constraint 

Note: `HydroCapUpTable` and `HydroCapDownTable` in `DetailedResult` have an extra module dimension: `NHSys × NMaxMod × NScen × NStage × NK`.

### `src/readinput.jl`

Added function `ReadOperatingReserves(...)`:

Called from `load.jl`. Reads three CSV files and `model.h5` to build the `OperatingReserves`
object. Returns a dummy empty struct immediately when `LOperatingReserves = false`.

Logic:
- Reads `ORData_zone_reqs.csv` — one row per zone with RI and NI values
- Reads `ORData_area_zones.csv` — Maps area indices to price zones. Areas in the same zone jointly satisfy the reserve requirement. Listing an area is optional — unlisted areas are excluded
- Computes `MaxLoad` per zone as the maximum of the sum of all load series in that zone
- Identifies offshore wind areas by `-OWP` suffix in `AreaName`
- Reads `ORData_mark.csv` — defines which market power categories can provide reserves (by zone, pos/neg flag)
- Reads `model.h5 / market_data / power_type` to match market step names to area indices
- Builds `pos_by_area` and `neg_by_area`: local `iMark` indices per area eligible for reserve provision, separated by the sign of their capacity
- Safety check: no `iMark` may appear in both sets for the same area
- Sets empirical ENTSO-E constants: `a = 10`, `b = 150`


### `src/stageprob_full.jl` (aggregated model)

`LOperatingReserves` and `ORData` added as input arguments to the Build function, allowing the stage problem to conditionally include reserve variables and constraints based on the `LOperatingReserves` flag.

Entire reserve block guarded by are added `if LOperatingReserves`.

New decision variables:

| Variable | Description 
| `cap_zone_up[z,k]` | Total upward reserve in zone `z` at step `k` 
| `cap_zone_down[z,k]` | Total downward reserve in zone `z` at step `k` 
| `cap_hydro_up[iSys,k]` | Hydro upward reserve per aggregated system 
| `cap_hydro_down[iSys,k]` | Hydro downward reserve per aggregated system 
| `cap_wind_down[iArea,k]` | Wind downward reserve 
| `cap_mark_up_pos[a,iMark,k]` | Market upward reserve from positive-capacity steps 
| `cap_mark_down_pos[a,iMark,k]` | Market downward reserve from positive-capacity steps 
| `cap_mark_up_neg[a,iMark,k]` | Market upward reserve from negative-capacity steps 
| `cap_mark_down_neg[a,iMark,k]` | Market downward reserve from negative-capacity steps 
| `wp_avail[a,k]` | Available wind (RHS overwritten per scenario in `simulate.jl`) 

Reserve requirement expressions using the ENTSO-E formula:
```julia
cap_up_amount[z,k] = RI_up · DT
    + NI_up     · Σ wp_avail[a,k]  for onshore a in zone z
    + NI_up_OWP · Σ wp_avail[a,k]  for offshore a in zone z
    + (√(a · MaxLoad[z] / (MW2GW · DT) + b²) − b) · MW2GW · DT
```
Same form for `cap_down_amount` using `RI_down`, `NI_down`, `NI_down_OWP`.

New constraints:

| Constraint | Description 
| `reserve_req_up[z,k]` | `cap_zone_up[z,k] ≥ cap_up_amount[z,k]`, Ensures that the total upward reserve capacity provided in zone z at time step k meets or exceeds the required amount
| `reserve_req_down[z,k]` | `cap_zone_down[z,k] ≥ cap_down_amount[z,k]`, Ensures that the total downward reserve capacity provided in zone z at time step k meets or exceeds the required amount
| `reserve_split_up[z,k]` | The total upward reserve in a zone equals the sum of upward reserve contributions from all hydro systems and eligible market steps within that zone.
| `reserve_split_down[z,k]` | The total downward reserve in a zone equals the sum of downward reserve contributions from all hydro systems, wind and eligible market steps within that zone.
| `hydro_up[iSys,k]` | The sum of current hydro production and upward reserve capacity for a system cannot exceed its maximum production capacity.
| `hydro_dn[iSys,k]` | Current hydro production minus downward reserve capacity must remain at or above the minimum production level, ensuring the unit can actually reduce output by the reserved amount.
| `hydroRes_cap_up[iSys,k]` | The upward reserve offered by a hydro system cannot exceed its current reservoir level, ensuring there is enough stored water to follow through on the upward regulation.
| `wind_dn[iArea,k]` | Wind downward reserve cannot exceed actual wind production, so only power already being generated can be offered as downward regulation.
| `mark_up_pos / mark_dn_neg / mark_up_neg / mark_dn_pos` | Limits the reserve capacity that each market step can offer, based on its remaining headroom relative to its capacity bound.
| `wp_avail_fix[a,k]` | Fixes the available wind variable to zero by default. The right-hand side is overwritten with the realized wind output for each scenario and stage during simulation, making the reserve requirement wind-dependent.



### `src/stageprob_det.jl` (detailed model)

`LOperatingReserves`, `ORData` and `H2Data` are added as input arguments to the Build function, allowing the stage problem to conditionally include reserve variables and constraints based on the `LOperatingReserves` flag.

H2 variables are added (always present, regardless of `LOperatingReserves`):

| Variable | Description 
| `h2dis[iArea,k]` | H2 discharge [GWh/step] 
| `h2chg[iArea,k]` | H2 charge [GWh/step] 
| `h2res[iSys,k]` | H2 reservoir level [GWh] 
| `h2init[iArea]` | Initial H2 level (RHS set in simulate) 

H2 constraints added:

| Constraint | Description 
| `h2storage0[iArea, k=1]` | Hydrogen storage balance at the first time step where the reservoir level equals the initial storage minus net charging losses plus discharge.
| `h2storage[iArea, k=2:NK]` | Hydrogen storage balance for all subsequent time steps the reservoir level at each step equals the level from the previous step, adjusted for charging losses and discharge.
| `h2state[iArea]` | Sets the initial hydrogen storage level to zero by default; the right-hand side is overwritten with the actual carry-over storage from the previous stage during simulation


Power balance constraints `pbalHyd` and `pbalTerm` are removed and `pbal` are added:
H2 charge/discharge terms added to `pbalHyd` and it is extended to iterate over NArea instead of NHSys:
```julia
- (H2Data.Ind[iArea] > 0 ? h2chg[H2Data.Ind[iArea],k] : 0.0)   #ADDED
+ (H2Data.Ind[iArea] > 0 ? h2dis[H2Data.Ind[iArea],k] : 0.0)   #ADDED
```

Slack variables are added (always declared, penalised only when `LOperatingReserves = true`) and added to the objective function:
```julia
+ (LOperatingReserves ?
       sum(CNS.CRat·slackUp[z,k]   for z=1:ORData.NZ for k=1:NK)
     + sum(CNS.CRat·slackDown[z,k] for z=1:ORData.NZ for k=1:NK)
   : 0.0)   #ADDED
```
A reserve block guarded by `if LOperatingReserves` is added (module-level hydro resolution):

New decision variables:

| Variable | Description |
| `cap_zone_up[z,k]` | Total upward reserve in zone `z` at step `k` 
| `cap_zone_down[z,k]` | Total downward reserve in zone `z` at step `k` 
| `cap_hydro_up_mod[iArea,iMod,k]` | Hydro upward reserve per module 
| `cap_hydro_down_mod[iArea,iMod,k]` | Hydro downward reserve per module 
| `cap_wind_down[iArea,k]` | Wind downward reserve 
| `cap_mark_up_pos / down_pos / up_neg / down_neg` | Market step reserve variables (same as aggregated) 
| `wp_avail[a,k]` | Available wind (RHS overwritten in simulate) 

New constraints:

| Constraint | Description |
| `reserve_req_up[z,k]` | The total upward reserve capacity in a zone, plus a slack variable, must meet or exceed the required amount. The slack variable makes this a soft constraint, so the model remains feasible even if the requirement cannot be fully met, but slack is penalised in the objective.
| `reserve_req_down[z,k]` | Same as above for downward reserves.
| `reserve_split_up[z,k]` | The total upward reserve in a zone equals the sum of upward reserve contributions from all hydro modules and eligible market steps within that zone.
| `reserve_split_down[z,k]` | The total downward reserve in a zone equals the sum of downward reserve contributions from all hydro modules, wind areas, and eligible market steps within that zone.
| `hydro_up_mod[iArea,iMod,k]` | The upward reserve offered by a hydro module cannot exceed its remaining production headroom, i.e. the difference between its maximum discharge capacity and its current dispatch.
| `hydro_down_mod[iArea,iMod,k]` | The downward reserve offered by a hydro module cannot exceed its current power output, ensuring the unit can actually reduce production by the reserved amount.
| `hydro_res_up_mod[iArea,iMod,k]` | The upward reserve offered by a hydro module cannot exceed the energy equivalent of its current reservoir level, ensuring sufficient stored water is available to sustain the upward regulation.
| `wind_dn[iArea,k]` | Wind downward reserve cannot exceed actual wind production, so only power already being generated can be offered as downward regulation.
| `wp_avail_fix[a,k]` | Fixes the available wind variable to zero by default. The right-hand side is overwritten with the realized wind output for each scenario and stage during simulation, making the reserve requirement wind-dependent.

The cut constraint is changed: H2 end-state term added:
```julia
alpha - Σ CCR·rstate - Σ CCH·h2res[iArea,end] ≥ 0   #h2res Added
```

### `src/load.jl`

Call to `ReadOperatingReserves` is added:
```julia
ORData = ReadOperatingReserves(dataset, NArea, NHSys, NAreaSys, AreaSys,
                               AMData, AreaName, CTR.LOperatingReserves)  #Added
```

`Model(...)` constructor call is changed to include: `ORData` inserted after `DRData`:
```julia
model = Model(AHData, AMData, DMData, H2Data, USMod, WPData, DRData, ORData, ...)  #ORData Added
```

`AHData`is removed from the H2Data call.

### `src/init.jl`

`init_result` signature: `NZ` added as parameter:

```julia
function init_result(NArea, NHSys, NMaxMStep, NScen, NStage, NK, NLine, NZ)::Result  #NZ Added
```

Reserve result arrays in `init_result` are added (all below H2 tables):

```julia
CapZoneUpTable    = zeros(Float64, NZ,    NScen, NStage, NK)   #ALL Below ADDED
CapZoneDownTable  = zeros(Float64, NZ,    NScen, NStage, NK)
HydroCapUpTable   = zeros(Float64, NHSys, NScen, NStage, NK)
HydroCapDownTable = zeros(Float64, NHSys, NScen, NStage, NK)
WindCapDownTable  = zeros(Float64, NArea, NScen, NStage, NK)
CapDualUpTable    = zeros(Float64, NZ,    NScen, NStage, NK)
CapDualDownTable  = zeros(Float64, NZ,    NScen, NStage, NK)
ObjTable          = zeros(Float64, NScen, NStage)
WaterValueTable   = zeros(Float64, NHSys, NScen, NStage)
MarkCapUpTablePos    = zeros(Float64, NArea, NScen, NStage, NK)
MarkCapDownTablePos  = zeros(Float64, NArea, NScen, NStage, NK)
MarkCapUpTableNeg    = zeros(Float64, NArea, NScen, NStage, NK)
MarkCapDownTableNeg  = zeros(Float64, NArea, NScen, NStage, NK)
```

`init_detailed_result` signature: `NZ` added as parameter. The Same reserve arrays added are added, plus:

```julia
SlackUpTable   = zeros(Float64, NZ, NScen, NStage, NK)
SlackDownTable = zeros(Float64, NZ, NScen, NStage, NK)
```

Note: `HydroCapUpTable` and `HydroCapDownTable` in the detailed result have an extra module
dimension: `NHSys × NMaxMod × NScen × NStage × NK`.


### `src/solvebwd.jl` and `src/solvefwd.jl`

`ORData` is passed to stage problem builder during SDDP training, with `LOperatingReserves` explicitly set to `false`:

```julia
StageProbFull.Build(..., model.H2Data, false, model.ORData, optimizer)
```

Operating reserves are intentionally disabled during backward/forward training passes and only
active during simulation.



### `src/reademps.jl`

An issue was identified in how hydropower areas are constructed in the detailed simulation. The dataset contains 24 hydropower areas. Rather than creating hydro variables for the correct price area IDs, the code created variables for the first 24 (lenght of the hydropower areas) price areas by index.

To fix the area indexing issue, the detailed simulation was extended to create hydropower variables for all 78 price areas rather than only the 
first 24. This ensures that all areas with hydropower are correctly included. Since batteries are modeled as hydropower in the current framework, areas containing only battery storage are included in the hydropower variables. However, creating hydropower variables for hydrogen areas is 
undesirable, as hydrogen storage is handled through separate variables.

To resolve this, the way data is read into the optimization problem was modified. Rather than reading hydrogen characteristics into `AHData` in `ReadHDF5`, dummy objects are now assigned to hydrogen areas in that struct by excluding these areas:

```julia
if AreaName[iArea] in MyKeys && !contains(AreaName[iArea], "H2")
```

The construction of `H2Data` in `ReadH2`was also changed to read hydrogen characteristics directly from the input dataset instead of extracting them from `AHData`:

 This ensures that hydrogen storage is represented 
only through the dedicated hydrogen variables and not incorrectly duplicated through the hydropower variables.


### `src/simulate.jl`

Canges in `simulate_detailed`:

`initial_values::InitialValues` is added as an explicit argument to `simulate_detailed`, allowing the function to receive the correct initial reservoir and hydrogen storage levels rather than relying on hardcoded or default values.

H2 trajectory tracking:
```julia
SimulatedH2Traj = zeros(Float64, model.H2Data.NArea, NScenSim, NStageSim)  
H2Init  = zeros(Float64, model.H2Data.NArea)                              
```
Initial H2 levels set from `initial_values.H2Init`, which are initialized from `ResInitFrac` and the maximum reservoir capacity of each H2 area:

```julia
H2Init[1:model.H2Data.NArea] = initial_values.H2Init[1:model.H2Data.NArea]
```
and H2 storage is updated:

```julia
H2Init[1:model.H2Data.NArea] = SimulatedH2Traj[1:model.H2Data.NArea,iScen,t-1,end]
```

H2 initial constraint RHS updated each stage/scenario:

```julia
for iH2a = 1:model.H2Data.NArea 
    JuMP.set_normalized_rhs(SP_FORW[:h2storage0][iH2a,1], H2Init[iH2a])
end
```

H2 end-state saved after each solve:

```julia
for iH2a = 1:model.H2Data.NArea   #Added
    SimulatedH2Traj[iH2a, iScen, t] = JuMP.value(SP_FORW[:h2res][iH2a, end])
end
```

`NZ` is added and passed to `init_detailed_result` and `init_result`:

```julia
DetailedResultTable = init_detailed_result(..., model.ORData.NZ)  #ADDED NZ
ResultTable         = init_result(..., model.ORData.NZ)           #NZ ADDED
```

`LOperatingReserves` and `ORData` are added and passed to stage problem builder in both
`simulate_detailed` and `simulate_aggregated`:

```julia
StageProbDet.Build(..., model.H2Data,
                   parameters.Control.LOperatingReserves, model.ORData, optimizer) 
```

`wp_avail_fix` RHS updated per scenario/stage** when reserves are active is added for both `simulate_detailed` and `simulate_aggregated`:

```julia
if parameters.Control.LOperatingReserves   
    JuMP.set_normalized_rhs(SP_FORW[:wp_avail_fix][iArea,k],
                            max(model.WPData[...], 0.0))
end
```

Wind scenario sampling is also changed in `simulate_aggregated` to sample the same wind scenarios as detailed simulation: each simulation scenario now uses a distinct wind year instead of always using year 1:

```julia
# SampledWindYears = fill(1, parameters.Control.NScenSim)   #changed to the next
SampledWindYears = collect(1:parameters.Control.NScenSim)   
...
wYear = SampledWindYears[iScen] 
```

Added bug fix preventing negative initial reservoir levels in `simulate_detailed`:

```julia
JuMP.set_normalized_rhs(SP_FORW[:resbalReg0][iArea,iMod],
    max(0.0, ResInit[i,iMod] + CurrInf))
#ADDED max(0.0,...) to avoid negative reservoir levels
```

### `src/save.jl`

Changed `save!` signature: `LOperatingReserves` and `ORData` added:

```julia
function save!(RT::Result, SP_FORW, AMData, H2Data, InflowSys,
               NArea, NHSys, NK, NLine, s, t,
               LOperatingReserves, ORData)   #ORData, LOperatingReserves 
```

Added water value extraction in `save!`:

```julia
RT.WaterValueTable[iSys,s,t] = JuMP.shadow_price(SP_FORW[:rstate][iSys]) 
```

Added reserve result extraction block in `save!`, guarded by `if LOperatingReserves` where the results used for the analysis of the thesis are extracted.


Changed `save_detailed!` signature: `LOperatingReserves`, `HSys` amd `H2Data` are added:

```julia
function save_detailed!(DRT::DetailedResult, SP_FORW, AMData, H2Data, AHData, NArea, NHSys, NK, NLine, s, t, LOperatingReserves, HSys) 
```

Added stage objective extraction in `save_detailed!`:

```julia
DRT.ObjTable[s,t] = JuMP.objective_value(SP_FORW) - JuMP.value(SP_FORW[:alpha])
```

Added water value extraction in `save_detailed!` (uses end-volume constraint):

```julia
DRT.WaterValueTable[i,s,t] = JuMP.shadow_price(SP_FORW[:endvol][iSys])  
```

Added H2 result extraction in `save_detailed!`, guarded by `H2Data.Ind[iArea] > 0`:

```julia
DRT.H2StoreTable[iArea,s,t,k] = JuMP.value(SP_FORW[:h2res][H2Data.Ind[iArea],k])
DRT.H2DisTable[iArea,s,t,k]   = -(1-CompLoss)*h2chg + h2dis
```

Added detailed reserve extraction block in `save_detailed!`, guarded by
`if LOperatingReserves`, where the results used for the analysis of the thesis are extracted.


### `src/print.jl`

The two functions that are changed are `print_results_h5` and `print_detailed_results_h5`to extract results for the analysis of the thesis. 

`ObjectiveValue` is added and written to HDF5 in both `print_results_h5` and `print_detailed_results_h5`:

```julia
write(file, "ObjectiveValue", RT.ObjTable)  
```

`WaterValue` is added and  written to HDF5 per hydro area in both print functions:

```julia
write(hydroGroup, "WaterValue", RT.WaterValueTable[hydro_idx,:,:])  
attrs(hydroGroup["WaterValue"])["Dim 1"] = "NScen"                  
attrs(hydroGroup["WaterValue"])["Dim 2"] = "NStage"                 
```

Changed market step HDF5 output: duplicate name handling added to avoid HDF5 dataset
name collisions when multiple market steps share the same name:

```julia
for iMark = 1:model.AMData[iArea].NMStep   #ADDED (replaces old loop)
    base_name = model.AMData[iArea].MSData[iMark].Name
    name = base_name
    if haskey(marketStepGroup, name)
        name = string(base_name, "_", iMark)
    end
    write(marketStepGroup, name, RT.MarkTable[iArea,iMark,:,:,:])
end
```

H2 results are added and  written to HDF5 in `print_detailed_results_h5`,
guarded by `model.H2Data.Ind[iArea] > 0`:

```julia
if model.H2Data.Ind[iArea] > 0  
    H2Group = create_group(areaGroup, "H2")
    write(H2Group, "Storage",   DRT.H2StoreTable[iArea,:,:,:])
    write(H2Group, "Discharge", DRT.H2DisTable[iArea,:,:,:])
end
```

Reserve HDF5 group is added and written when `LOperatingReserves = true` in both print functions. Group `XNordic_Reserve_req` 


Changed `print_results` text output: market sum uses `init=0.0` to handle areas with no market steps:

```julia
val = sum(RT.MarkTable[iArea,iMark,...] for iMark=1:...; init=0.0)  
@printf(out, "%16.6f ", val)                                         
```


