# ReSDDP Extensions — Operating Reserves

This document describes the extensions added to the ReSDDP model as part of a master's thesis. All additions are marked with `#ADDED` in the source code.

## Overview

The main extension introduces **operating reserve constraints** into the SDDP optimization. When enabled, the model requires that each price zone holds sufficient upward and downward reserve capacity at every time step, supplied by hydropower and wind power.

---

## Enabling the Extension

Set `LOperatingReserves = true` in [`examples/params.jl`](examples/params.jl) before building the `Control` struct:

```julia
LOperatingReserves = true
CTR = ReSDDP.Control(..., LOperatingReserves, ...)
```

When `false` (the default), the model runs exactly as before — no reserve variables or constraints are added.

---

## File-by-File Summary

### `src/structures.jl`

**New structs:**

`ReserveZoneReq` — reserve requirements for one price zone:

| Field | Description |
|---|---|
| `zone_name` | Zone identifier (e.g. `"NO1"`) |
| `RI_up / RI_down` | Reference incident up/down [GW] |
| `NI_up / NI_down` | Normal imbalance factor for onshore wind |
| `NI_up_OWP / NI_down_OWP` | Normal imbalance factor for offshore wind |
| `MaxLoad` | Maximum load in zone [GW] |
| `owp_areas_in_zone` | Indices of offshore wind areas in this zone |

`OperatingReserves` — top-level reserve data object stored on `Model`:

| Field | Description |
|---|---|
| `NZ` | Total number of price zones |
| `price_zones` | Zone name strings |
| `zone_reqs` | `Vector{ReserveZoneReq}` |
| `area_to_zone` | Maps each model area index to a zone index |
| `areas_in_zone` | Maps each zone to its list of model area indices |
| `hydrosys_to_area` | Maps each aggregated hydro system to its model area |
| `LMarkReserves` | Flag for market-based reserves |
| `a`, `b` | Empirical ENTSO-E constants for the load-dependent reserve term |

**Extended structs:**

`Control` — new field `LOperatingReserves::Bool`.

`Model` — new field `ORData::OperatingReserves`.

`Result` and `DetailedResult` — new output tables (all added together):

| Table | Description |
|---|---|
| `CapZoneUpTable` | Upward reserve capacity provided per zone [GWh] |
| `CapZoneDownTable` | Downward reserve capacity provided per zone [GWh] |
| `HydroCapUpTable` | Hydro contribution to upward reserves |
| `HydroCapDownTable` | Hydro contribution to downward reserves |
| `WindCapDownTable` | Wind contribution to downward reserves |
| `CapDualUpTable` | Shadow price on upward reserve constraint |
| `CapDualDownTable` | Shadow price on downward reserve constraint |
| `ObjTable` | Stage-wise operational cost (objective minus future-cost alpha) |
| `WaterValueTable` | Water values (shadow price on reservoir state constraint) |
| `MarkCapUpTablePos/Neg` | Market reserve capacity contributions (positive/negative) |
| `MarkCapDownTablePos/Neg` | Market downward reserve contributions |
| `SlackUpTable / SlackDownTable` | Reserve slack (DetailedResult only) |

---

### `src/stageprob_full.jl` and `src/stageprob_det.jl` 

When `LOperatingReserves = true`, the following is added inside each stage problem build function.

**New decision variables:**

| Variable | Description |
|---|---|
| `cap_zone_up[z,k]` | Total upward reserve in zone `z` at step `k` |
| `cap_zone_down[z,k]` | Total downward reserve in zone `z` at step `k` |
| `cap_hydro_up[sys,k]` | Hydro upward reserve (aggregated model) |
| `cap_hydro_down[sys,k]` | Hydro downward reserve (aggregated model) |
| `cap_hydro_up_mod[sys,mod,k]` | Hydro upward reserve per module (detailed model) |
| `cap_hydro_down_mod[sys,mod,k]` | Hydro downward reserve per module (detailed model) |
| `cap_wind_down[area,k]` | Wind downward reserve |
| `wp_avail[area,k]` | Available wind production (RHS updated per scenario) |
| `slackUp[z,k]` | Slack on upward reserve requirement |
| `slackDown[z,k]` | Slack on downward reserve requirement |

**Reserve requirement expressions** (`cap_up_amount`, `cap_down_amount`) implement the ENTSO-E formula:

```
R_req[z,k] = RI * DT
           + NI_onshore * sum(wp_avail[a,k] for onshore a in zone z)
           + NI_offshore * sum(wp_avail[a,k] for offshore a in zone z)
           + (√(a · MaxLoad[z] / DT + b²) − b) · DT
```

where `DT` is the time-step length in hours, `RI` is the reference incident, `NI` is the nominal imbalance factor, and the last term is a load-dependent stochastic component.

**New constraints:**

| Constraint | Description |
|---|---|
| `reserve_req_up/down` | Zone capacity + slack ≥ required amount |
| `reserve_split_up/down` | Zone total = sum of hydro + wind contributions |
| `hydro_up` | Hydro upward reserve ≤ remaining capacity above current dispatch |
| `hydro_down` | Hydro downward reserve ≤ current dispatch above minimum |
| `hydroRes_cap_up` | Hydro upward reserve ≤ current reservoir level |
| `wind_dn` | Wind downward reserve ≤ current wind dispatch |

Slack variables are penalized in the objective using coefficient `CRat`, so the model prefers to meet reserve requirements but is not infeasible if unable.

---

### `src/simulate.jl`

- `LOperatingReserves` and `ORData` are passed to the stage problem builder in both `simulate` and `simulate_detailed`
- When reserves are active, `wp_avail_fix` is updated each scenario/stage with the realized wind output, making reserve requirements dynamic
- `NZ` is passed to `init_detailed_result` to correctly size reserve output tables
- Bug fix: `max(0.0, ResInit + CurrInf)` prevents negative initial reservoir levels being passed to the solver

---

### `src/save.jl`

`save!` and `save_detailed!` were extended to extract reserve results after each solve:

- Zone reserve capacities (`cap_zone_up/down`)
- Shadow prices on reserve constraints (`reserve_req_up/down`)
- Hydro and wind reserve contributions
- Market reserve contributions (if available)
- Stage-wise objective value (excluding future-cost term `alpha`)
- Water values as shadow prices on the reservoir state or end-volume constraints

---

### `src/print.jl`

When `LOperatingReserves = true`, a group `XNordic_Reserve_req` is written to the HDF5 output file containing:

- Zone names and the area-to-zone mapping
- Per zone: `ReserveUp`, `ReserveDown`, `ReserveDualUp`, `ReserveDualDown` (dimensions: NScen × NStage × NK)
- Per area within each zone: hydro, wind, and market reserve contributions
- Water values are written alongside existing hydro results

---
## Data Flow

```
params.jl
  └─ LOperatingReserves = true
       │
       ▼
load.jl → ReadOperatingReserves() → ORData stored on Model
       │
       ▼
simulate / train
  └─ StageProbDet.Build(..., LOperatingReserves, ORData)
       └─ Adds reserve variables + constraints to JuMP model
       └ ─ Updates wp_avail RHS from realized wind
       │
       ▼
save! / save_detailed!
  └─ Extracts reserve quantities and duals into Result tables
       │
       ▼
print_results / print_detailed_results
  └─ Writes XNordic_Reserve_req group to HDF5 output
```
