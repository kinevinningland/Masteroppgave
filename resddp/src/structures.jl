struct Constants
    MAGEFF2GWH::Float64
    GWH2MAGEFF::Float64
    MW2GWHWEEK::Float64
    MW2GW::Float64
    M3S2MM3::Float64
    CEnPen::Float64
    CResPen::Float64
    CRampPen::Float64
    CCapPen::Float64
    CAuxPen::Float64
    CSpi::Float64
    CByp::Float64
    CNegRes::Float64
    CMinRes::Float64
    CRat::Float64
    CFeas::Float64
    Big::Float64
    AlphaMax::Float64
    InfUB::Float64
    FeasTol::Float64
end
 
struct Aggregation
    ModCutoff::Int
    ResCutoff::Float64
    ProdCutoff::Float64
    DeplCutoff::Float64
end
 
struct Time
    NSecHour::Float64
    NHoursWeek::Float64
    NInflowYear::Int
    NWeek::Int
    NK::Int
    DT::Float64
    WeekFrac::Float64
end
 
struct Control
    LFeasSpace::Bool
    LFeasPerStage::Bool
    LCostApprox::Bool
    LCostApproxNewCuts::Bool
    LWindStoch::Bool
    LDemandResponse::Bool
    LExtreme::Bool
    LIgnoreCrossCorr::Bool
    MaxIter::Int
    CCMaxIter::Int
    ConvEps::Float64
    NScen::Int
    NWindScen::Int
    NScenSim::Int
    NBranch::Int
    NStage::Int
    NStageSim::Int
    ResInitFrac::Float64
    ResMinFrac::Float64
    MaxResScale::Float64
    LoadScale::Float64
    LineCapScale::Float64
    CapReqFrac::Float64
    H2CompLoss::Float64
end

struct Discrete
   NInfPt::Int 
   NResInitPt::Int
   NResEndPt::Int
   NEnPt::Int
   NRampPt::Int
   NCapPt::Int
end

struct DemandResponse
    NLoadRecStep::Int
    LoadRec::Array{Int}
    MaxUpShift::Array{Float64}
    MaxDnShift::Array{Float64}
    LIncludeExtraConstr::Array{Bool}
    ExtraConstrFilter::Array{Int}
    ExtraConstrSigma::Float64
end


struct Parameters
    Control::Control
    Aggregation::Aggregation
    Time::Time
    Constants::Constants
    Discrete::Discrete
end
 
struct WaterMarkData
    AreaName::String
    NData::Int
    Scaling::Array{Float64}
    AverageInflow::Array{Float64}  
    SeriesNames::Array{String}
end
 
struct ModuleData
    AreaName::String
    UInfName::String
    RInfName::String
    PlantName::String
    ResName::String

    UInfVol::Float64
    RInfVol::Float64
    RegDeg::Float64
    ResBottom::Float64
    MaxBypass::Float64
    MaxFlow::Float64
    NomElevation::Float64
    SubmersionLevel::Float64
    MaxRes::Float64
    ProdCap::Float64
    ConvFac::Float64

    ModId::Int
    ModCnt::Int
    InfFirstYr::Int
    CoplNumber::Int
    InfLastYr::Int
    SpillTo::Int
    BypassTo::Int
    DischargeTo::Int
end
 
struct DynModData
    MaMax::Array{Float64}
    MaMin::Array{Float64}
    QMax::Array{Float64}
    QMin::Array{Float64}
    QfoMax::Array{Float64}
    QfoMin::Array{Float64}
end
    
 
struct UpstreamMod
    NDis::Int
    NByp::Int
    NSpi::Int
    Dis::Array{Int}
    Byp::Array{Int}
    Spi::Array{Int}
end
 
struct AreaUpstreamMod
    USModA::Array{UpstreamMod}
end
 
struct PQCurveData
    NSeg::Int
    MyId::Int
    M3S::Array{Float64}
    DMax::Array{Float64}
    MW::Array{Float64}
    Eff::Array{Float64}
    BestEff::Float64
    PMax::Float64
    PMin::Float64
    QMax::Float64
    MaxDis::Float64
    FirstWeek::Int
    LastWeek::Int
    iCurve::Int
end
 
struct ResCurveData
    NPkt::Int
    Elevation::Array{Float64}
    Volume::Array{Float64}
end
 
struct PumpData
    PumpCap1::Float64
    PumpCap2::Float64
    PumpCap3::Float64
    ToId::Int
end
 
struct AreaDataHydro
    NMod::Int
    EffSea::Array{Float64}
    RegDegComp::Array{Float64}
    MagShare::Array{Float64}
    RegShare::Array{Float64}
    URegShare::Array{Float64}
    RampFrac::Array{Float64}
    MData::Array{ModuleData}
    PQData::Array{PQCurveData}
    RCData::Array{ResCurveData}
end
 
struct HydroSystemData
    AreaNo::Int
    NMod::Int
    ModNo::Array{Int} #Running from 1:AHData[iArea].NMod
    ModId::Array{Int} #Modeled nummer
    ModCnt::Array{Int}#Running from 1:NModTot
    MinRes::Float64
    MaxRes::Float64
    MinProd::Array{Float64}
    MaxProd::Float64
    AveInflow::Float64
end
 
struct HydroSystemArea
    NSys::Int
    Systems::Array{HydroSystemData}
end
 
struct AgrHydro
    MinRes::Float64
    MaxRes::Float64
    MinProd::Array{Float64}
    MaxProd::Float64
end

struct H2Area
    AreaNo::Int
    MaxDis::Float64
    MaxRes::Float64
    CompLoss::Float64
    LStrategic::Bool
end

struct H2System
    NArea::Int
    Ind::Array{Int}
    Areas::Array{H2Area}
end
 
struct MarketStepData
    Name::String
    Capacity::Array{Float64}
    Price::Array{Float64}
end
 
struct MarketLoadData
    Name::String
    Load::Array{Float64}
end
 
struct AreaDataMarket
    NMStep::Int
    NLoad::Int
    MSData::Array{MarketStepData}
    MLData::Array{MarketLoadData}
end
 
struct Connections
    NCon::Int
    LIndxOut::Array{Int}
    LIndxIn::Array{Int}
end
 
struct MarketConnections
    ConArea::Array{Connections}
end

struct EndValuation
    EVType::String
    NEndCut::Int
    EndCutCoef::Array{Float64}
    EndCutRHS::Array{Float64}
    NEndSeg::Int
    EndSegCoef::Array{Float64}
end

struct InflowModel
    NSer::Int
    NResid::Int
    InflowAv::Array{Float64}
    InflowMean::Array{Float64}
    InflowSDev::Array{Float64}
    Resid::Array{Float64}
    CorrMat::Array{Float64}
end

struct SampleScen
    SScen::Array{Int}
    Zscen::Array{Float64}
end

struct Model
    # Read HDF5
    AHData::Array{AreaDataHydro} # Aread hydro data - Array{AreaDataHydro}
    AMData::Array{AreaDataMarket} # Aread market data - Array{AreaDataMarket}
    DMData::DynModData #
    H2Data::H2System #
    USMod::Array{AreaUpstreamMod} # Array{AreaUpstreamMod}
    WPData # Wind power data
    DRData::DemandResponse #Demand response data
    NArea # Number of areas
    ModInfReg # Module regulated inflow - Array{Float64, n_mod, n_week, n_inflow_years}
    ModInfUReg # Module unregulated inflow
    AreaName # Dict{string,int}
    MyKeys # Array{string}
    MaxModArea # Max number of modules in any area
    # Read maske
    MCon::Array{Connections}
    LineCap # Array with line capacities 
    LineLoss # Array with line losses 
    NLine # number of lines 
    # Identify cascades
    NCascade # Total number of hydro modules
    HCascade # Array{HydroSystemArea}
    # Aggregated systems
    NHSys # Number of aggregated hydro systems
    HSys # Array of hydro systems
    NAreaSys # Number of areas in aggregated system
    AreaSys
    USModSys # Area upstream module
    # Inflow
    InfReg # Regulated inflow
    RegFrac # Average fraction of regulated inflow
    # End value
    EV # type EndValuation
end

struct InitialValues
    ResInit::Array{Float64}
    H2Init::Array{Float64}
end

mutable struct Strategy
    StateTraj::Array{Float64, 3} #Hydro storage at the end of the period
    H2Traj::Array{Float64, 3} #H2 storage at the end of the period
    WindYears::Array{Int64} #Wind years used when finding strategy 
    CCR::Array{Float64, 3} #Cut Coefficients for Reservoir
    CCH::Array{Float64, 3} #Cut Coefficients for H2 storage
    CCI::Array{Float64, 3} #Cut Coefficients for Inflow 
    CRHS::Array{Float64, 2} #Cut RHS
    CFP::Array{Float64, 3} #Cost Function Coefficient Production
    CFR::Array{Float64, 3} #Cost Function Coefficient Ramping
    CFC::Array{Float64, 3} #Cost Function Coefficient Reserve
    CFHD::Array{Float64, 3} #Cost Function Coefficient H2 discharge
    CFHS::Array{Float64, 3} #Cost Function Coefficient H2 storage
    CFRHS::Array{Float64, 2} #Cost Function RHS
    NCut::Int # Numbers of cuts
    NCostCut::Int # Numbers of cost cuts
    MaxIter::Int # Itertion number to continue from
    NWindScen::Int # Number of VRES scenarios
    LWindStoch::Bool # If VRES uncertainty in strategy computation
end

struct FeasibilitySpace
    iStage::Int #stage created for
    FCC::Array{Float64, 3} #Cut coefficients 
    NFeasCut::Vector{Int} #Numbers of feasibility cuts 
    NFeasCutMax::Int #Max dimension 
end

struct Result
    ReservoirTable::Array{Float64}
    HProdTable::Array{Float64}
    HRampTable::Array{Float64}
    HCapTable::Array{Float64}
    MarkTable::Array{Float64}
    FlowTable::Array{Float64}
    SpillTable::Array{Float64}
    InflowTable::Array{Float64}
    LoadTable::Array{Float64}
    WindTable::Array{Float64}
    PriceTable::Array{Float64}
    RationingTable::Array{Float64}
    DemandUpTable::Array{Float64}
    DemandDnTable::Array{Float64}
    H2StoreTable::Array{Float64}
    H2DisTable::Array{Float64}
end

struct DetailedResult
    ReservoirTable::Array{Float64}
    HProdTable::Array{Float64}
    MarkTable::Array{Float64}
    FlowTable::Array{Float64}
    DischargeTable::Array{Float64}
    SpillTable::Array{Float64}
    BypassTable::Array{Float64}
    LoadTable::Array{Float64}
    WindTable::Array{Float64}
    PriceTable::Array{Float64}
    RationingTable::Array{Float64}
    DemandUpTable::Array{Float64}
    DemandDnTable::Array{Float64}
    H2StoreTable::Array{Float64}
    H2DisTable::Array{Float64}
end
