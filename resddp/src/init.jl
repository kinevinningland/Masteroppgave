function init_strategy(model::Model, parameters::Parameters)::Strategy
    NHSys = model.NHSys
    NScen = parameters.Control.NScen
    NWindScen = parameters.Control.NWindScen
    NStage = parameters.Control.NStage
    MaxIter = parameters.Control.MaxIter
    CCMaxIter = parameters.Control.CCMaxIter
    NH2Area = model.H2Data.NArea

    WindStoch = true
    if !parameters.Control.LWindStoch
        NWindScen = 1
        WindStoch = false
    end

    StateTraj = zeros(NHSys,NScen,NStage)   #Hydro storage at the end of the period
    H2Traj = zeros(NH2Area,NScen,NStage)    #H2 storage at the end of the period
    WindYears = zeros(Int64,NWindScen) #Wind years used when creating strategy

    CCR = zeros(NHSys,NStage,NScen*MaxIter) #Cut Coefficients for Reservoir
    CCH = zeros(NH2Area,NStage,NScen*MaxIter) #Cut Coefficients for H2 storage
    CCI = zeros(NHSys,NStage,NScen*MaxIter) #Cut Coefficients for Inflow 
    CRHS = zeros(NStage,NScen*MaxIter)      #Cut RHS
    CFP = zeros(NHSys,NStage,NScen*MaxIter*CCMaxIter)    #Cost Function Coefficient Production
    CFR = zeros(NHSys,NStage,NScen*MaxIter*CCMaxIter)    #Cost Function Coefficient Ramping
    CFC = zeros(NHSys,NStage,NScen*MaxIter*CCMaxIter)    #Cost Function Coefficient Reserve
    CFHD = zeros(NH2Area,NStage,NScen*MaxIter*CCMaxIter)  #Cost Function Coefficient H2 Discharge
    CFHS = zeros(NH2Area,NStage,NScen*MaxIter*CCMaxIter)  #Cost Function Coefficient H2 Storage
    CFRHS = zeros(NStage,NScen*MaxIter*CCMaxIter)        #Cost Function RHS

    return Strategy(StateTraj, H2Traj, WindYears, CCR, CCH, CCI, CRHS, CFP, CFR, CFC, CFHD, CFHS, CFRHS, 0, 0, 1, NWindScen, WindStoch)
end

function extend_strategy(model::Model, parameters::Parameters, oldstrategy::Strategy)::Strategy
    NHSys = model.NHSys
    NH2Area = model.H2Data.NArea
    NScen = parameters.Control.NScen
    NStage = parameters.Control.NStage
    MaxIter = parameters.Control.MaxIter
    CCMaxIter = parameters.Control.CCMaxIter
    WindStoch = false
    NWindScen = length(oldstrategy.WindYears)
    if NWindScen > 1
        WindStoch = true
    end

    NOldCut = oldstrategy.NCut
    NOldCostCut = oldstrategy.NCostCut
    OldMaxIter = oldstrategy.MaxIter

    StateTraj = zeros(NHSys,NScen,NStage)   #Hydro storage at the end of the period
    H2Traj = zeros(NH2Area,NScen,NStage)    #H2 storage at the end of the period
    WindYears = zeros(NWindScen)

    CCR = zeros(NHSys,NStage,NOldCut+NScen*MaxIter) #Cut Coefficients for Reservoir
    CCH = zeros(NH2Area,NStage,NOldCut+NScen*MaxIter) #Cut Coefficients for H2-Storage
    CCI = zeros(NHSys,NStage,NOldCut+NScen*MaxIter) #Cut Coefficients for Inflow 
    CRHS = zeros(NStage,NOldCut+NScen*MaxIter)      #Cut RHS

    if parameters.Control.LCostApproxNewCuts
        CFP = zeros(NHSys,NStage,NOldCostCut+NScen*MaxIter*CCMaxIter)    #Cost Function Coefficient Production
        CFR = zeros(NHSys,NStage,NOldCostCut+NScen*MaxIter*CCMaxIter)    #Cost Function Coefficient Ramping
        CFC = zeros(NHSys,NStage,NOldCostCut+NScen*MaxIter*CCMaxIter)    #Cost Function Coefficient Reserve
        CFHD = zeros(NH2Area,NStage,NOldCostCut+NScen*MaxIter*CCMaxIter) #Cost Function Coefficient H2-Discharge
        CFHS = zeros(NH2Area,NStage,NOldCostCut+NScen*MaxIter*CCMaxIter) #Cost Function Coefficient H2-Storage
        CFRHS = zeros(NStage,NOldCostCut+NScen*MaxIter*CCMaxIter)        #Cost Function RHS
    else
        CFP = zeros(NHSys,NStage,NOldCostCut)    #Cost Function Coefficient Production
        CFR = zeros(NHSys,NStage,NOldCostCut)    #Cost Function Coefficient Ramping
        CFC = zeros(NHSys,NStage,NOldCostCut)    #Cost Function Coefficient Reserve
        CFHD = zeros(NH2Area,NStage,NOldCostCut) #Cost Function Coefficient H2-Discharge
        CFHS = zeros(NH2Area,NStage,NOldCostCut) #Cost Function Coefficient H2-Storage
        CFRHS = zeros(NStage,NOldCostCut)        #Cost Function RHS
    end

    #Copy values
    CCR[1:NHSys,1:NStage,1:NOldCut] = oldstrategy.CCR[1:NHSys,1:NStage,1:NOldCut] 
    CCH[1:NH2Area,1:NStage,1:NOldCut] = oldstrategy.CCH[1:NH2Area,1:NStage,1:NOldCut] 
    CCI[1:NHSys,1:NStage,1:NOldCut] = oldstrategy.CCI[1:NHSys,1:NStage,1:NOldCut] 
    CRHS[1:NStage,1:NOldCut] = oldstrategy.CRHS[1:NStage,1:NOldCut] 
    CFP[1:NHSys,1:NStage,1:NOldCostCut] = oldstrategy.CFP[1:NHSys,1:NStage,1:NOldCostCut] 
    CFR[1:NHSys,1:NStage,1:NOldCostCut] = oldstrategy.CFR[1:NHSys,1:NStage,1:NOldCostCut] 
    CFC[1:NHSys,1:NStage,1:NOldCostCut] = oldstrategy.CFC[1:NHSys,1:NStage,1:NOldCostCut] 
    CFHD[1:NH2Area,1:NStage,1:NOldCostCut] = oldstrategy.CFHD[1:NH2Area,1:NStage,1:NOldCostCut] 
    CFHS[1:NH2Area,1:NStage,1:NOldCostCut] = oldstrategy.CFHS[1:NH2Area,1:NStage,1:NOldCostCut] 
    CFRHS[1:NStage,1:NOldCostCut] = oldstrategy.CFRHS[1:NStage,1:NOldCostCut] 
    WindYears[1:NWindScen] = oldstrategy.WindYears[1:NWindScen]

    println("..reading ",NOldCut," old cuts and ",NOldCostCut," cost cuts")
    return Strategy(StateTraj, H2Traj, WindYears, CCR, CCH, CCI, CRHS, CFP, CFR, CFC, CFHD, CFHS, CFRHS, NOldCut, NOldCostCut, OldMaxIter, NWindScen, WindStoch)
end

function init_system(model::Model, parameters::Parameters)::InitialValues
    ResInit0 = [parameters.Control.ResInitFrac * model.HSys[iSys].MaxRes for iSys in 1:model.NHSys]
    H2Init0 = [parameters.Control.ResInitFrac * model.H2Data.Areas[iSys].MaxRes for iSys in 1:model.H2Data.NArea]
    return InitialValues(ResInit0,H2Init0)
end

function init_result(NArea,NHSys,NMaxMStep,NScen,NStage,NK,NLine)::Result
    ReservoirTable = zeros(Float64,NHSys,NScen,NStage,NK) 
    HProdTable = zeros(Float64,NHSys,NScen,NStage,NK) 
    HRampTable = zeros(Float64,NHSys,NScen,NStage)  
    HCapTable = zeros(Float64,NHSys,NScen,NStage)   
    MarkTable = zeros(Float64,NArea,NMaxMStep,NScen,NStage,NK) 
    FlowTable = zeros(Float64,NLine,NScen,NStage,NK) 
    SpillTable = zeros(Float64,NHSys,NScen,NStage,NK) 
    InflowTable = zeros(Float64,NHSys,NScen,NStage) 
    LoadTable = zeros(Float64,NArea,NScen,NStage,NK) 
    WindTable = zeros(Float64,NArea,NScen,NStage,NK) 
    PriceTable = zeros(Float64,NArea,NScen,NStage,NK)
    RationingTable = zeros(Float64,NArea,NScen,NStage,NK)
    DemandUpTable = zeros(Float64,NArea,NScen,NStage,NK) 
    DemandDnTable = zeros(Float64,NArea,NScen,NStage,NK)
    H2StoreTable = zeros(Float64,NArea,NScen,NStage,NK)
    H2DisTable = zeros(Float64,NArea,NScen,NStage,NK)
    
    return Result(ReservoirTable,HProdTable,HRampTable,HCapTable,MarkTable,FlowTable,SpillTable,InflowTable,LoadTable,WindTable,PriceTable,
                  RationingTable,DemandUpTable,DemandDnTable,H2StoreTable,H2DisTable)
 end

 function init_detailed_result(NArea,NHSys,NMaxMStep,NScen,NStage,NK,NLine,NMaxMod)::DetailedResult

    ReservoirTable = zeros(Float64,NHSys,NMaxMod,NScen,NStage,NK)
    HProdTable = zeros(Float64,NHSys,NMaxMod,NScen,NStage,NK) 
    MarkTable = zeros(Float64,NArea,NMaxMStep,NScen,NStage,NK) 
    FlowTable = zeros(Float64,NLine,NScen,NStage,NK)
    DischargeTable = zeros(Float64,NHSys,NMaxMod,NScen,NStage,NK)
    SpillTable = zeros(Float64,NHSys,NMaxMod,NScen,NStage,NK)
    BypassTable = zeros(Float64,NHSys,NMaxMod,NScen,NStage,NK)
    LoadTable = zeros(Float64,NArea,NScen,NStage,NK) 
    WindTable = zeros(Float64,NArea,NScen,NStage,NK) 
    PriceTable = zeros(Float64,NArea,NScen,NStage,NK)
    RationingTable = zeros(Float64,NArea,NScen,NStage,NK)
    DemandUpTable = zeros(Float64,NArea,NScen,NStage,NK) 
    DemandDnTable = zeros(Float64,NArea,NScen,NStage,NK)
    H2StoreTable = zeros(Float64,NArea,NScen,NStage,NK)
    H2DisTable = zeros(Float64,NArea,NScen,NStage,NK)

    return DetailedResult(ReservoirTable,HProdTable,MarkTable,FlowTable,DischargeTable,SpillTable,BypassTable,LoadTable,WindTable,PriceTable,
                          RationingTable,DemandUpTable,DemandDnTable,H2StoreTable,H2DisTable)
 end

