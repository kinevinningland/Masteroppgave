function load(dataset::String, parameters::Parameters)::Model
    CNS = parameters.Constants
    CTI = parameters.Time
    CTR = parameters.Control
    CAGR = parameters.Aggregation

    AHData,AMData,USMod,WPData,NArea,ModInfReg,ModInfUReg,AreaName,MyKeys,MaxModArea = ReadHDF5(dataset,CNS,CTI,CTR)
    MCon,LineCap,LineLoss,NLine = ReadMaske(dataset,NArea,CNS,CTR)
    DMData = ReadDynmod(dataset,AHData,NArea,AreaName,MyKeys,MaxModArea,CTI.NWeek)
    DRData = ReadDemandResponse(dataset,NArea,CTI.NWeek,AreaName,CTR.LDemandResponse)
    H2Data = ReadH2(dataset,AHData,NArea,AreaName,CNS,CTI,CTR)
    PrintMarketSummary(dataset,NArea,AMData,WPData,MCon,LineCap,CTI)
    

    NCascade,HCascade = IdentifyCascades(AHData,NArea,USMod,ModInfReg,ModInfUReg,CTR,CTI,CNS)
    NHSys,HSys,NAreaSys,AreaSys,USModSys = AggrSystems(HCascade,AHData,NArea,ModInfReg,ModInfUReg,CTR,CTI,CNS,CAGR)
    InfReg,RegFrac = AggrInflow(NHSys,HSys,ModInfReg,ModInfUReg,AHData,CTR,CTI,CNS)

    ORData = ReadOperatingReserves(NArea, NHSys, NAreaSys,AreaSys, H2Data, AMData,CTR.LOperatingReserves) #Added

    PrintCascadeStats(dataset,NHSys,HSys,CAGR)
    PrintInflow(dataset,InfReg,NHSys,CTI.NWeek,CTR.NBranch,CTI.NInflowYear)

    EV = ReadEndValue(dataset,NHSys)

    model = Model(AHData, AMData, DMData, H2Data, USMod, WPData, DRData, ORData,NArea, ModInfReg,
        ModInfUReg, AreaName, MyKeys, MaxModArea, MCon, LineCap, LineLoss, NLine,
        NCascade, HCascade, NHSys,HSys, NAreaSys, AreaSys, USModSys,
        InfReg, RegFrac, EV) #ORData Added

    return model
end

function load_inflow(dataset::String, model::Model, parameters::Parameters)::InflowModel

    inflow_model = ReadInflow(dataset, model.NHSys, parameters.Time.NWeek, parameters.Control.NScen, 
                              parameters.Control.NBranch, parameters.Control.LIgnoreCrossCorr)

    return inflow_model
end

function load_feas_cut(dataset::String, model::Model, parameters::Parameters)::FeasibilityCuts
    NFeasCut,FCC = ReadFeasCut(parameters.Control.LFeasSpace,model.NHSys,dataset)
    feas_cut = FeasibilityCuts(NFeasCut, FCC)
    return feas_cut
end
