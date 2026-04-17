function simulate_detailed(model::Model, inflow_model::InflowModel, parameters::Parameters, strategy::Strategy; optimizer=JuMP.optimizer_with_attributes(Clp.Optimizer, "SolveType" => 0, "PresolveType" => 1, "LogLevel" => 0))::DetailedResult

    NMaxMod = maximum([model.AHData[iSys].NMod for iSys in 1:model.NHSys])

    SimulatedStateTraj = zeros(Float64,model.NHSys, NMaxMod, parameters.Control.NScenSim, parameters.Control.NStageSim)
    SimulatedCost = zeros(Float64, parameters.Control.NScenSim)

    NMaxMStep = maximum([model.AMData[iArea].NMStep for iArea in 1:model.NArea])

    DetailedResultTable = init_detailed_result(model.NArea, model.NHSys, NMaxMStep, parameters.Control.NScenSim, parameters.Control.NStageSim, parameters.Time.NK, model.NLine, NMaxMod)

    ResInit0 = zeros(Float64,model.NHSys,NMaxMod)

    for iSys = 1:model.NHSys
       for iMod = 1:model.AHData[iSys].NMod
          ResInit0[iSys,iMod] = parameters.Control.ResInitFrac*parameters.Control.MaxResScale*model.AHData[iSys].MData[iMod].MaxRes
       end
    end

    NCluster = min(Threads.nthreads(), parameters.Control.NScenSim) # Never more threads than scenarios in simulation
    NScenPerCluster = Int(ceil(parameters.Control.NScenSim/NCluster)) # Maximum number of scenario per thread
    
    dTS1 = dTS2 = 0.0
    t1 = time_ns()
    for t = 1:parameters.Control.NStageSim
        println("Stage ",t)
        sWeek = mod1(t,parameters.Time.NWeek)

        Threads.@threads for sCluster = 1:NCluster
            start_scen = (sCluster-1) * NScenPerCluster + 1
            if start_scen <= parameters.Control.NScenSim

                end_scen = min(sCluster * NScenPerCluster, parameters.Control.NScenSim)
                ResInit = zeros(Float64,model.NHSys,NMaxMod)

                ResInit = zeros(Float64,model.NHSys,NMaxMod)

                SP_FORW = StageProbDet.Build(t,sWeek,model.USMod,model.AHData,model.AMData,model.HSys,model.MCon,
                model.EV,strategy.CCR,strategy.CCH,parameters.Constants,strategy.NCut,model.NHSys,model.NArea,
                model.NLine,model.LineCap,model.LineLoss,parameters.Time,t==parameters.Control.NStageSim,
                parameters.Control.LDemandResponse,model.DRData,optimizer)

                for iScen = start_scen:end_scen
                    if t > 1
                        ResInit[1:model.NHSys,1:NMaxMod] = SimulatedStateTraj[1:model.NHSys,1:NMaxMod,iScen,t-1] 
                    else
                        ResInit[1:model.NHSys,1:NMaxMod] = ResInit0[1:model.NHSys,1:NMaxMod]
                    end

                    for iArea=1:model.NHSys
                        for iMod=1:model.AHData[iArea].NMod
                            CurrInf = parameters.Time.WeekFrac*(model.ModInfReg[model.AHData[iArea].MData[iMod].ModCnt,sWeek,iScen] + model.ModInfUReg[model.AHData[iArea].MData[iMod].ModCnt,sWeek,iScen])

                            JuMP.set_normalized_rhs(SP_FORW[:resbalReg0][iArea,iMod], ResInit[iArea,iMod] + CurrInf)
                            for k=2:parameters.Time.NK
                                JuMP.set_normalized_rhs(SP_FORW[:resbalReg][iArea,iMod,k], CurrInf)
                            end
                        end
                    end
                    for iArea=1:model.NArea
                        for k=1:parameters.Time.NK
                            JuMP.set_normalized_rhs(SP_FORW[:wptarget][iArea,k], max(model.WPData[iArea,iScen,sWeek,k],0.0))
                        end
                    end
                
                    if t < parameters.Control.NStageSim
                        Ztilst = zeros(Float64,model.NHSys)
                        for iSys = 1:model.NHSys
                            qSys = 0.0
                            for iMod = 1:model.AHData[iSys].NMod
                                myNr = model.AHData[iSys].MData[iMod].ModCnt
                                qSys += (model.ModInfReg[myNr,sWeek,iScen]+ model.ModInfUReg[myNr,sWeek,iScen])*model.AHData[iSys].EffSea[iMod]*parameters.Constants.MAGEFF2GWH
                            end
                            Ztilst[iSys] = (qSys-inflow_model.InflowMean[iSys,sWeek])/inflow_model.InflowSDev[iSys,sWeek]
                        end
                        adjust = transpose(strategy.CCI[1:inflow_model.NSer,t,1:strategy.NCut])*Ztilst[1:inflow_model.NSer]
                        for iCut = 1:strategy.NCut
                            JuMP.set_normalized_rhs(SP_FORW[:cut][iCut],strategy.CRHS[t,iCut]+adjust[iCut])
                        end
                    end

                    #Solve problem
                    t2 = time_ns()
                    optimize!(SP_FORW)
                    dTS2 = dTS2+(time_ns()-t2)/NCluster
                    if primal_status(SP_FORW) != MOI.FEASIBLE_POINT 
                        write_to_file(SP_FORW, "SPF_err.lp")
                        termstat = termination_status(SP_FORW)
                        error(println("Solver terminated with status $termstat in forward iteration (stage,scen): ",t," ",iScen))
                    end
                    save_detailed!(DetailedResultTable, SP_FORW, model.AMData, model.AHData, model.NArea, model.NHSys, parameters.Time.NK, model.NLine, iScen, t)

                    for iSys = 1:model.NHSys
                        for iMod = 1:model.AHData[iSys].NMod
                            SimulatedStateTraj[iSys,iMod,iScen,t] = JuMP.value(SP_FORW[:res][iSys,iMod,parameters.Time.NK])
                        end
                    end

                    if t < parameters.Control.NStageSim
                        SimulatedCost[iScen] += (JuMP.objective_value(SP_FORW)-JuMP.value(SP_FORW[:alpha]))
                    else
                        SimulatedCost[iScen] += JuMP.objective_value(SP_FORW)
                    end
                end
            end
        end
    end
    @printf("%s %6.2f \n"," Time - Total : ",(time_ns()-t1)*1.0E-9)
    @printf("%s %6.2f \n"," Time - Solver: ",dTS2*1.0E-9)

    return DetailedResultTable
end

function simulate_aggregated(model::Model, inflow_model::InflowModel, parameters::Parameters, strategy::Strategy, feas_spaces::Vector{FeasibilitySpace}, 
        initial_values::InitialValues; optimizer=JuMP.optimizer_with_attributes(Clp.Optimizer, "SolveType" => 0, "PresolveType" => 1, "LogLevel" => 0), fixed_seed = false)::Result

    SimulatedStateTraj = zeros(Float64,model.NHSys, parameters.Control.NScenSim, parameters.Control.NStageSim)
    SimulatedH2Traj = zeros(Float64,model.H2Data.NArea, parameters.Control.NScenSim, parameters.Control.NStageSim)
    SimulatedCost = zeros(Float64, parameters.Control.NScenSim)

    SS = SampleScenario(parameters.Control.NScenSim, parameters.Control.NStageSim, parameters.Time.NWeek, inflow_model, parameters.Control.LExtreme; fixed_seed = fixed_seed)

    NMaxMStep = maximum([model.AMData[iArea].NMStep for iArea in 1:model.NArea])

    ResultTable = init_result(model.NArea, model.NHSys, NMaxMStep, parameters.Control.NScenSim, parameters.Control.NStageSim, parameters.Time.NK, model.NLine)

    #START SIMULATION
    H2Init = zeros(Float64,model.H2Data.NArea)

    NCluster = min(Threads.nthreads(), parameters.Control.NScenSim) # Never more threads than scenarios in simulation
    NScenPerCluster = Int(ceil(parameters.Control.NScenSim/NCluster)) # Maximum number of scenario per thread

    dTS1 = dTS2 = 0.0
    t1 = time_ns()
    println("VRES-years simulated: ",strategy.WindYears)
    for t = 1:parameters.Control.NStageSim
        println("Stage ",t)
        sWeek = mod1(t, parameters.Time.NWeek)
        fWeek = sWeek; if !parameters.Control.LFeasPerStage; fWeek = 1; end

        Threads.@threads for sCluster = 1:NCluster
            start_scen = (sCluster-1) * NScenPerCluster + 1
            if start_scen <= parameters.Control.NScenSim

                end_scen = min(sCluster * NScenPerCluster, parameters.Control.NScenSim)
                ResInit = zeros(Float64, model.NHSys)

                ResInit = zeros(Float64, model.NHSys)

                wYear = strategy.WindYears[1]
                MyWPData = model.WPData[1:model.NArea, wYear, sWeek, 1:parameters.Time.NK]
                SP_FORW = StageProbFull.Build(t, sWeek, inflow_model, model.NHSys,
                    model.HSys, model.NAreaSys, model.AreaSys, model.AMData,
                    model.MCon, model.EV, strategy.CCR, strategy.CCH, strategy.CCI,
                    parameters.Constants, strategy.NCut, model.NArea, model.NLine,
                    model.LineCap,model.LineLoss, parameters.Time, t==parameters.Control.NStage,
                    parameters.Control.LFeasSpace, feas_spaces[fWeek].NFeasCut, feas_spaces[fWeek].FCC,
                    parameters.Control.CapReqFrac,MyWPData,parameters.Control.LDemandResponse,model.DRData,
                    model.H2Data,optimizer)

                for iScen = start_scen:end_scen
                    wYear = sample(strategy.WindYears)
                    #Update constraint right-hand sides
                    Zstate = zeros(Float64,inflow_model.NSer)
                    Eps = inflow_model.Resid[1:inflow_model.NSer,sWeek,SS.SScen[iScen,t]]
                    if t > 1
                        ResInit[1:model.NHSys] = SimulatedStateTraj[1:model.NHSys,iScen,t-1,end]
                        H2Init[1:model.H2Data.NArea] = SimulatedH2Traj[1:model.H2Data.NArea,iScen,t-1,end]
                        Zstate = SS.Zscen[iScen,t-1,1:inflow_model.NSer]
                    else
                        ResInit[1:model.NHSys] = initial_values.ResInit[1:model.NHSys]
                        H2Init[1:model.H2Data.NArea] = initial_values.H2Init[1:model.H2Data.NArea]
                    end

                    DevEpsVec = inflow_model.InflowSDev[1:inflow_model.NSer,sWeek].*Eps             #Dim(NSer*1)
                    CorrStateVec = inflow_model.CorrMat*Zstate                                      #Dim(NSer*1)
                    DevStateVec = inflow_model.InflowSDev[1:inflow_model.NSer,sWeek].*CorrStateVec  #Dim(NSer*1)
                    InflowSys = DevStateVec+DevEpsVec+inflow_model.InflowMean[1:inflow_model.NSer,sWeek]

                    for iSys = 1:model.NHSys
                        JuMP.set_normalized_rhs(SP_FORW[:rstate][iSys],ResInit[iSys])
                        JuMP.set_normalized_rhs(SP_FORW[:zstate][iSys],Zstate[iSys])
                        JuMP.set_normalized_rhs(SP_FORW[:inflow][iSys],DevEpsVec[iSys]+inflow_model.InflowMean[iSys,sWeek])
                    end
                    for iH2a = 1:model.H2Data.NArea
                        JuMP.set_normalized_rhs(SP_FORW[:h2storage0][iH2a,1],H2Init[iH2a])
                    end
                    if t < parameters.Control.NStageSim
                        adjust = transpose(strategy.CCI[1:inflow_model.NSer, t, 1:strategy.NCut]) * Eps[1:inflow_model.NSer]
                        for iCut = 1:strategy.NCut
                            JuMP.set_normalized_rhs(SP_FORW[:cut][iCut], strategy.CRHS[t,iCut]+adjust[iCut])
                        end
                    end
                    for iArea = 1:model.NArea
                        for k = 1:parameters.Time.NK
                            JuMP.set_normalized_rhs(SP_FORW[:wptarget][iArea,k], max(model.WPData[iArea,wYear,sWeek,k],0.0))
                        end
                    end

                    #Solve problem
                    t2 = time_ns()
                    optimize!(SP_FORW)
                    #write_to_file(SP_FORW, "SPF.lp")
                    #exit()
                    dTS2 = dTS2+(time_ns()-t2)/NCluster
                    if primal_status(SP_FORW) != MOI.FEASIBLE_POINT ||  dual_status(SP_FORW) != MOI.FEASIBLE_POINT
                        write_to_file(SP_FORW, "SPF_err.lp")
                        write_to_file(SP_FORW, "SPF_err.mps")
                        termstat = termination_status(SP_FORW)
                        error(println("Solver terminated with status $termstat in forward iteration (stage,scen): ",t," ",iScen))
                    end
                    save!(ResultTable, SP_FORW, model.AMData,model.H2Data, InflowSys, model.NArea, model.NHSys, parameters.Time.NK, model.NLine, iScen, t)

                    for iSys = 1:model.NHSys
                        SimulatedStateTraj[iSys,iScen,t] = JuMP.value(SP_FORW[:res][iSys,end])
                    end
                    for iH2a = 1:model.H2Data.NArea
                        SimulatedH2Traj[iH2a,iScen,t] = JuMP.value(SP_FORW[:h2res][iH2a,end])
                    end

                    if t < parameters.Control.NStageSim
                        SimulatedCost[iScen] += (JuMP.objective_value(SP_FORW)-JuMP.value(SP_FORW[:alpha]))
                    else
                        SimulatedCost[iScen] += JuMP.objective_value(SP_FORW)
                    end
                end
            end
        end
    end
    @printf("%s %6.2f \n"," Time - Total : ",(time_ns()-t1)*1.0E-9)
    @printf("%s %6.2f \n"," Time - Solver: ",dTS2*1.0E-9)

    return ResultTable
end
