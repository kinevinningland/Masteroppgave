function SolveFWD(model::Model, inflow_model::InflowModel, feas_spaces::Vector{FeasibilitySpace},
    parameters::Parameters, strategy::Strategy, init_val::InitialValues, SS::SampleScen,  scenario_range::UnitRange{Int}, optimizer)

    NScenPerCluster = length(scenario_range)
    NStage = parameters.Control.NStage
    NHSys = model.NHSys
    NH2Area = model.H2Data.NArea 
    NWeek = parameters.Time.NWeek
    HSys = model.HSys
    
    StartScen = scenario_range[begin]
    EndScen = scenario_range[end]
 
    NDim = (parameters.Control.CCMaxIter*(3*NHSys+2*NH2Area+1)+NHSys+NH2Area)*NScenPerCluster*NStage+NScenPerCluster+1+1
    ResArray = zeros(Float64,NDim)
    SimulatedStateTraj = zeros(Float64,NHSys,NScenPerCluster,NStage)
    SimulatedH2Traj = zeros(Float64,NHSys,NScenPerCluster,NStage)
    SimulatedCost = zeros(Float64,NScenPerCluster)
    SumBackward = 0.0
 
    dTF = 0.0
    ResInit = zeros(Float64,NHSys)
    H2Init = zeros(Float64,NH2Area)
 
 
    jCnt = 1
    for t = 1:NStage
        iWeek = mod1(t,NWeek)
        fWeek = iWeek; if !parameters.Control.LFeasPerStage; fWeek = 1; end

        wYear = strategy.WindYears[1]
        MyWPData = model.WPData[1:model.NArea, wYear, iWeek, 1:parameters.Time.NK]
        if parameters.Control.LCostApprox
            SP_FORW = StageProb.Build(t, iWeek, inflow_model, NHSys, HSys,
                model.AreaSys, model.AMData, model.EV,
                strategy.CCR, strategy.CCH, strategy.CCI, strategy.CFP, strategy.CFR,
                strategy.CFC, strategy.CFHD, strategy.CFHS, strategy.CFRHS, parameters.Constants, strategy.NCut, strategy.NCostCut,
                model.NArea, parameters.Time, t==NStage, parameters.Control.LFeasSpace, 
                feas_spaces[fWeek].NFeasCut, feas_spaces[fWeek].FCC, model.H2Data, optimizer)
            if parameters.Control.LCostApproxNewCuts
                CP = CostProb.Build(t, iWeek, NHSys, HSys, model.NAreaSys,
                    model.AreaSys, model.AMData, model.MCon, parameters.Constants,
                    model.NArea, model.NLine, model.LineCap,model.LineLoss, parameters.Time, MyWPData, 
                    parameters.Control.LDemandResponse, model.DRData, model.H2Data, optimizer)
            end
        else
            SP_FORW = StageProbFull.Build(t, iWeek, inflow_model, NHSys, HSys,
                model.NAreaSys, model.AreaSys, model.AMData,
                model.MCon, model.EV, strategy.CCR,strategy.CCH, strategy.CCI,
                parameters.Constants, strategy.NCut, model.NArea, model.NLine,
                model.LineCap,model.LineLoss, parameters.Time, t==NStage,
                parameters.Control.LFeasSpace, feas_spaces[fWeek].NFeasCut, feas_spaces[fWeek].FCC,
                parameters.Control.CapReqFrac,MyWPData, parameters.Control.LDemandResponse, model.DRData, model.H2Data, optimizer)
        end
    
        iScen = 0
        for aScen = StartScen:EndScen
            iScen += 1
    
            #Update constraints according to state
            Zstate = zeros(Float64,inflow_model.NSer)
            Eps = inflow_model.Resid[1:inflow_model.NSer,iWeek,SS.SScen[aScen,t]]
            if t > 1
                ResInit[1:NHSys] = SimulatedStateTraj[1:NHSys,iScen,t-1,end]
                H2Init[1:NH2Area] = SimulatedH2Traj[1:NH2Area,iScen,t-1,end]
                Zstate = SS.Zscen[aScen,t-1,1:inflow_model.NSer]
            else
                ResInit[1:NHSys] = init_val.ResInit[1:NHSys]
                H2Init[1:NH2Area] = init_val.H2Init[1:NH2Area]
            end
            
            DevEpsVec = inflow_model.InflowSDev[1:inflow_model.NSer,iWeek].*Eps          #Dim(NSer*1)
            #CorrStateVec = inflow_model.CorrMat*Zstate                            #Dim(NSer*1)
            #DevStateVec = inflow_model.InflowSDev[1:inflow_model.NSer,iWeek].*CorrStateVec  #Dim(NSer*1)
            #InflowSys = DevStateVec+DevEpsVec+inflow_model.InflowMean[1:inflow_model.NSer,iWeek]
    
            for iSys = 1:NHSys
                JuMP.set_normalized_rhs(SP_FORW[:rstate][iSys],ResInit[iSys])
                JuMP.set_normalized_rhs(SP_FORW[:zstate][iSys],Zstate[iSys])
                JuMP.set_normalized_rhs(SP_FORW[:inflow][iSys],DevEpsVec[iSys]+inflow_model.InflowMean[iSys,iWeek])
            end
            for iArea = 1:NH2Area
                JuMP.set_normalized_rhs(SP_FORW[:h2state][iArea],H2Init[iArea])
            end
            
            if t < NStage
                adjust = transpose(strategy.CCI[1:inflow_model.NSer,t,1:strategy.NCut])*Eps[1:inflow_model.NSer]
                for iCut = 1:strategy.NCut
                    JuMP.set_normalized_rhs(SP_FORW[:cut][iCut],strategy.CRHS[t,iCut]+adjust[iCut])
                end
            end
    
            if parameters.Control.LCostApprox
                HydroProd = zeros(Float64,NHSys)
                RampRate = zeros(Float64,NHSys)
                CapProc = zeros(Float64,NHSys)
                H2Discharge = zeros(Float64,NH2Area)
                #myCostCons = Vector{ConstraintRef}(undef,CTR.CCMaxIter)

                
                t1 = time_ns()

                if parameters.Control.LCostApproxNewCuts
                    for c = 1:parameters.Control.CCMaxIter
                        MyCFP = zeros(Float64,NHSys)
                        MyCFR = zeros(Float64,NHSys)
                        MyCFC = zeros(Float64,NHSys)
                        MyCFHD = zeros(Float64,NH2Area)
                        MyCFHS = zeros(Float64,NH2Area)
                        MyCFRHS = 0.0
                        ObjCost = 0.0

                        optimize!(SP_FORW)
                        if primal_status(SP_FORW) != MOI.FEASIBLE_POINT 
                            write_to_file(SP_FORW, "SP_FORW_err.lp")
                            termstat = termination_status(SP_FORW)
                            error(println("Solver terminated with status $termstat for stageprob in forward iteration (stage,scen): ",t," ",aScen))
                        end
                        ObjRel = JuMP.objective_value(SP_FORW)
                        Beta = JuMP.value(SP_FORW[:beta])
                        for iSys = 1:NHSys
                           HydroProd[iSys] = max(0.0,JuMP.value(SP_FORW[:prod][iSys]))
                           #RampRate[iSys] = max(0.0,JuMP.value(SP_FORW[:ramp][iSys]))
                           #CapProc[iSys] = max(0.0,JuMP.value(SP_FORW[:cap][iSys]))
                           JuMP.set_normalized_rhs(CP[:hydprod][iSys],HydroProd[iSys])
                           #JuMP.set_normalized_rhs(CP[:rampcap][iSys],RampRate[iSys])
                           #JuMP.set_normalized_rhs(CP[:capproc][iSys],CapProc[iSys])
                        end
                        for iArea = 1:NH2Area
                           H2Discharge[iArea] = JuMP.value(SP_FORW[:h2dis][iArea])
                           JuMP.set_normalized_rhs(CP[:h2discharge][iArea],H2Discharge[iArea])
                           JuMP.set_normalized_rhs(CP[:h2storage0][iArea,1],H2Init[iArea])
                        end
        
                        #Wind power uncertainty (fan of scenarios) 
                        if parameters.Control.LWindStoch
                            NWind = strategy.NWindScen
                            for wScen = 1:NWind
                                wYear = strategy.WindYears[wScen]
                                for iArea = 1:model.NArea
                                    for k = 1:parameters.Time.NK
                                        JuMP.set_normalized_rhs(CP[:wptarget][iArea,k],max(model.WPData[iArea,wYear,iWeek,k],0.0))
                                    end
                                end
            
                                optimize!(CP)
                                if primal_status(CP) != MOI.FEASIBLE_POINT 
                                    write_to_file(CP, "CP_err.lp")
                                    termstat = termination_status(CP)
                                    error(println("Solver terminated with status $termstat for costprob in forward iteration (stage,scen): ",t," ",aScen))
                                end
            
                                ObjCost += JuMP.objective_value(CP)
                                for iSys = 1:NHSys
                                    MyCFP[iSys] += -JuMP.shadow_price(CP[:hydprod][iSys])
                                    #MyCFR[iSys] += -JuMP.shadow_price(CP[:rampcap][iSys])
                                    #MyCFC[iSys] += -JuMP.shadow_price(CP[:capproc][iSys])
                                end
                                for iArea = 1:NH2Area
                                    MyCFHD[iArea] += -JuMP.shadow_price(CP[:h2discharge][iArea])
                                    MyCFHS[iArea] += -JuMP.shadow_price(CP[:h2storage0][iArea,1])
                                end
                            end
                            MyCFP = MyCFP/Float64(NWind)
                            #MyCFR = MyCFR/Float64(NWind)
                            #MyCFC = MyCFC/Float64(NWind)
                            MyCFHD = MyCFHD/Float64(NWind)
                            MyCFHS = MyCFHS/Float64(NWind)
                            ObjCost = ObjCost/Float64(NWind)
                        else
                            optimize!(CP)
                            if primal_status(CP) != MOI.FEASIBLE_POINT 
                                write_to_file(CP, "CP_err.lp")
                                termstat = termination_status(CP)
                                error(println("Solver terminated with status $termstat for costprob in forward iteration (stage,scen): ",t," ",aScen))
                            end
                            ObjCost += JuMP.objective_value(CP)
                            for iSys = 1:NHSys
                                MyCFP[iSys] += -JuMP.shadow_price(CP[:hydprod][iSys])
                                #MyCFR[iSys] += -JuMP.shadow_price(CP[:rampcap][iSys])
                                #MyCFC[iSys] += -JuMP.shadow_price(CP[:capproc][iSys])
                            end
                            for iArea = 1:NH2Area
                                MyCFHD[iArea] += -JuMP.shadow_price(CP[:h2discharge][iArea])
                                MyCFHS[iArea] += -JuMP.shadow_price(CP[:h2storage0][iArea,1])
                            end
                        end
        
                        for iSys = 1:NHSys
                            ResArray[jCnt] = MyCFP[iSys]; jCnt += 1
                            ResArray[jCnt] = MyCFR[iSys]; jCnt += 1
                            ResArray[jCnt] = MyCFC[iSys]; jCnt += 1
                        end
                        for iArea = 1:NH2Area
                            ResArray[jCnt] = MyCFHD[iArea]; jCnt += 1
                            ResArray[jCnt] = MyCFHS[iArea]; jCnt += 1
                        end
                        MyCFRHS = ObjCost+sum(MyCFP[iSys]*HydroProd[iSys] for iSys=1:NHSys)
                        if NH2Area >= 1
                            MyCFRHS = MyCFRHS + sum(MyCFHD[iArea]*H2Discharge[iArea]+MyCFHS[iArea]*H2Init[iArea] for iArea=1:NH2Area)
                        end
                        #MyCFRHS = ObjCost+sum((MyCFP[iSys]*HydroProd[iSys]+MyCFR[iSys]*RampRate[iSys]+MyCFC[iSys]*CapProc[iSys]) for iSys=1:NHSys)
                        ResArray[jCnt] = MyCFRHS; jCnt += 1
        
                        Gap = ObjCost-Beta
                        if c < parameters.Control.CCMaxIter 
                            @constraint(SP_FORW,SP_FORW[:beta]
                                        +sum(MyCFP[iSys]*SP_FORW[:prod][iSys] for iSys=1:NHSys)+sum(MyCFR[iSys]*SP_FORW[:ramp][iSys] for iSys=1:NHSys)+sum(MyCFC[iSys]*SP_FORW[:cap][iSys] for iSys=1:NHSys) 
                                        +sum(MyCFHD[iArea]*SP_FORW[:h2dis][iArea] for iArea=1:NH2Area) >= MyCFRHS-sum(MyCFHS[iArea]*H2Init[iArea] for iArea=1:NH2Area))
                        end
        
                        #@printf("%s %4.0f %12.2f %12.2f %12.2f %12.2f \n","(Obj,Cost,Beta,Gap): ",c,ObjRel,ObjCost,Beta,Gap)
                    end
                else 
                    #Use previously computed cost function approximation
                    jCnt += parameters.Control.CCMaxIter*(3*NHSys+2*NH2Area+1)

                    optimize!(SP_FORW)
                    if primal_status(SP_FORW) != MOI.FEASIBLE_POINT 
                        write_to_file(SP_FORW, "SP_FORW_err.lp")
                        termstat = termination_status(SP_FORW)
                        error(println("Solver terminated with status $termstat for stageprob in forward iteration (stage,scen): ",t," ",aScen))
                    end

                end
                dTF = dTF+(time_ns()-t1)

            else #full stage problem formulation
                jCnt += parameters.Control.CCMaxIter*(3*NHSys+2*NH2Area+1)
                t1 = time_ns()

                optimize!(SP_FORW)
                if primal_status(SP_FORW) != MOI.FEASIBLE_POINT 
                    write_to_file(SP_FORW, "SP_FORW_FULL_err.lp")
                    termstat = termination_status(SP_FORW)
                    error(println("Solver terminated with status $termstat for stageprob in forward iteration (stage,scen): ",t," ",aScen))
                end
                dTF = dTF+(time_ns()-t1)
            end
    
            for iSys = 1:NHSys
                SimulatedStateTraj[iSys,iScen,t] = JuMP.value(SP_FORW[:res][iSys,end])
                ResArray[jCnt] = SimulatedStateTraj[iSys,iScen,t]; jCnt += 1
            end
            for iArea = 1:NH2Area
                if model.H2Data.Areas[iArea].LStrategic 
                    SimulatedH2Traj[iArea,iScen,t] = JuMP.value(SP_FORW[:h2res][iArea,end])
                else
                    SimulatedH2Traj[iArea,iScen,t] = 0.0
                end

                ResArray[jCnt] = SimulatedH2Traj[iArea,iScen,t]; jCnt += 1
            end
    
            if t < NStage
                SimulatedCost[iScen] += (JuMP.objective_value(SP_FORW)-JuMP.value(SP_FORW[:alpha]))
            else
                SimulatedCost[iScen] += JuMP.objective_value(SP_FORW)
            end
    
            if  t == 1
                SumBackward += (1/parameters.Control.NScen)*JuMP.objective_value(SP_FORW)
            end
        end 
    end
    for iScen = 1:NScenPerCluster
       ResArray[jCnt] = SimulatedCost[iScen]; jCnt += 1
    end
    ResArray[jCnt] = SumBackward
    
    ResArray[NDim] = dTF
 
    return ResArray
 end
 
