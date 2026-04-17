function SolveBWD(model::Model, inflow_model::InflowModel, feas_space::FeasibilitySpace,
    parameters::Parameters, strategy::Strategy, init_val::InitialValues, SS::SampleScen,
    t_bwd, iWeek, bProb, scenario_range::UnitRange{Int}, optimizer)

    NScenPerCluster = length(scenario_range)
    NStage = parameters.Control.NStage
    NHSys = model.NHSys
    NWeek = parameters.Time.NWeek
    NH2Area = model.H2Data.NArea 
    HSys = model.HSys

    StartScen = scenario_range[begin]
    EndScen = scenario_range[end]
 
    NDimCut = 2*NHSys+NH2Area+1
    NDim = NScenPerCluster*NDimCut+1
    ResArray = zeros(Float64,NDim)
    ResInitBwd = zeros(Float64,NHSys)
    H2InitBwd = zeros(Float64,NH2Area)
 
    if parameters.Control.LCostApprox
        SP_BWD = StageProb.Build(t_bwd, iWeek, inflow_model, NHSys, HSys,
            model.AreaSys, model.AMData, model.EV,
            strategy.CCR, strategy.CCH, strategy.CCI, strategy.CFP, strategy.CFR,
            strategy.CFC, strategy.CFHD, strategy.CFHS, strategy.CFRHS, parameters.Constants, strategy.NCut, strategy.NCostCut,
            model.NArea, parameters.Time, t_bwd==NStage,
            parameters.Control.LFeasSpace, feas_space.NFeasCut, feas_space.FCC, model.H2Data, optimizer)
    else
        wYear = 1
        MyWPData = model.WPData[1:model.NArea, wYear, iWeek, 1:parameters.Time.NK]
        SP_BWD = StageProbFull.Build(t_bwd, iWeek, inflow_model, NHSys, HSys,
            model.NAreaSys, model.AreaSys, model.AMData,
            model.MCon, model.EV, strategy.CCR,strategy.CCH,strategy.CCI,
            parameters.Constants, strategy.NCut, model.NArea, model.NLine,
            model.LineCap,model.LineLoss, parameters.Time, t_bwd==NStage,
            parameters.Control.LFeasSpace, feas_space.NFeasCut, feas_space.FCC,
            parameters.Control.CapReqFrac, MyWPData, parameters.Control.LDemandResponse, model.DRData, model.H2Data, optimizer)
    end
       
    cCnt = 1
    dTB = 0.0
    for  aScen = StartScen:EndScen
        ResInitBwd[1:NHSys] = strategy.StateTraj[1:NHSys,aScen,t_bwd-1]
        H2InitBwd[1:NH2Area] = strategy.H2Traj[1:NH2Area,aScen,t_bwd-1]
        Zstate = SS.Zscen[aScen,t_bwd-1,1:inflow_model.NSer]
        for iSys = 1:NHSys
            JuMP.set_normalized_rhs(SP_BWD[:rstate][iSys],ResInitBwd[iSys])
            JuMP.set_normalized_rhs(SP_BWD[:zstate][iSys],Zstate[iSys])
        end
        for iArea = 1:NH2Area
            JuMP.set_normalized_rhs(SP_BWD[:h2state][iArea],H2InitBwd[iArea])
        end
 
        #Consider sampling outside scenario-loop
        SampleYears = collect(1:inflow_model.NResid)
        sRes = sample(SampleYears,parameters.Control.NBranch,replace=false)
        for iBranch = 1:parameters.Control.NBranch
            jCnt = 1
            Eps = inflow_model.Resid[1:inflow_model.NSer,iWeek,sRes[iBranch]]
            DevEpsVec = inflow_model.InflowSDev[1:inflow_model.NSer,iWeek].*Eps
            for iSys = 1:NHSys
                JuMP.set_normalized_rhs(SP_BWD[:inflow][iSys],DevEpsVec[iSys]+inflow_model.InflowMean[iSys,iWeek])
            end
            if t_bwd < NStage
                adjust = transpose(strategy.CCI[1:inflow_model.NSer,t_bwd,1:strategy.NCut])*Eps[1:inflow_model.NSer]
                for iCut = 1:strategy.NCut
                    JuMP.set_normalized_rhs(SP_BWD[:cut][iCut],strategy.CRHS[t_bwd,iCut]+adjust[iCut])
                end
            end
    
            #Solve LP problem
            t1 = time_ns()
            optimize!(SP_BWD)
            dTB = dTB+(time_ns()-t1)

            #ncbind = 0
            #sumbind = 0.0
            #for c=1:strategy.NCostCut
                #dual = JuMP.shadow_price(SP_BWD[:costcut][c])
                #if (abs(dual) > 1.0E-3)
                    ##@printf("%4.0f %8.3f \n",c,dual)
                    #ncbind += 1
                    #sumbind += dual
                #end
            #end
            #@printf("%4.0f %4.0f %4.0f %4.0f %4.0f %8.3f \n",t_bwd,aScen,iBranch,strategy.NCostCut,ncbind,sumbind)
    
            if dual_status(SP_BWD) != MOI.FEASIBLE_POINT
                write_to_file(SP_BWD, "SPB_err.lp")
                termstat = termination_status(SP_BWD)
                error(println("Solver terminated with status $termstat in backward iteration (stage,scen,branch,sRes): ",t_bwd," ",aScen," ",iBranch," ",sRes[iBranch]))
            end
    
            #Find contributions to cut parameters
            for iSys = 1:NHSys
                ResArray[(cCnt-1)*NDimCut+jCnt] += bProb*JuMP.shadow_price(SP_BWD[:rstate][iSys]); jCnt += 1
                ResArray[(cCnt-1)*NDimCut+jCnt] += bProb*JuMP.shadow_price(SP_BWD[:zstate][iSys]); jCnt += 1
            end
            for iArea = 1:NH2Area
                if model.H2Data.Areas[iArea].LStrategic
                    ResArray[(cCnt-1)*NDimCut+jCnt] += bProb*JuMP.shadow_price(SP_BWD[:h2state][iArea])
                end 
                jCnt += 1
            end
    
            rhs = JuMP.objective_value(SP_BWD)
            rhsDelta = 0.0
            for iSys = 1:NHSys
                rhsDelta += ResInitBwd[iSys]*JuMP.shadow_price(SP_BWD[:rstate][iSys]) 
                rhsDelta += Zstate[iSys]*JuMP.shadow_price(SP_BWD[:zstate][iSys])
            end
            for iArea = 1:NH2Area
                if model.H2Data.Areas[iArea].LStrategic
                    rhsDelta += H2InitBwd[iArea]*JuMP.shadow_price(SP_BWD[:h2state][iArea]) 
                end
            end
            ResArray[(cCnt-1)*NDimCut+jCnt] += bProb*(rhs-rhsDelta)
        end
        cCnt += 1
    end
    ResArray[NDim] = dTB
    return ResArray
 end
 
