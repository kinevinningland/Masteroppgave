function train!(strategy::Strategy, init_val::InitialValues, model::Model, inflow_model::InflowModel, feas_spaces::Vector{FeasibilitySpace}, parameters::Parameters; datapath::String=".",
    optimizer=JuMP.optimizer_with_attributes(Clp.Optimizer, "SolveType" => 0, "PresolveType" => 1, "LogLevel" => 0))
    
    #Start SDDP algorithm
    SumForward = Array{Float64}(undef, 0)
    SumBackward = Array{Float64}(undef, 0)
    MaxIter = strategy.MaxIter+parameters.Control.MaxIter-1

    open(joinpath(datapath, "clog.dat"), "w") do clog
        open(joinpath(datapath, "tlog.dat"), "w") do tlog
            NCluster = min(Threads.nthreads(), parameters.Control.NScen) # Never more threads than scenarios
            NScenPerCluster = Int(ceil(parameters.Control.NScen/NCluster)) # Maximum number of scenario per thread

            NH2Area = model.H2Data.NArea 
            NScen = parameters.Control.NScen
            NStage = parameters.Control.NStage
            NWeek = parameters.Time.NWeek

            NCut = strategy.NCut
            NCostCut = strategy.NCostCut

            if NCut == 0
                if parameters.Control.LWindStoch
                    strategy.WindYears = sample(collect(1:parameters.Time.NInflowYear),parameters.Control.NWindScen,replace=false)
                else
                    strategy.WindYears = [1]
                end
            end

            println("VRES-years in strategy: ",strategy.WindYears)

            for it = strategy.MaxIter:MaxIter
                dTF1 = dTF2 = dTB1 = dTB2 = 0.0
                println("Iteration ",it)
                t_iter_fwd = time_ns()
                
                # Forward iteration
                push!(SumBackward, 0)
                ForwardProfit = zeros(Float64,NScen)
                SS = SampleScenario(NScen, NStage, NWeek, inflow_model, parameters.Control.LExtreme)
            
                t1 = time_ns()

                results = Array{Array{Float64}}(undef, NCluster)
                clusterScenarios = Array{UnitRange}(undef, NCluster)
                Threads.@threads for sCluster = 1:NCluster
                    start_scen = (sCluster-1) * NScenPerCluster + 1
                    if start_scen <= NScen
                        end_scen = min(sCluster * NScenPerCluster, NScen)
                        clusterScenarios[sCluster] = start_scen:end_scen
                        results[sCluster] = SolveFWD(model, inflow_model, feas_spaces, parameters, strategy, init_val, SS, clusterScenarios[sCluster],optimizer)
                    else
                        clusterScenarios[sCluster] = 1:0
                    end
                end
                dTF2 = dTF2+(time_ns()-t1)
        
                for iCluster = 1:NCluster
                    jCnt = 1
                    if length(clusterScenarios[iCluster]) > 0 # only run if scenarios has been assigned to corresponding thread
                        for t = 1:NStage
                            for aScen = clusterScenarios[iCluster]
                                for cc = 1:parameters.Control.CCMaxIter
                                    cloc = parameters.Control.CCMaxIter*((it-1)*NScen + (aScen-1)) + cc
                                    if parameters.Control.LCostApproxNewCuts
                                        for iSys = 1:model.NHSys
                                            strategy.CFP[iSys,t,cloc] = results[iCluster][jCnt]; jCnt += 1
                                            strategy.CFR[iSys,t,cloc] = results[iCluster][jCnt]; jCnt += 1
                                            strategy.CFC[iSys,t,cloc] = results[iCluster][jCnt]; jCnt += 1
                                        end
                                        for iArea = 1:NH2Area
                                            strategy.CFHD[iArea,t,cloc] = results[iCluster][jCnt]; jCnt += 1
                                            strategy.CFHS[iArea,t,cloc] = results[iCluster][jCnt]; jCnt += 1
                                        end
                                        strategy.CFRHS[t,cloc] = results[iCluster][jCnt]; jCnt += 1
                                    else
                                        jCnt += 3*model.NHSys+2*NH2Area+1
                                    end
                                end
                                for iSys = 1:model.NHSys
                                    strategy.StateTraj[iSys,aScen,t] = results[iCluster][jCnt]; jCnt += 1
                                end
                                for iArea = 1:NH2Area
                                    strategy.H2Traj[iArea,aScen,t] = results[iCluster][jCnt]; jCnt += 1
                                end
                            end
                        end
                        for aScen = clusterScenarios[iCluster]
                            ForwardProfit[aScen] = results[iCluster][jCnt]; jCnt += 1
                        end
                        SumBackward[end] += results[iCluster][jCnt]
                        dTF1 += results[iCluster][jCnt+1]
                    end
                end
        
                push!(SumForward, mean(ForwardProfit))
                sumd = std(ForwardProfit)
                sigma = sqrt(sumd/(NScen-1))
            
                gap = SumForward[end]-SumBackward[end]
                println(" SumBack $(@sprintf("%.0f",SumBackward[end])) SumForw $(@sprintf("%.0f",SumForward[end])) gap $(@sprintf("%.2f",gap))")
                @printf(clog,"%.0f %12.2f %12.2f \n",it,SumBackward[end],SumForward[end])
                if it == MaxIter || abs(gap) < parameters.Control.ConvEps
                    if it == MaxIter println("Maximum no. iterations") end
                    if abs(gap) < parameters.Control.ConvEps println("Converged!") end
                    break
                end
                time_iter_fwd = (time_ns()-t_iter_fwd)*1.0E-9
            
                NCut = NCut+NScen
                strategy.NCut = NCut
                if parameters.Control.LCostApproxNewCuts
                    NCostCut = NCostCut+parameters.Control.CCMaxIter*NScen
                end
                strategy.NCostCut = NCostCut
        
                # Backward iteration
                t_iter_bwd = time_ns()
                NDimCut = 2*model.NHSys+NH2Area+1
                bProb = 1.0/Float64(parameters.Control.NBranch)
                for t = NStage:-1:2
                    iWeek = mod1(t,NWeek)
                    fWeek = iWeek; if !parameters.Control.LFeasPerStage; fWeek = 1; end
                    t_bwd = t
                    t1 = time_ns()
                    Threads.@threads for sCluster = 1:NCluster
                        if length(clusterScenarios[sCluster]) > 0
                            results[sCluster] = SolveBWD(model, inflow_model, feas_spaces[fWeek], parameters, strategy, init_val, SS, t_bwd, iWeek, bProb, clusterScenarios[sCluster], optimizer)
                        end
                    end
                    dTB2 = dTB2+(time_ns()-t1)
                
                    for iCluster = 1:NCluster
                        if length(clusterScenarios[iCluster]) > 0
                            for (iScen, aScen) = enumerate(clusterScenarios[iCluster])
                                jCnt = 1
                                cloc = (it-1)*NScen + aScen
                                for iSys = 1:model.NHSys
                                    strategy.CCR[iSys,t-1,cloc] = results[iCluster][(iScen-1)*NDimCut+jCnt]; jCnt += 1
                                    strategy.CCI[iSys,t-1,cloc] = results[iCluster][(iScen-1)*NDimCut+jCnt]; jCnt += 1
                                end
                                for iArea = 1:NH2Area
                                    strategy.CCH[iArea,t-1,cloc] = results[iCluster][(iScen-1)*NDimCut+jCnt]; jCnt += 1
                                end
                                strategy.CRHS[t-1,cloc] = results[iCluster][(iScen-1)*NDimCut+jCnt]
                            end
                            dTB1 += results[iCluster][length(clusterScenarios[iCluster])*NDimCut+1]
                        end
                    end
                end
                time_iter_bwd = (time_ns()-t_iter_bwd)*1.0E-9
                @printf("%s %6.2f %6.2f %6.2f %s %6.2f %6.2f %6.2f \n",
                        " Time (FWD): ",dTF1/Float64(NCluster)*1.0E-9,dTF2*1.0E-9,time_iter_fwd,
                        " Time (BWD): ",dTB1/Float64(NCluster)*1.0E-9,dTB2*1.0E-9,time_iter_bwd)
                @printf(tlog,"%4.0f %6.2f %6.2f %6.2f %6.2f \n",it,dTF2*1.0E-9,time_iter_fwd,dTB2*1.0E-9,time_iter_bwd)
                #time_iter = (time_ns()-t_iter)*1.0E-9
            end
        end # close tlog
    end # close clog
    
    strategy.MaxIter = MaxIter
    return SumForward, SumBackward
end
