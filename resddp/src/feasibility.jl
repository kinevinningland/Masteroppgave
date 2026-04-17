#This function computes feasibility cuts for aggregated hydropower systems
#- Evaluate if combinations of discrete states for the aggregated system are feasible.  
#- Create feasibility cuts for non-feasible combinations

function feasibility(model::Model, inflow_model::InflowModel, parameters::Parameters, datapath::String; optimizer=JuMP.optimizer_with_attributes(Clp.Optimizer, "SolveType" => 0, "PresolveType" => 1, "LogLevel" => 0))::Vector{FeasibilitySpace}

   CTI = parameters.Time
   CNS = parameters.Constants
   CDI = parameters.Discrete
   CTR = parameters.Control
   CAGR = parameters.Aggregation

   SumProd = sum(model.HSys[iSys].MaxProd for iSys=1:model.NHSys)
   ProdShare = zeros(Float64,model.NHSys)
   for iSys = 1:model.NHSys
      ProdShare[iSys] = model.HSys[iSys].MaxProd/SumProd
   end

   NCutDim = CDI.NInfPt*CDI.NResInitPt*CDI.NResEndPt*CDI.NEnPt*CDI.NRampPt*CDI.NCapPt
   NEntries = 7 #gamv, game, gamr, gamc, kappav, kappai, RHS

   NFeasWeeks = 1
   if CTR.LFeasPerStage
      NFeasWeeks = CTI.NWeek
   end
   FeasSpaces = Vector{FeasibilitySpace}(undef,NFeasWeeks)

   #Define hydroprob ramping schedule
   StartRamp = EndRamp = 0
   if CTI.NK == 28; StartRamp = 2; EndRamp = 3; end 
   if CTI.NK == 42; StartRamp = 2; EndRamp = 5; end 
   if CTI.NK == 56; StartRamp = 3; EndRamp = 6; end 
   if CTI.NK == 84; StartRamp = 4; EndRamp = 9; end 
   if CTI.NK == 168; StartRamp = 8; EndRamp = 18; end 

   for iWeek = 1:NFeasWeeks
      AMAT = zeros(Float64,model.NHSys,NCutDim,NEntries)
      NFeasCut = zeros(Int,model.NHSys)
      println("week ",iWeek)

      for iSys = 1:model.NHSys
         println("   system ",iSys)

         MyHSys = model.HSys[iSys]
         if MyHSys.MaxRes < CAGR.ResCutoff || MyHSys.MaxProd < CAGR.ProdCutoff 
             continue
         end

         iArea = MyHSys.AreaNo
         MyAHData = model.AHData[iArea]
         MyUSMod = model.USModSys[iSys]

         NMod = MyHSys.NMod
         #println("System ",iSys, " with ",NMod," modules")

         MinResCurrWeek = sum(MyAHData.EffSea[MyHSys.ModNo[iMod]]*CNS.MAGEFF2GWH*model.DMData.MaMin[iArea,MyHSys.ModNo[iMod],iWeek] for iMod=1:NMod)

         deltaInflow = 2*mean(model.InfReg[iSys,iWeek,1:CTI.NInflowYear])/Float64(CDI.NInfPt-1)      # [GWh]
         deltaResInit = (MyHSys.MaxRes-MinResCurrWeek)/Float64(CDI.NResInitPt-1)                     # [GWh]
         deltaResEnd =  (MyHSys.MaxProd-MyHSys.MinProd[iWeek])/Float64(CDI.NResEndPt-1)              # [GWh/stage]
         deltaEnergy = (MyHSys.MaxProd-MyHSys.MinProd[iWeek])/Float64(CDI.NEnPt-1)                   # [GWh/stage]
         deltaRamp =  CTI.WeekFrac*(MyHSys.MaxProd-MyHSys.MinProd[iWeek])/Float64(CDI.NRampPt-1)     # [GWh/step]
         deltaCap = 0.30*CTI.WeekFrac*(MyHSys.MaxProd-MyHSys.MinProd[iWeek])/Float64(CDI.NCapPt-1)   # [GWh/step/week]

         MyRegFrac = model.RegFrac[iSys,iWeek]
         DetProb = HydProb.Build(iWeek,iSys,iArea,MyHSys,NMod,MyAHData,MyUSMod,model.DMData,CTI.NK,CNS.InfUB,CTI.DT,CNS,CNS.M3S2MM3,CNS.GWH2MAGEFF,
                                 CNS.MAGEFF2GWH,CNS.MW2GWHWEEK,CTI.WeekFrac,MyRegFrac,StartRamp,EndRamp,optimizer)

         NPtTot = 0
         rCnt = 0
         for infCnt = 1:(CDI.NInfPt-1)
            InflowSys = infCnt*deltaInflow
            JuMP.set_normalized_rhs(DetProb[:inflow],InflowSys)

            for resInitCnt = 1:(CDI.NResInitPt-1)
               ResStart = MinResCurrWeek + (resInitCnt-0.5)*deltaResInit
               JuMP.set_normalized_rhs(DetProb[:resinit],ResStart)

               for resEndCnt = 1:CDI.NResEndPt
                  ResEnd = max(ResStart-(MyHSys.MaxProd-MyHSys.MinProd[iWeek]),0.0)+(resEndCnt-0.5)*deltaResEnd
                  JuMP.set_normalized_rhs(DetProb[:resreq],ResEnd)

                  for enCnt = 1:(CDI.NEnPt-1)
                     EnergyLevel = (enCnt-0.5)*deltaEnergy
                     JuMP.set_normalized_rhs(DetProb[:enreq],EnergyLevel)

                     for rampCnt = 1:(CDI.NRampPt-1)
                        RampLevel = (rampCnt-0.5)*deltaRamp
                        JuMP.set_normalized_rhs(DetProb[:rampreq],RampLevel)

                        for capCnt = 1:(CDI.NCapPt-1)
                           CapLevel = (capCnt-0.5)*deltaCap
                           JuMP.set_normalized_rhs(DetProb[:capreq],CapLevel)

                           #println(infCnt,"  ",resInitCnt,"  ",resEndCnt,"  ",enCnt,"  ",rampCnt,"  ",capCnt)
                           optimize!(DetProb)
                           if primal_status(DetProb) != MOI.FEASIBLE_POINT 
                              write_to_file(DetProb, "DetProb_err.lp")
                              termstat = termination_status(DetProb)
                              error(println("Solver terminated with status $termstat for system: ",iSys))
                           end
                           obj = JuMP.objective_value(DetProb)
                           NPtTot += 1
                           

                           #shadow_price gives negative values, adjust signs to reflect gain/loss
                           if obj > CNS.FeasTol
                              #write_to_file(DetProb, string(string("DetProb",rCnt),".lp"))
                              rCnt += 1

                              gamv = max(-JuMP.shadow_price(DetProb[:resreq]),0.0) #Higher is worse [mu/GWh]
                              game = max(-JuMP.shadow_price(DetProb[:enreq]),0.0)      #Higher is worse [mu/GWh/week]
                              gamr = max(-JuMP.shadow_price(DetProb[:rampreq]),0.0)    #Higher is worse [mu/GWh/step] 
                              gamc = max(-JuMP.shadow_price(DetProb[:capreq]),0.0)     #Higher is worse [mu/GWh/step/week]
                              kappav = JuMP.shadow_price(DetProb[:resinit])   #Higher is better [mu/GWh]
                              kappai = JuMP.shadow_price(DetProb[:inflow])    #Higher is better [mu/GWh]

                              if abs(gamv) > 1.001 println("ERROR: gamv > 1") end
                              if abs(game) > 1.001 println("ERROR: game > 1") end
                              if abs(gamr) > 1.001 println("ERROR: gamr > 1") end
                              if abs(gamc) > 1.001 println("ERROR: gamc > 1") end
                              if abs(kappai) > 1.001 println("ERROR: kappav > 1") end
                              if abs(kappav) > 1.001 println("ERROR: kappai > 1") end

                              CRHSFeas = -obj+gamv*ResEnd+game*EnergyLevel+gamr*RampLevel+gamc*CapLevel+kappav*ResStart+kappai*InflowSys

                              NDigs = 4
                              AMAT[iSys,rCnt,1] = round(gamv,digits=NDigs) 
                              AMAT[iSys,rCnt,2] = round(game,digits=NDigs) 
                              AMAT[iSys,rCnt,3] = round(gamr,digits=NDigs) 
                              AMAT[iSys,rCnt,4] = round(gamc,digits=NDigs) 
                              AMAT[iSys,rCnt,5] = round(kappav,digits=NDigs)
                              AMAT[iSys,rCnt,6] = round(kappai,digits=NDigs)
                              AMAT[iSys,rCnt,7] = round(CRHSFeas,digits=NDigs)

                              #if gamv+game > 1.99 && gamr > 0.0019 && CRHSFeas > -120.51 && CRHSFeas < -120.50
                                 #ResTab = Results.Initialize(NMod,NK)
                                 #Results.Store(ResTab,DetProb,AHData,iSys,NMod,NK,MAGEFF2GWH)
                                 #Results.Print(dataset,ResTab,NMod,NK,CTI)
                                 #println("Inflow:   ",InflowSys)
                                 #println("ResStart: ",ResStart)
                                 #println("ResEnd:   ",ResEnd)
                                 #println("Energy:   ",EnergyLevel)
                                 #println("Ramp:     ",RampLevel)
                                 #println("Capacity: ",CapLevel)
                                 #println("CRHSFeas: ",CRHSFeas)
                                 #println("gamv:     ",gamv)
                                 #println("game:     ",game)
                                 #println("gamr:     ",gamr)
                                 #println("gamc:     ",gamc)
                                 #println("kappav:   ",kappav)
                                 #println("kappai:   ",kappai)
                                 #write_to_file(DetProb, string(string("DetProb",rCnt),".lp"))
                                 #@goto exit_label
                              #end
                           end
                        end
                     end
                  end
               end
            end
         end
            
         AMATRED = unique(round.(AMAT[iSys,1:rCnt,1:NEntries],digits=3),dims=1) #or dims=2?
         NFeasCut[iSys] = size(AMATRED)[1]
         fill!(AMAT[iSys,:,:],0.0)
         AMAT[iSys,1:NFeasCut[iSys],1:NEntries] = AMATRED
         
         #println("... ",NPtTot," evaluated, ",rCnt," binding, " ,NFeasCut[iSys]," unique cuts")
      end

      NFeasCutMax = maximum(NFeasCut)
      FCC = zeros(Float64,model.NHSys,NFeasCutMax,NEntries)
      FCC = AMAT[:,1:NFeasCutMax,:]

      FeasSpaces[iWeek] = FeasibilitySpace(iWeek,FCC,NFeasCut,NFeasCutMax)
   end

   return FeasSpaces

   #@label exit_label
end

