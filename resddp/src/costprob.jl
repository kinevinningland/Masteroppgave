#Aggregated stage problem 
module CostProb
   using JuMP 
   using MathOptInterface
   const MOI = MathOptInterface

   function Build(t,iWeek,NHSys,HSys,NAreaSys,AreaSys,AMData,MCon,CNS,NArea,NLine,LineCap,LineLoss,CTI,MyWPData,LDemandResponse,DR,H2Data,optimizer)
      M = Model(optimizer)

      NK = CTI.NK
      NH2Area = H2Data.NArea
      WeekFrac = CTI.WeekFrac

      NLoadRecStep = DR.NLoadRecStep
      MaxLoadRec = DR.LoadRec[NLoadRecStep]

      @variable(M,WeekFrac*HSys[iSys].MinProd[iWeek] <= prod[iSys=1:NHSys,k=1:NK] <= WeekFrac*HSys[iSys].MaxProd,base_name="prod")  # GWh/step
      @variable(M,0.0 <= cap[iSys=1:NHSys] <= 0.5*WeekFrac*(HSys[iSys].MaxProd-HSys[iSys].MinProd[iWeek]),base_name="cap")          # GWh/step/stage, reserve capacity from hydro
      @variable(M,0.0 <= ramp[iSys=1:NHSys] <= WeekFrac*(HSys[iSys].MaxProd-HSys[iSys].MinProd[iWeek]),base_name="ramp")            # GWh/step/stage, ramping capability for hydro
      @variable(M,min(0.0,WeekFrac*AMData[iArea].MSData[iMark].Capacity[iWeek]) 
                <= mark[iArea=1:NArea,iMark=1:AMData[iArea].NMStep,k=1:NK] 
                <= WeekFrac*max(0.0,AMData[iArea].MSData[iMark].Capacity[iWeek]), base_name="mark")                                 # GWh/step Purchase(+) or sales(-) from market 
      @variable(M,0.0 <= etran[iLine=1:NLine,k=1:NK] <= WeekFrac*LineCap[iLine], base_name="etran")                                 # GWh Exchange 
      @variable(M,0.0 <= wprod[iArea=1:NArea,k=1:NK] <= CNS.Big,base_name="wprod")                                                  # GWh/step
      @variable(M,0.0 <= rat[iArea=1:NArea,k=1:NK] <= CNS.Big,base_name="rat")                                                      # GWh/step
      @variable(M,0.0 <= h2dis[iArea=1:NH2Area,k=1:NK] <= H2Data.Areas[iArea].MaxDis, base_name="h2dis")                            # GWh/step
      @variable(M,0.0 <= h2chg[iArea=1:NH2Area,k=1:NK] <= H2Data.Areas[iArea].MaxDis, base_name="h2chg")                            # GWh/step
      @variable(M,0.0 <= h2res[iArea=1:NH2Area,k=1:NK] <= H2Data.Areas[iArea].MaxRes, base_name="h2res")                            # GWh

      #Demand response variables
      if LDemandResponse      
         @variable(M,-CNS.Big <= dr_tot[iArea=1:NArea, k=1:NK] <= CNS.Big, base_name="dr_tot")
         @variable(M,0.0 <= dr_up[iArea=1:NArea, k=1:NK] <= CNS.Big, base_name="dr_up")
         @variable(M,0.0 <= dr_dn[iArea=1:NArea, kk=1:(2*MaxLoadRec+1), k=1:NK] <= CNS.Big, base_name="dr_dn")
         @variable(M,0.0 <= dr_dir[iArea=1:NArea, k=1:NK] <= 1.0, base_name="dr_dir")
      end

      #Min Cost [10E3 EUR]
      @objective(M,MathOptInterface.MIN_SENSE,
                 sum(AMData[iArea].MSData[iMark].Price[iWeek]*mark[iArea,iMark,k] for iArea=1:NArea for iMark=1:AMData[iArea].NMStep for k=1:NK)
                 +sum(CNS.CRat*rat[iArea,k] for iArea=1:NArea for k=1:NK))
     
      #POWER BALANCE [GWh/step]
      @constraint(M,pbal[iArea=1:NArea,k=1:NK],
                  + sum(prod[AreaSys[iArea,iSys],k] for iSys=1:NAreaSys[iArea])+rat[iArea,k]+sum(mark[iArea,iMark,k] for iMark=1:AMData[iArea].NMStep) + wprod[iArea,k] 
                  -sum(etran[MCon[iArea].LIndxOut[iLine],k] for iLine=1:MCon[iArea].NCon)
                  +sum((1.0-LineLoss[MCon[iArea].LIndxIn[iLine]])*etran[MCon[iArea].LIndxIn[iLine],k] for iLine=1:MCon[iArea].NCon)
                  - (LDemandResponse ? dr_tot[iArea,k] : 0.0) #Include dr_tot in power balance only if LDemandResponse is true
                  - (H2Data.Ind[iArea] > 0 ? h2chg[H2Data.Ind[iArea],k] : 0.0) 
                  + (H2Data.Ind[iArea] > 0 ? h2dis[H2Data.Ind[iArea],k] : 0.0) 
                  == sum(AMData[iArea].MLData[iLoad].Load[iWeek,k] for iLoad=1:AMData[iArea].NLoad))

      #RAMP UP [GWh/step]
      @constraint(M,rampup[iSys=1:NHSys,k=2:NK],prod[iSys,k] - prod[iSys,k-1] <= ramp[iSys])

      #RAMP DOWN [GWh/step]
      @constraint(M,rampdn[iSys=1:NHSys,k=2:NK],prod[iSys,k] - prod[iSys,k-1] >= -ramp[iSys])

      #PRODUCTION REQUIREMENT (MU-E) [GWh]
      @constraint(M,hydprod[iSys=1:NHSys],sum(prod[iSys,k] for k=1:NK) <= 0.0)
      
      #RAMPING REQUIREMENT (MU-R) [GWh/step/stage]
      @constraint(M,rampcap[iSys=1:NHSys], ramp[iSys] <= 0.0)

      #CAPACITY REQUIREMENT (MU-C) [GWh/step/stage]
      @constraint(M,capproc[iSys=1:NHSys], cap[iSys] <= 0.0)

      #H2 INITIAL STORAGE (MU-HS) [GWh]
      @constraint(M,h2storage0[iArea=1:NH2Area,k=[1]],h2res[iArea,k]-(1.0-H2Data.Areas[iArea].CompLoss)*h2chg[iArea,k]+h2dis[iArea,k] == 0.0)

      #H2 STORAGE [GWh]
      @constraint(M,h2storage[iArea=1:NH2Area,k=2:NK],h2res[iArea,k]-(1.0-H2Data.Areas[iArea].CompLoss)*h2chg[iArea,k]+h2dis[iArea,k]-h2res[iArea,k-1] == 0.0)

      #H2 DISCHARGE REQUIREMENT (MU-HD) [GWh]
      @constraint(M,h2discharge[iArea=1:NH2Area],sum(h2dis[iArea,k]-h2chg[iArea,k] for k=1:NK) <= 0.0)


      if LDemandResponse      
         #Demand balance per time step
         @constraint(M,drbal[iArea=1:NArea, k=1:NK],dr_tot[iArea, k]-dr_up[iArea, k]+sum(dr_dn[iArea, kk,k] for kk = max(1,MaxLoadRec+2-k):min(2*MaxLoadRec+1,MaxLoadRec+1+NK-k)) == 0.0)
         #Demand balance between up and down with load recovery
         @constraint(M,drbal2[iArea=1:NArea, k=1:NK],dr_up[iArea, k]-sum(dr_dn[iArea, MaxLoadRec+1+k-kk, kk] for kk = max(1,k-MaxLoadRec):min(NK,k+MaxLoadRec)) == 0.0)
         #Maximum demand down constraint
         @constraint(M,drmaxdn[iStep=1:NLoadRecStep, iArea=1:NArea, k=1:NK],
                     sum(dr_dn[iArea, kk,k] for kk = max(MaxLoadRec-DR.LoadRec[iStep]+1,2*MaxLoadRec-DR.LoadRec[iStep]+2-k):min(MaxLoadRec+DR.LoadRec[iStep]+1,DR.LoadRec[iStep]+1+NK-k))
                     <= dr_dir[iArea,k]*sum(DR.MaxDnShift[iArea,iWeek,jStep] for jStep=1:iStep))
         #Maximum demand up constraint
         @constraint(M,drmaxup[iArea=1:NArea, k=1:NK],dr_up[iArea, k] <= sum(DR.MaxUpShift[iArea,iWeek,iStep] for iStep=1:NLoadRecStep)*(1.0 - dr_dir[iArea,k]))

         #Set dr_dn variables on edges to zero (not used in the constraints)
         @constraint(M,drdnzeroleft[iArea=1:NArea, k=1:MaxLoadRec], sum(dr_dn[iArea, kk, k] for kk=1:(MaxLoadRec+1-k)) == 0.0)
         @constraint(M,drdnzeroright[iArea=1:NArea, k=(NK+1-MaxLoadRec):NK], sum(dr_dn[iArea, kk, k] for kk=(MaxLoadRec+2+NK-k):(2*MaxLoadRec+1)) == 0.0)

         #Additional constraints to limit simultaneous load shifting
         #Only add for up shifting if extra constraints should be added for longest load recovery period
         if DR.LIncludeExtraConstr[NLoadRecStep]
            @constraint(M,drmaxupacc[iArea=1:NArea, k=DR.ExtraConstrFilter[NLoadRecStep]:NK],
                        sum(dr_up[iArea, l] for l=(k-DR.ExtraConstrFilter[NLoadRecStep]+1):k) <= DR.ExtraConstrFilter[NLoadRecStep]*DR.ExtraConstrSigma*sum(DR.MaxUpShift[iArea,iWeek,iStep] for iStep=1:NLoadRecStep))
         end

         ExtraConstrSteps = findall(DR.LIncludeExtraConstr)
         @constraint(M, drmaxdnacc[iStep in ExtraConstrSteps, iArea=1:NArea, k=DR.ExtraConstrFilter[iStep]:NK],
                     sum(sum(dr_dn[iArea, ll,l] for ll = max(MaxLoadRec-DR.LoadRec[iStep]+1,2*MaxLoadRec-DR.LoadRec[iStep]+2-l):min(MaxLoadRec+DR.LoadRec[iStep]+1,DR.LoadRec[iStep]+1+NK-l)) for l=(k-DR.ExtraConstrFilter[iStep]+1):k)
                     <= DR.ExtraConstrSigma*DR.ExtraConstrFilter[iStep]*sum(DR.MaxDnShift[iArea,iWeek,jStep] for jStep=1:iStep))
      end

      #WIND POWER TARGET [GWh/step]
      @constraint(M,wptarget[iArea=1:NArea,k=1:NK], wprod[iArea,k] <= max(MyWPData[iArea,k],0.0))


      return M
   end
end
