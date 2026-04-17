#Aggregated stage problem 
module StageProbFull
   using JuMP 
   using MathOptInterface
   const MOI = MathOptInterface

   function Build(t,iWeek,IM,NHSys,HSys,NAreaSys,AreaSys,AMData,MCon,EV,CCR,CCH,CCI,CNS,NCut,NArea,NLine,LineCap,LineLoss,CTI,LEndVal,
                  LFeasSpace,NFeasCut,FCC,CapReq,MyWPData,LDemandResponse,DR,H2Data,optimizer)

      M = Model(optimizer)

      NK = CTI.NK
      WeekFrac = CTI.WeekFrac
      NH2Area = H2Data.NArea
      
      NLoadRecStep = DR.NLoadRecStep
      MaxLoadRec = DR.LoadRec[NLoadRecStep]
      
      @variable(M,0.0 <= ner[iSys=1:NHSys,k=1:NK] <= CNS.Big,base_name="ner")                                                             # GWh
      @variable(M,0.0 <= mir[iSys=1:NHSys,k=1:NK] <= HSys[iSys].MinRes, base_name="mir")                                           # GWh
      @variable(M,0.0 <= res[iSys=1:NHSys,k=1:NK] <= HSys[iSys].MaxRes, base_name="res")                                           # GWh
      @variable(M,0.0 <= spi[iSys=1:NHSys,k=1:NK] <= CNS.Big, base_name="spi")                                                     # GWh/step
      @variable(M,WeekFrac*HSys[iSys].MinProd[iWeek] <= prod[iSys=1:NHSys,k=1:NK] <= WeekFrac*HSys[iSys].MaxProd,base_name="prod") # GWh/step
      @variable(M, 0.0 <= cap[iSys=1:NHSys] <=  0.50*WeekFrac*(HSys[iSys].MaxProd-HSys[iSys].MinProd[iWeek]),base_name="cap")      # GWh/step/stage
      @variable(M, 0.0 <= ramp[iSys=1:NHSys] <=  WeekFrac*(HSys[iSys].MaxProd-HSys[iSys].MinProd[iWeek]),base_name="ramp")         # GWh/step
      @variable(M,0.0 <= fsl[iSys=1:NHSys,c=1:NFeasCut[iSys]] <= CNS.Big, base_name="fsl")                                         # GWh
      @variable(M,zinit[iSys=1:NHSys],base_name="zinit")                                                    
      @variable(M,inf[iSys=1:NHSys],base_name="inf")                                                    
      @variable(M,rinit[iSys=1:NHSys],base_name="rinit")                                                    
      @variable(M,min(0.0,WeekFrac*AMData[iArea].MSData[iMark].Capacity[iWeek]) 
                <= mark[iArea=1:NArea,iMark=1:AMData[iArea].NMStep,k=1:NK] 
                <= WeekFrac*max(0.0,AMData[iArea].MSData[iMark].Capacity[iWeek]), base_name="mark")            # GWh Purchase(+) or sales(-) from market 
      @variable(M,0.0 <= etran[iLine=1:NLine,k=1:NK] <= WeekFrac*LineCap[iLine], base_name="etran")            # GWh Exchange 
      @variable(M,0.0 <= wprod[iArea=1:NArea,k=1:NK] <= CNS.Big,base_name="wprod")                             # GWh/step
      @variable(M,0.0 <= rat[iArea=1:NArea,k=1:NK] <= CNS.Big,base_name="rat")                                 # GWh
      @variable(M,0.0 <= h2dis[iArea=1:NH2Area,k=1:NK] <= H2Data.Areas[iArea].MaxDis, base_name="h2dis")       # GWh/step
      @variable(M,0.0 <= h2chg[iArea=1:NH2Area,k=1:NK] <= H2Data.Areas[iArea].MaxDis, base_name="h2chg")       # GWh/step
      @variable(M,0.0 <= h2res[iSys=1:NH2Area,k=1:NK] <= H2Data.Areas[iSys].MaxRes, base_name="h2res")         # GWh
      @variable(M,h2init[iArea=1:NH2Area],base_name="h2init")                                                  # GWh
      @variable(M,-CNS.AlphaMax <= alpha <= CNS.AlphaMax,base_name="alp")                                      # 10E3 EUR

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
                 +sum(CNS.CNegRes*ner[iSys,k] for iSys=1:NHSys for k=1:NK)
                 +sum(WeekFrac*CNS.CMinRes*mir[iSys,k] for iSys=1:NHSys for k=1:NK)
                 +sum(CNS.CSpi*spi[iSys,k] for iSys=1:NHSys for k=1:NK)
                 +sum(CNS.CRat*rat[iArea,k] for iArea=1:NArea for k=1:NK)
                 +sum(CNS.CFeas*fsl[iSys,c] for iSys=1:NHSys for c=1:NFeasCut[iSys])+alpha)
     
      #Inflow = Phi*Sigma*zinit+Sigma*Eps(t) + MU. First term on LHS. DevCorrMat is the coefficient Phi(NHSys*NHSys).*Sigma(NHSys*1).
      DevCorrMat = IM.InflowSDev[1:NHSys,iWeek].*IM.CorrMat #Dim(NSer*NSer)
      MuCorrVec = zeros(Float64,NHSys,NCut)
      for c = 1:NCut
         MuCorrVec[1:NHSys,c] = IM.CorrMat*CCI[1:NHSys,t,c]
      end

      #Initial reservoir balance [GWh] 
      @constraint(M,resbal0[iSys=1:NHSys,k=[1]],res[iSys,k]+prod[iSys,k]+spi[iSys,k]-ner[iSys,k]-rinit[iSys]-WeekFrac*inf[iSys] == 0.0) 

      #Reservoir balance [GWh]
      @constraint(M,resbal[iSys=1:NHSys,k=2:NK],res[iSys,k]-res[iSys,k-1]+prod[iSys,k]+spi[iSys,k]-ner[iSys,k]-WeekFrac*inf[iSys] == 0.0) 

      #Minimum reservoir limit soft [GWh] 
      #@constraint(M,resminSoft[iSys=1:NHSys,k=1:NK],res[iSys,k]+mir[iSys,k] >= HSys[iSys].MinRes) 
   
      #Reservoir state: V_{t-1}
      @constraint(M,rstate[iSys=1:NHSys],rinit[iSys] == 0.0)

      #Normalized inflow: Z_{t-1} 
      @constraint(M,zstate[iSys=1:NHSys],zinit[iSys] == 0.0)

      #Physical inflow [GWh]
      @constraint(M,inflow[iSys=1:NHSys],inf[iSys] - sum(DevCorrMat[iSys,jSys]*zinit[jSys] for jSys=1:NHSys) <= 0.0)

      #POWER BALANCE [GWh/step]
      @constraint(M,pbal[iArea=1:NArea,k=1:NK],
                  +sum(prod[AreaSys[iArea,iSys],k] for iSys=1:NAreaSys[iArea])+rat[iArea,k]+sum(mark[iArea,iMark,k] for iMark=1:AMData[iArea].NMStep)+wprod[iArea,k] 
                  -sum(etran[MCon[iArea].LIndxOut[iLine],k] for iLine=1:MCon[iArea].NCon)
                  +sum((1.0-LineLoss[MCon[iArea].LIndxIn[iLine]])*etran[MCon[iArea].LIndxIn[iLine],k] for iLine=1:MCon[iArea].NCon)
                  - (LDemandResponse ? dr_tot[iArea,k] : 0) #Include dr_tot in power balance only if LDemandResponse is true
                  - (H2Data.Ind[iArea] > 0 ? h2chg[H2Data.Ind[iArea],k] : 0.0) 
                  + (H2Data.Ind[iArea] > 0 ? h2dis[H2Data.Ind[iArea],k] : 0.0) 
                  ==sum(AMData[iArea].MLData[iLoad].Load[iWeek,k] for iLoad=1:AMData[iArea].NLoad))

      #RAMP UP [GWh/step]
      #@constraint(M,rampup[iSys=1:NHSys,k=2:NK],prod[iSys,k] - prod[iSys,k-1] <= ramp[iSys])

      #RAMP DOWN [GWh/step]
      #@constraint(M,rampdn[iSys=1:NHSys,k=2:NK],prod[iSys,k-1] - prod[iSys,k] <= ramp[iSys])

      #CAPACITY CONSTRAINTS PER AREA [GWh/step]
      #   -Bounds minimum and maximum production when cap=0
      @constraint(M,capup[iSys=1:NHSys,k=1:NK],prod[iSys,k]+cap[iSys] <= WeekFrac*HSys[iSys].MaxProd)
      @constraint(M,capdn[iSys=1:NHSys,k=1:NK],prod[iSys,k]-cap[iSys] >= WeekFrac*HSys[iSys].MinProd[iWeek])

      #H2 STORAGE INIT [GWh]
      @constraint(M,h2storage0[iArea=1:NH2Area,k=[1]],h2res[iArea,k]-(1.0-H2Data.Areas[iArea].CompLoss)*h2chg[iArea,k]+h2dis[iArea,k]-h2init[iArea]== 0.0)

      #H2 STORAGE [GWh]
      @constraint(M,h2storage[iArea=1:NH2Area,k=2:NK],h2res[iArea,k]-(1.0-H2Data.Areas[iArea].CompLoss)*h2chg[iArea,k]+h2dis[iArea,k]-h2res[iArea,k-1] == 0.0)

      #H2 storage state: s_{t-1}
      @constraint(M,h2state[iArea=1:NH2Area],h2init[iArea] == 0.0)

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



      #CAPACITY CONSTRAINTS FOR SYSTEM [GWh/step]
      #@constraint(M,capbal[k=1:NK],sum(cap[iSys] for iArea=1:NArea for iSys=1:NAreaSys[iArea]) >= WeekFrac*CapReq)

      #WIND POWER TARGET [GWh/step]
      @constraint(M,wptarget[iArea=1:NArea,k=1:NK], wprod[iArea,k] <= max(MyWPData[iArea,k],0.0))

      if LEndVal
         if EV.NEndCut > 0
            @constraint(M,endc[c=1:EV.NEndCut],alpha+sum(EV.EndCutCoef[c,iSys]*res[iSys,end] for iSys=1:NHSys) >= EV.EndCutRHS[c])
         elseif EV.NEndSeg > 0
             SegFrac = 1.0/Float64(EV.NEndSeg)
             @variable(M,0.0 <= res_seg[iSys=1:NHSys,iSeg=1:EV.NEndSeg] <= SegFrac*HSys[iSys].MaxRes, base_name="res_seg") 
             @constraint(M,endsegval,alpha+sum(EV.EndSegCoef[iSeg]*res_seg[iSys,iSeg] for iSys=1:NHSys for iSeg=1:EV.NEndSeg) >= 0.0)
             @constraint(M,endsegbal[iSys=1:NHSys],sum(res_seg[iSys,iSeg] for iSeg=1:EV.NEndSeg) - res[iSys,end] == 0.0)
         else
            @constraint(M,endset,alpha == 0)
         end
      else
         @constraint(M,cut[c=1:NCut],alpha-sum(CCR[iSys,t,c]*res[iSys,end] for iSys=1:NHSys)
                     -sum(CCH[iArea,t,c]*h2res[iArea,end] for iArea=1:NH2Area)
                     -sum(MuCorrVec[iSys,c]*zinit[iSys] for iSys=1:NHSys) >= 0.0)
      end
      

      if LFeasSpace
         @constraint(M,feas[iSys=1:NHSys,c=1:NFeasCut[iSys]], FCC[iSys,c,1]*res[iSys,end] + FCC[iSys,c,2]*sum(prod[iSys,k] for k=1:NK)
                     +FCC[iSys,c,3]*ramp[iSys] +FCC[iSys,c,4]*cap[iSys] + FCC[iSys,c,5]*rinit[iSys] + FCC[iSys,c,6]*inf[iSys] -fsl[iSys,c] <= FCC[iSys,c,7])
      end

      return M
   end
end

