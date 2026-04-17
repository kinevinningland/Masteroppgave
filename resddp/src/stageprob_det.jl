#Detailed stage problem

module StageProbDet
   using JuMP 
   using MathOptInterface
   const MOI = MathOptInterface

   function Build(t,iWeek,USMod,AHData,AMData,HSys,MCon,EV,CCR,CCH,CNS,NCut,NHSys,NArea,NLine,LineCap,LineLoss,CTI,LEndVal,LDemandResponse,DR,optimizer)

      M = Model(optimizer)

      NK = CTI.NK
      DT = CTI.DT
      WeekFrac = CTI.WeekFrac

      M3S2MM3 = CNS.M3S2MM3
      MAGEFF2GWH = CNS.MAGEFF2GWH
      MW2GWHWEEK = CNS.MW2GWHWEEK
      InfUB = CNS.InfUB

      NLoadRecStep = DR.NLoadRecStep
      MaxLoadRec = DR.LoadRec[NLoadRecStep]

      @variable(M,0.0 <= res[iArea=1:NHSys,iMod=1:AHData[iArea].NMod,k=1:NK] <= AHData[iArea].MData[iMod].MaxRes, base_name="res")  #Mm3
      @variable(M,0.0 <= disSeg[iArea=1:NHSys,iMod=1:AHData[iArea].NMod,iSeg=1:AHData[iArea].PQData[iMod].NSeg,k=1:NK] <= AHData[iArea].PQData[iMod].DMax[iSeg], base_name="disSeg") #m3s
      @variable(M,0.0 <= dis[iArea=1:NHSys,iMod=1:AHData[iArea].NMod,k=1:NK] <= InfUB, base_name="dis")                             #m3s                                
      @variable(M,0.0 <= spi[iArea=1:NHSys,iMod=1:AHData[iArea].NMod,k=1:NK] <= InfUB, base_name="spi")                             #m3s
      @variable(M,0.0 <= byp[iArea=1:NHSys,iMod=1:AHData[iArea].NMod,k=1:NK] <= InfUB, base_name="byp")                             #m3s
      @variable(M,0.0 <= ghy[iArea=1:NHSys,iMod=1:AHData[iArea].NMod,k=1:NK] <= InfUB, base_name="ghy")                             #GWh/step
      @variable(M,0.0 <= rstate[iArea=1:NHSys] <= InfUB, base_name="rstate")                                                        #GWh/step
      @variable(M,0.0 <= wprod[iArea=1:NArea,k=1:NK] <= CNS.Big,base_name="wprod")                                                  # GWh/step


      #Demand response variables
      if LDemandResponse      
         @variable(M,-InfUB <= dr_tot[iArea=1:NArea, k=1:NK] <= InfUB, base_name="dr_tot")
         @variable(M,0.0 <= dr_up[iArea=1:NArea, k=1:NK] <= InfUB, base_name="dr_up")
         @variable(M,0.0 <= dr_dn[iArea=1:NArea, kk=1:(2*MaxLoadRec+1), k=1:NK] <= InfUB, base_name="dr_dn")
         @variable(M,0.0 <= dr_dir[iArea=1:NArea, k=1:NK] <= 1.0, base_name="dr_dir")
      end

      @variable(M,min(0.0,WeekFrac*AMData[iArea].MSData[iMark].Capacity[iWeek]) 
                <= mark[iArea=1:NArea,iMark=1:AMData[iArea].NMStep,k=1:NK] 
                <= WeekFrac*max(0.0,AMData[iArea].MSData[iMark].Capacity[iWeek]), base_name="mark")                                 # GWh Purchase(+) or sales(-) from market 
      @variable(M,0.0 <= etran[iLine=1:NLine,k=1:NK] <= WeekFrac*LineCap[iLine], base_name="etran")                                 # GWh Exchange 
      @variable(M,0.0 <= rat[iArea=1:NArea,k=1:NK] <= CNS.Big,base_name="rat")                                                      # GWh
      @variable(M,-CNS.AlphaMax <= alpha <= CNS.AlphaMax,base_name="alp")                                                           # 10E3 EUR

      #Min Cost [10E3 EUR]
      @objective(M,MathOptInterface.MIN_SENSE,
                 sum(AMData[iArea].MSData[iMark].Price[iWeek]*mark[iArea,iMark,k] for iArea=1:NArea for iMark=1:AMData[iArea].NMStep for k=1:NK)
                 +sum(CNS.CSpi*spi[iArea,iMod,k] for iArea=1:NHSys for iMod=1:AHData[iArea].NMod for k=1:NK) 
                 +sum(CNS.CByp*byp[iArea,iMod,k] for iArea=1:NHSys for iMod=1:AHData[iArea].NMod for k=1:NK) 
                 +sum(CNS.CRat*rat[iArea,k] for iArea=1:NArea for k=1:NK) +alpha)

      #INITIAL RESERVOIR BALANCE [MM3/DT]
      @constraint(M,resbalReg0[iArea=1:NHSys,iMod=1:AHData[iArea].NMod],res[iArea,iMod,1]
                +DT*M3S2MM3*(dis[iArea,iMod,1]+byp[iArea,iMod,1]+spi[iArea,iMod,1]) 
                -DT*M3S2MM3*(sum(dis[iArea,USMod[iArea].USModA[iMod].Dis[iDisUS],1] for iDisUS=1:USMod[iArea].USModA[iMod].NDis))
                -DT*M3S2MM3*(sum(byp[iArea,USMod[iArea].USModA[iMod].Byp[iBypUS],1] for iBypUS=1:USMod[iArea].USModA[iMod].NByp))
                -DT*M3S2MM3*(sum(spi[iArea,USMod[iArea].USModA[iMod].Spi[iSpiUS],1] for iSpiUS=1:USMod[iArea].USModA[iMod].NSpi))
                == 0.0)

      #RESERVOIR BALANCE [MM3/DT]
      @constraint(M,resbalReg[iArea=1:NHSys,iMod=1:AHData[iArea].NMod,k=2:NK],res[iArea,iMod,k]-res[iArea,iMod,k-1]
                +DT*M3S2MM3*(dis[iArea,iMod,k]+byp[iArea,iMod,k]+spi[iArea,iMod,k]) 
                -DT*M3S2MM3*(sum(dis[iArea,USMod[iArea].USModA[iMod].Dis[iDisUS],k] for iDisUS=1:USMod[iArea].USModA[iMod].NDis))
                -DT*M3S2MM3*(sum(byp[iArea,USMod[iArea].USModA[iMod].Byp[iBypUS],k] for iBypUS=1:USMod[iArea].USModA[iMod].NByp))
                -DT*M3S2MM3*(sum(spi[iArea,USMod[iArea].USModA[iMod].Spi[iSpiUS],k] for iSpiUS=1:USMod[iArea].USModA[iMod].NSpi))
                == 0.0)

      #DEFINE DISCHARGE [m3s]
      @constraint(M,discharge[iArea=1:NHSys,iMod=1:AHData[iArea].NMod,k=1:NK], 
                  dis[iArea,iMod,k]-sum(disSeg[iArea,iMod,iSeg,k] for iSeg=1:AHData[iArea].PQData[iMod].NSeg) == 0.0)

      #HYDROPOWER GENERATION [GWh/Step]
      @constraint(M,prodfunc[iArea=1:NHSys,iMod=1:AHData[iArea].NMod,k=1:NK],ghy[iArea,iMod,k]
                  -WeekFrac*MW2GWHWEEK*sum(AHData[iArea].PQData[iMod].Eff[iSeg]*disSeg[iArea,iMod,iSeg,k] for iSeg=1:AHData[iArea].PQData[iMod].NSeg) == 0.0)

      #Reservoir state: V_{t}
      @constraint(M,endvol[iArea=1:NHSys], rstate[iArea] - sum(AHData[iArea].EffSea[iMod]*MAGEFF2GWH*res[iArea,iMod,NK] for iMod=1:AHData[iArea].NMod) == 0.0)


      #POWER BALANCE [GWh/step]
      @constraint(M,pbalHyd[iArea=1:NHSys,k=1:NK],sum(ghy[iArea,iMod,k] for iMod=1:AHData[iArea].NMod) 
                  +rat[iArea,k]+sum(mark[iArea,iMark,k] for iMark=1:AMData[iArea].NMStep)+wprod[iArea,k]
                  -sum(etran[MCon[iArea].LIndxOut[iLine],k] for iLine=1:MCon[iArea].NCon)
                  +sum((1.0-LineLoss[MCon[iArea].LIndxIn[iLine]])*etran[MCon[iArea].LIndxIn[iLine],k] for iLine=1:MCon[iArea].NCon)
                  - (LDemandResponse ? dr_tot[iArea,k] : 0) #Include dr_tot in power balance only if LDemandResponse is true
                  == (AMData[iArea].NLoad > 0 ? sum(AMData[iArea].MLData[iLoad].Load[iWeek,k] for iLoad=1:AMData[iArea].NLoad) : 0.0))

      #POWER BALANCE [GWh/step]
      @constraint(M,pbalTerm[iArea=(NHSys+1):NArea,k=1:NK],rat[iArea,k]+sum(mark[iArea,iMark,k] for iMark=1:AMData[iArea].NMStep)+wprod[iArea,k]
                  -sum(etran[MCon[iArea].LIndxOut[iLine],k] for iLine=1:MCon[iArea].NCon)
                  +sum((1.0-LineLoss[MCon[iArea].LIndxIn[iLine]])*etran[MCon[iArea].LIndxIn[iLine],k] for iLine=1:MCon[iArea].NCon)
                  - (LDemandResponse ? dr_tot[iArea,k] : 0) #Include dr_tot in power balance only if LDemandResponse is true
                  == (AMData[iArea].NLoad > 0 ? sum(AMData[iArea].MLData[iLoad].Load[iWeek,k] for iLoad=1:AMData[iArea].NLoad) : 0.0))

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
      @constraint(M,wptarget[iArea=1:NArea,k=1:NK], wprod[iArea,k] <= 0.0)

      if LEndVal
         if EV.NEndCut > 0
            @constraint(M,endc[c=1:EV.NEndCut],alpha+sum(EV.EndCutCoef[c,iSys]*rstate[iSys] for iSys=1:NHSys) >= EV.EndCutRHS[c])
         elseif EV.NEndSeg > 0
             SegFrac = 1.0/Float64(EV.NEndSeg)
             @variable(M,0.0 <= res_seg[iSys=1:NHSys,iSeg=1:EV.NEndSeg] <= SegFrac*HSys[iSys].MaxRes, base_name="res_seg") 
             @constraint(M,endsegval,alpha+sum(EV.EndSegCoef[iSeg]*res_seg[iSys,iSeg] for iSys=1:NHSys for iSeg=1:EV.NEndSeg) >= 0.0)
             @constraint(M,endsegbal[iSys=1:NHSys],sum(res_seg[iSys,iSeg] for iSeg=1:EV.NEndSeg) - rstate[iSys] == 0.0)
         else
            @constraint(M,endset,alpha == 0)
         end
      else
         @constraint(M,cut[c=1:NCut],alpha-sum(CCR[iSys,t,c]*rstate[iSys] for iSys=1:NHSys) >= 0.0)
      end

      return M
   end
end
