#Detailed stage problem

module StageProbDet
   using JuMP 
   using MathOptInterface
   const MOI = MathOptInterface

   function Build(t,iWeek,USMod,AHData,AMData,HSys,MCon,EV,CCR,CCH,CNS,NCut,NHSys,NArea,NLine,LineCap,LineLoss,CTI,LEndVal,LDemandResponse,DR,H2Data,LOperatingReserves,ORData,optimizer) #H2Data Added

      M = Model(optimizer)

      NK = CTI.NK
      DT = CTI.DT
      WeekFrac = CTI.WeekFrac
      NH2Area = H2Data.NArea #ADDED

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
      @variable(M,0.0 <= h2dis[iArea=1:NH2Area,k=1:NK] <= H2Data.Areas[iArea].MaxDis, base_name="h2dis")       # GWh/step ADDED
      @variable(M,0.0 <= h2chg[iArea=1:NH2Area,k=1:NK] <= H2Data.Areas[iArea].MaxDis, base_name="h2chg")       # GWh/step ADDED
      @variable(M,0.0 <= h2res[iSys=1:NH2Area,k=1:NK] <= H2Data.Areas[iSys].MaxRes, base_name="h2res")         # GWh ADDED
      @variable(M,h2init[iArea=1:NH2Area],base_name="h2init")                                                  # GWh ADDED

      println("StageProbDet: NHSys=$NHSys, NArea=$NArea")
      for iSys in 1:NHSys
         println("  HydroSys $iSys: MaxRes=$(HSys[iSys].MaxRes), MaxProd=$(HSys[iSys].MaxProd)")
      end

      for iArea in 1:NH2Area
         println("  H2Area $iArea: MaxRes=$(H2Data.Areas[iArea].MaxRes), MaxDis=$(H2Data.Areas[iArea].MaxDis)")
      end

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

      @variable(M, 0 <= slackDown[z=1:ORData.NZ, k=1:NK], base_name="slackDown")#added
      @variable(M, 0 <= slackUp[z=1:ORData.NZ, k=1:NK],   base_name="slackUp")#added

      #Min Cost [10E3 EUR]
      @objective(M,MathOptInterface.MIN_SENSE,
                 sum(AMData[iArea].MSData[iMark].Price[iWeek]*mark[iArea,iMark,k] for iArea=1:NArea for iMark=1:AMData[iArea].NMStep for k=1:NK)
                 +sum(CNS.CSpi*spi[iArea,iMod,k] for iArea=1:NHSys for iMod=1:AHData[iArea].NMod for k=1:NK) 
                 +sum(CNS.CByp*byp[iArea,iMod,k] for iArea=1:NHSys for iMod=1:AHData[iArea].NMod for k=1:NK) 
                 +sum(CNS.CRat*rat[iArea,k] for iArea=1:NArea for k=1:NK) +alpha
                 + (LOperatingReserves ? sum(CNS.CRat*slackUp[z,k] for z=1:ORData.NZ for k=1:NK) + sum(CNS.CRat*slackDown[z,k] for z=1:ORData.NZ for k=1:NK) : 0.0))#ADDED

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
                  - (H2Data.Ind[iArea] > 0 ? h2chg[H2Data.Ind[iArea],k] : 0.0) #ADDED
                  + (H2Data.Ind[iArea] > 0 ? h2dis[H2Data.Ind[iArea],k] : 0.0) #ADDED
                  - (LDemandResponse ? dr_tot[iArea,k] : 0) #Include dr_tot in power balance only if LDemandResponse is true
                  == (AMData[iArea].NLoad > 0 ? sum(AMData[iArea].MLData[iLoad].Load[iWeek,k] for iLoad=1:AMData[iArea].NLoad) : 0.0))

      #POWER BALANCE [GWh/step]
      @constraint(M,pbalTerm[iArea=(NHSys+1):NArea,k=1:NK],rat[iArea,k]+sum(mark[iArea,iMark,k] for iMark=1:AMData[iArea].NMStep)+wprod[iArea,k]
                  -sum(etran[MCon[iArea].LIndxOut[iLine],k] for iLine=1:MCon[iArea].NCon)
                  +sum((1.0-LineLoss[MCon[iArea].LIndxIn[iLine]])*etran[MCon[iArea].LIndxIn[iLine],k] for iLine=1:MCon[iArea].NCon)
                  - (H2Data.Ind[iArea] > 0 ? h2chg[H2Data.Ind[iArea],k] : 0.0) #ADDED
                  + (H2Data.Ind[iArea] > 0 ? h2dis[H2Data.Ind[iArea],k] : 0.0) #ADDED
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
      
      if LOperatingReserves #ADDED
         NZ = ORData.NZ #antall prissoner
         zone_reqs = ORData.zone_reqs #reservemengder
         areas_in_zone = ORData.areas_in_zone #map over områdene i hver sone
         hydrosys_to_area = ORData.hydrosys_to_area #map over 
         a = ORData.a #empirisk konstant fra Entso-E
         b = ORData.b #empirisk konstant fra Entso-E
         pos_by_area = ORData.pos_by_area
         neg_by_area = ORData.neg_by_area

         
         #Opp- og nedreguleringsreserver
         @variable(M, 0 <= cap_zone_up[z=1:NZ, k=1:NK], base_name="cap_zone_up")
         @variable(M, 0 <= cap_zone_down[z=1:NZ, k=1:NK], base_name="cap_zone_down")
         @variable(M, 0 <= cap_hydro_up_mod[iArea=1:NHSys,iMod=1:AHData[iArea].NMod,k=1:NK], base_name="cap_hydro_up_mod")
         @variable(M, 0 <= cap_hydro_down_mod[iArea=1:NHSys,iMod=1:AHData[iArea].NMod,k=1:NK],base_name="cap_hydro_down_mod")
         @variable(M, 0 <= cap_wind_down[iArea=1:NArea,k=1:NK], base_name="cap_wind_down")

         if ORData.LMarkReserves
            @variable(M, 0 <= cap_mark_up_pos[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK],   base_name="cap_mark_up_pos")
            @variable(M, 0 <= cap_mark_down_pos[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK], base_name="cap_mark_down_pos")
            @variable(M, 0 <= cap_mark_up_neg[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK],   base_name="cap_mark_up_neg")
            @variable(M, 0 <= cap_mark_down_neg[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK], base_name="cap_mark_down_neg")
            @constraint(M, mark_up_pos[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK; iMark in get(pos_by_area, a, Set{Int}())], mark[a,iMark,k] + cap_mark_up_pos[a,iMark,k] <= WeekFrac * max(0.0, AMData[a].MSData[iMark].Capacity[iWeek])) #OK
            @constraint(M, mark_dn_neg[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK; iMark in get(neg_by_area, a, Set{Int}())], mark[a,iMark,k] - cap_mark_down_neg[a,iMark,k] >= WeekFrac * min(0.0, AMData[a].MSData[iMark].Capacity[iWeek]))#OK
            @constraint(M, mark_up_neg[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK; iMark in get(neg_by_area, a, Set{Int}())], -mark[a,iMark,k] >= cap_mark_up_neg[a,iMark,k]) #OK
            @constraint(M, mark_dn_pos[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK; iMark in get(pos_by_area, a, Set{Int}())], mark[a,iMark,k] >= cap_mark_down_pos[a,iMark,k]) #OK     
         end

         #Diverse
         @variable(M, wp_avail[a=1:NArea, k=1:NK] >= 0, base_name="wp_avail")
         @constraint(M, wp_avail_fix[a=1:NArea, k=1:NK],wp_avail[a,k] == 0.0) #RHS endres i simulate


         #Opp- og nedreguleringskrav
         @expression(M, cap_up_amount[z = 1:NZ, k=1:NK], 
            zone_reqs[z].RI_up * CTI.DT                                                                                                   #RI ledd
            + zone_reqs[z].NI_up * sum(wp_avail[a,k] for a in areas_in_zone[z] if !(a in zone_reqs[z].owp_areas_in_zone); init=0.0)       #NI onshore ledd
            + zone_reqs[z].NI_up_OWP * sum(wp_avail[a,k] for a in zone_reqs[z].owp_areas_in_zone; init=0.0)                               #NI offshore ledd
            + (sqrt(a*zone_reqs[z].MaxLoad/(CNS.MW2GW*CTI.DT)+b^2)-b)*CNS.MW2GW*CTI.DT)                                                   #Last ledd
         
         @expression(M, cap_down_amount[z = 1:NZ, k=1:NK], 
            zone_reqs[z].RI_down * CTI.DT                                                                                                 #RI ledd
            + zone_reqs[z].NI_down * sum(wp_avail[a,k] for a in areas_in_zone[z] if !(a in zone_reqs[z].owp_areas_in_zone); init=0.0)     #NI onshore ledd
            + zone_reqs[z].NI_down_OWP * sum(wp_avail[a,k] for a in zone_reqs[z].owp_areas_in_zone; init=0.0)                             #NI offshore ledd
            + (sqrt(a*zone_reqs[z].MaxLoad/(CNS.MW2GW*CTI.DT)+b^2)-b)*CNS.MW2GW*CTI.DT)                                                   #Last ledd
      
         @constraint(M, reserve_req_up[z=1:NZ, k=1:NK], cap_zone_up[z,k] + slackUp[z,k] >= cap_up_amount[z,k])
         @constraint(M, reserve_req_down[z=1:NZ, k=1:NK], cap_zone_down[z,k] + slackDown[z,k] >= cap_down_amount[z,k])


         #Sammenhengen mellom total sonesum og bidrag fra kraftverkene i sonen
         @constraint(M, reserve_split_down[z=1:NZ, k=1:NK], cap_zone_down[z,k] ==
            sum(sum(cap_hydro_down_mod[iSys, iMod, k] for iMod in 1:AHData[iSys].NMod) for iSys in 1:NHSys if (hydrosys_to_area[iSys] in areas_in_zone[z]); init=0.0)
            + sum(cap_wind_down[a,k] for a in areas_in_zone[z]; init=0.0) 
            + (ORData.LMarkReserves ? sum(cap_mark_down_pos[a, iMark, k] for a in areas_in_zone[z] for iMark in get(ORData.pos_by_area, a, Set{Int}()); init=0.0) : 0.0)
            + (ORData.LMarkReserves ? sum(cap_mark_down_neg[a, iMark, k] for a in areas_in_zone[z] for iMark in get(ORData.neg_by_area, a, Set{Int}()); init=0.0) : 0.0)
         
         )

         @constraint(M, reserve_split_up[z=1:NZ, k=1:NK], cap_zone_up[z,k] ==
            sum(sum(cap_hydro_up_mod[iSys, iMod, k] for iMod in 1:AHData[iSys].NMod) for iSys in 1:NHSys if (hydrosys_to_area[iSys] in areas_in_zone[z]); init=0.0)
            + (ORData.LMarkReserves ? sum(cap_mark_up_pos[a, iMark, k] for a in areas_in_zone[z] for iMark in get(ORData.pos_by_area, a, Set{Int}()); init=0.0) : 0.0)
            + (ORData.LMarkReserves ? sum(cap_mark_up_neg[a, iMark, k] for a in areas_in_zone[z] for iMark in get(ORData.neg_by_area, a, Set{Int}()); init=0.0) : 0.0)
         )

         #koble hver teknologis cap-variabel til dens egne fysiske grenser
         @constraint(M, hydro_up_mod[iArea=1:NHSys, iMod=1:AHData[iArea].NMod, k=1:NK], cap_hydro_up_mod[iArea,iMod,k] <= WeekFrac * MW2GWHWEEK * sum(AHData[iArea].PQData[iMod].Eff[iSeg] * (AHData[iArea].PQData[iMod].DMax[iSeg] - disSeg[iArea,iMod,iSeg,k]) for iSeg=1:AHData[iArea].PQData[iMod].NSeg))
         @constraint(M, hydro_down_mod[iArea=1:NHSys, iMod=1:AHData[iArea].NMod, k=1:NK], cap_hydro_down_mod[iArea,iMod,k] <= ghy[iArea,iMod,k])
         @constraint(M, hydro_res_up_mod[iArea=1:NHSys, iMod=1:AHData[iArea].NMod, k=1:NK], cap_hydro_up_mod[iArea,iMod,k] <= AHData[iArea].EffSea[iMod] * MAGEFF2GWH * res[iArea,iMod,k]) #skalere neD? 
         @constraint(M, wind_dn[iArea=1:NArea,k=1:NK], wprod[iArea,k] >= cap_wind_down[iArea,k]) 

      end   

      #H2 STORAGE INIT [GWh]
      @constraint(M,h2storage0[iArea=1:NH2Area,k=[1]],h2res[iArea,k]-(1.0-H2Data.Areas[iArea].CompLoss)*h2chg[iArea,k]+h2dis[iArea,k]-h2init[iArea]== 0.0) # ADDED
      #H2 STORAGE [GWh]
      @constraint(M,h2storage[iArea=1:NH2Area,k=2:NK],h2res[iArea,k]-(1.0-H2Data.Areas[iArea].CompLoss)*h2chg[iArea,k]+h2dis[iArea,k]-h2res[iArea,k-1] == 0.0) # ADDED
      #H2 storage state: s_{t-1}
      @constraint(M,h2state[iArea=1:NH2Area],h2init[iArea] == 0.0) # ADDED
            

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
         #@constraint(M,cut[c=1:NCut],alpha-sum(CCR[iSys,t,c]*rstate[iSys] for iSys=1:NHSys) >= 0.0)
         @constraint(M,cut[c=1:NCut],alpha- sum(CCR[iSys,t,c]*rstate[iSys] for iSys=1:NHSys) - sum(CCH[iArea,t,c]*h2res[iArea,end] for iArea=1:NH2Area)>= 0.0) #h2res Added   
      end

      return M
   end
end
