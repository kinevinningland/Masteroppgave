#Detailed stage problem

module StageProbDet
   using JuMP 
   using MathOptInterface
   const MOI = MathOptInterface

   function Build(t,iWeek,USMod,AHData,AMData,HSys,MCon,EV,CCR,CCH,CNS,NCut,NHSys,NArea,NLine,LineCap,LineLoss,CTI,LEndVal,LDemandResponse,DR,LOperatingReserves,ORData,optimizer)

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

      @variable(M, 0 <= slackDown[z=1:ORData.NZ, k=1:NK], base_name="slackDown")
      @variable(M, 0 <= slackUp[z=1:ORData.NZ, k=1:NK],   base_name="slackUp")

      #Min Cost [10E3 EUR]
      @objective(M,MathOptInterface.MIN_SENSE,
                 sum(AMData[iArea].MSData[iMark].Price[iWeek]*mark[iArea,iMark,k] for iArea=1:NArea for iMark=1:AMData[iArea].NMStep for k=1:NK)
                 +sum(CNS.CSpi*spi[iArea,iMod,k] for iArea=1:NHSys for iMod=1:AHData[iArea].NMod for k=1:NK) 
                 +sum(CNS.CByp*byp[iArea,iMod,k] for iArea=1:NHSys for iMod=1:AHData[iArea].NMod for k=1:NK) 
                 +sum(CNS.CRat*rat[iArea,k] for iArea=1:NArea for k=1:NK) +alpha
                 + (LOperatingReserves ? sum(CNS.CRat*slackUp[z,k]   for z=1:ORData.NZ-1 for k=1:NK) + sum(CNS.CRat*slackDown[z,k] for z=1:ORData.NZ-1 for k=1:NK) : 0.0))#ADDED

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
      
      if LOperatingReserves #ADDED
         NZ = ORData.NZ
         NZa = ORData.NZ_active
         price_zones = ORData.price_zones
         zone_reqs = ORData.zone_reqs
         area_to_zone = ORData.area_to_zone
         areas_in_zone = ORData.areas_in_zone
         hydrosys_to_area = ORData.hydrosys_to_area
         h2_to_area = ORData.h2_to_area
         pos_by_area = ORData.pos_by_area
         neg_by_area = ORData.neg_by_area

         #Opp- og nedreguleringsreserver
         @variable(M, 0 <= cap_zone_up[z=1:NZ, k=1:NK], base_name="cap_zone_up")
         @variable(M, 0 <= cap_zone_down[z=1:NZ, k=1:NK], base_name="cap_zone_down")
         @variable(M, 0 <= cap_hydro_up_mod[iArea=1:NHSys,iMod=1:AHData[iArea].NMod,k=1:NK], base_name="cap_hydro_up_mod")
         @variable(M, 0 <= cap_hydro_down_mod[iArea=1:NHSys,iMod=1:AHData[iArea].NMod,k=1:NK],base_name="cap_hydro_down_mod")
         #@expression(M, cap_hydro_up[iArea=1:NHSys,k=1:NK], sum(cap_hydro_up_mod[iArea,iMod,k] for iMod=1:AHData[iArea].NMod))
         #@expression(M, cap_hydro_down[iArea=1:NHSys,k=1:NK], sum(cap_hydro_down_mod[iArea,iMod,k] for iMod=1:AHData[iArea].NMod))
         @variable(M, 0 <= cap_wind_down[iArea=1:NArea,k=1:NK], base_name="cap_wind_down")
         if ORData.LMarkReserves
            @variable(M, 0 <= cap_mark_up_pos[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK],   base_name="cap_mark_up_pos")
            @variable(M, 0 <= cap_mark_down_pos[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK], base_name="cap_mark_down_pos")
            @variable(M, 0 <= cap_mark_up_neg[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK],   base_name="cap_mark_up_neg")
            @variable(M, 0 <= cap_mark_down_neg[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK], base_name="cap_mark_down_neg")
         end

         #Diverse
         @variable(M, wp_avail[a=1:NArea, k=1:NK] >= 0, base_name="wp_avail")
         @constraint(M, wp_avail_fix[a=1:NArea, k=1:NK],wp_avail[a,k] == 0.0)
         
         
         sharing = true #legge inn i ORData
         if sharing
            sharing_amount = 0.01
            @variable(M, 0 <= sharing_up[z1=1:NZa,z2=1:NZa, k=1:NK], base_name="sharing_up")
            #@variable(M, 0 <= sharing_down[z1=1:NZa,z2=1:NZa, k=1:NK], base_name="sharing_down")
            #constraint som passer på at sonene grenser
            neighboring_zones = Set{Tuple{Int, Int}}()
            for iArea in 1:NArea
               z1 = area_to_zone[iArea]
               for iLine in 1:MCon[iArea].NCon
                  lineIdx = MCon[iArea].LIndxOut[iLine]
                  for a in 1:NArea
                     if lineIdx in MCon[a].LIndxIn
                        z2 = area_to_zone[a]
                        if z1 != z2
                        push!(neighboring_zones, (min(z1,z2), max(z1,z2)))
                        end
                     end
                  end
               end
            end 
            @constraint(M, sharing_balance_up[z1=1:NZa, z2=1:NZa, k=1:NK; z1 != z2 && (min(z1,z2),max(z1,z2)) in neighboring_zones], sharing_up[z1,z2,k] <= sharing_amount * (zone_reqs[z2].RI_up * 3 + zone_reqs[z2].NI_up * sum(wp_avail[a,k] for a in areas_in_zone[z2]; init=0.0) + 0.03 * sum(AMData[iArea].MLData[iLoad].Load[iWeek,k] for iArea in areas_in_zone[z2] for iLoad in 1:AMData[iArea].NLoad; init=0.0)))
   
         end
         


         #Oppreguleringskrav
         @constraint(M, reserve_req_up[z=1:NZ-1, k=1:NK],
         cap_zone_up[z,k] + slackUp[z,k] >= (zone_reqs[z].RI_up * 3 #Fikse at 3 er DT
            + zone_reqs[z].NI_up * sum(wp_avail[a,k] for a in areas_in_zone[z]; init=0.0) + 0.03 * sum(AMData[iArea].MLData[iLoad].Load[iWeek,k]
                 for iArea in areas_in_zone[z]
                 for iLoad in 1:AMData[iArea].NLoad; init=0.0))
            )


         #Nedreguleringskrav
         @constraint(M, reserve_req_down[z=1:NZ-1, k=1:NK],
         cap_zone_down[z,k] + slackDown[z,k] >= (zone_reqs[z].RI_down * 3 #Fikse at 3 er DT
            + zone_reqs[z].NI_down * sum(wp_avail[a,k] for a in areas_in_zone[z]; init=0.0) + 0.03 * sum(AMData[iArea].MLData[iLoad].Load[iWeek,k]
                 for iArea in areas_in_zone[z]
                 for iLoad in 1:AMData[iArea].NLoad; init=0.0)))
         

         #Sammenhengen mellom zonesum av kapasiteter og sum av individuelle kapasiteter
         @constraint(M, reserve_split_down[z=1:NZ, k=1:NK],
         cap_zone_down[z,k] ==
            #sum(cap_hydro_down[iSys, k] for iSys in 1:NHSys if (hydrosys_to_area[iSys] in areas_in_zone[z]); init=0.0) 
            sum(sum(cap_hydro_down_mod[iSys, iMod, k] for iMod in 1:AHData[iSys].NMod) for iSys in 1:NHSys if (hydrosys_to_area[iSys] in areas_in_zone[z]); init=0.0)
            + sum(cap_wind_down[a,k] for a in areas_in_zone[z]; init=0.0) 
            + (ORData.LMarkReserves ? sum(cap_mark_down_pos[a, iMark, k] for a in areas_in_zone[z] for iMark in get(ORData.pos_by_area, a, Set{Int}()); init=0.0) : 0.0)
            + (ORData.LMarkReserves ? sum(cap_mark_down_neg[a, iMark, k] for a in areas_in_zone[z] for iMark in get(ORData.neg_by_area, a, Set{Int}()); init=0.0) : 0.0)
         )

         @constraint(M, reserve_split_up[z=1:NZ, k=1:NK],
         cap_zone_up[z,k] ==
            #sum(cap_hydro_up[iSys, k] for iSys in 1:NHSys if (hydrosys_to_area[iSys] in areas_in_zone[z]); init=0.0) 
            sum(sum(cap_hydro_up_mod[iSys, iMod, k] for iMod in 1:AHData[iSys].NMod) for iSys in 1:NHSys if (hydrosys_to_area[iSys] in areas_in_zone[z]); init=0.0)
            #+ (ORData.LMarkReserves ? sum(cap_mark_up[a, iMark, k] for a in areas_in_zone[z] for iMark in get(ORData.pos_by_area, a, Set{Int}()); init=0.0) : 0.0) #Fikse denne
            + (ORData.LMarkReserves ? sum(cap_mark_up_pos[a, iMark, k] for a in areas_in_zone[z] for iMark in get(ORData.pos_by_area, a, Set{Int}()); init=0.0) : 0.0)
            + (ORData.LMarkReserves ? sum(cap_mark_up_neg[a, iMark, k] for a in areas_in_zone[z] for iMark in get(ORData.neg_by_area, a, Set{Int}()); init=0.0) : 0.0)
            + (sharing && z <= NZa ? sum(sharing_up[z1, z, k] for z1 in 1:NZa if z1 != z && (min(z1,z),max(z1,z)) in neighboring_zones; init=0.0) : 0.0)
            - (sharing && z <= NZa ? sum(sharing_up[z, z2, k] for z2 in 1:NZa if z2 != z && (min(z,z2),max(z,z2)) in neighboring_zones; init=0.0) : 0.0)
         )


         #koble hver teknologis cap-variabel til dens egne fysiske grenser
         @constraint(M, hydro_up_mod[iArea=1:NHSys, iMod=1:AHData[iArea].NMod, k=1:NK], cap_hydro_up_mod[iArea,iMod,k] <= WeekFrac * MW2GWHWEEK * sum(AHData[iArea].PQData[iMod].Eff[iSeg] * (AHData[iArea].PQData[iMod].DMax[iSeg] - disSeg[iArea,iMod,iSeg,k]) for iSeg=1:AHData[iArea].PQData[iMod].NSeg))
         @constraint(M, hydro_down_mod[iArea=1:NHSys, iMod=1:AHData[iArea].NMod, k=1:NK], cap_hydro_down_mod[iArea,iMod,k] <= ghy[iArea,iMod,k])
         @constraint(M, hydro_res_up_mod[iArea=1:NHSys, iMod=1:AHData[iArea].NMod, k=1:NK], cap_hydro_up_mod[iArea,iMod,k] <= AHData[iArea].EffSea[iMod] * MAGEFF2GWH * res[iArea,iMod,k])
         @constraint(M, hydro_down_mod_res[iArea=1:NHSys, iMod=1:AHData[iArea].NMod, k=1:NK], cap_hydro_down_mod[iArea,iMod,k] <= AHData[iArea].EffSea[iMod] * MAGEFF2GWH * (AHData[iArea].MData[iMod].MaxRes - res[iArea,iMod,k]))
         
         @constraint(M, wind_dn[iArea=1:NArea,k=1:NK], wprod[iArea,k] >= cap_wind_down[iArea,k]) #OK
         @constraint(M, wind_dn2[iArea=1:NArea,k=1:NK], cap_wind_down[iArea,k] <= 0.0) #OK

         if ORData.LMarkReserves
            #@constraint(M, mark_up_pos[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK; iMark in get(pos_by_area, a, Set{Int}())], mark[a,iMark,k] + cap_mark_up_pos[a,iMark,k] <= WeekFrac * max(0.0, AMData[a].MSData[iMark].Capacity[iWeek])) #OK
            #@constraint(M, mark_dn_neg[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK; iMark in get(neg_by_area, a, Set{Int}())], mark[a,iMark,k] - cap_mark_down_neg[a,iMark,k] >= WeekFrac * min(0.0, AMData[a].MSData[iMark].Capacity[iWeek]))#OK
            #@constraint(M, mark_up_neg[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK; iMark in get(neg_by_area, a, Set{Int}())], -mark[a,iMark,k] >= cap_mark_up_neg[a,iMark,k]) #OK
            #@constraint(M, mark_dn_pos[a=1:NArea, iMark=1:AMData[a].NMStep, k=1:NK; iMark in get(pos_by_area, a, Set{Int}())], mark[a,iMark,k] >= cap_mark_down_pos[a,iMark,k]) #OK     
         end
         
         #Går ann å sette cap til null og ta bort dens constraints som allerede ligger i modellen
         #Ta else (hvis ikke OR: sette på dens constraints)
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

function ReadOperatingReserves(NArea, NHSys, NAreaSys, AreaSys, H2Data, AMData,AreaName,LOperatingReserves)
    #Return dummy object if OR is not included
    if !LOperatingReserves
        return OperatingReserves(0,0,String[],ReserveZoneReq[],Int[],Vector{Vector{Int}}(),Int[],Int[],false,false,Dict{Int, Set{Int}}(),Dict{Int, Set{Int}}())
    end

    LH2Reserves = false
    LMarkReserves = false

    price_zones = ["NO1", "NO2", "NO3", "NO4","Others"]
    NZ = length(price_zones)
    zone_reqs = [
        ReserveZoneReq("NO1", 0.344, 0.172, 0.58, 0.56),
        ReserveZoneReq("NO2", 1.4, 1.4, 0.37, 0.37),
        ReserveZoneReq("NO3", 0.29, 0.145, 0.33, 0.38),
        ReserveZoneReq("NO4", 0.35, 0.175, 0.31, 0.33),
        #ReserveZoneReq("NO5", 1.4, 1.4, 0.37, 0.37),
        ReserveZoneReq("Others", 0.0, 0.0, 0.0, 0.0)
    ]

    area_to_zone = fill(findfirst(==("Others"), price_zones),NArea) #Areas not explicitly listed default to "Others"
    area_to_zone[33] = findfirst(==("NO3"), price_zones) #OK
    area_to_zone[34] = findfirst(==("NO2"), price_zones) #OK Var NO5
    area_to_zone[35] = findfirst(==("NO2"), price_zones) #OK
    area_to_zone[36] = findfirst(==("NO2"), price_zones) #OK
    area_to_zone[39] = findfirst(==("Others"), price_zones) #OK
    area_to_zone[56] = findfirst(==("Others"), price_zones) # DK1 eller DK2? Hydrogen danmark
    area_to_zone[66] = findfirst(==("NO3"), price_zones) # H2-M?
    area_to_zone[67] = findfirst(==("NO2"), price_zones) #H2-S?
    area_to_zone[68] = findfirst(==("NO4"), price_zones) #H2-N?
    area_to_zone[73] = findfirst(==("NO4"), price_zones) #OK
    area_to_zone[74] = findfirst(==("NO4"), price_zones) #OK
    area_to_zone[75] = findfirst(==("NO4"), price_zones) #OK
    area_to_zone[57] = findfirst(==("Others"), price_zones) #OK
    area_to_zone[62] = findfirst(==("Others"), price_zones) #Sweden-H2?
       

    area_to_zone[1] = findfirst(==("NO1"), price_zones) #OK
    area_to_zone[2] = findfirst(==("NO1"), price_zones) #spørre om NO1 eller NO2, SørØst
    area_to_zone[3] = findfirst(==("NO1"), price_zones) #spørre, NO1,NO3,NO5? Hallingdal
    area_to_zone[4] = findfirst(==("NO2"), price_zones) #spørre NO1,NO2? Telemark

    area_to_zone[5] = findfirst(==("NO2"), price_zones) #OK
    area_to_zone[6] = findfirst(==("NO2"), price_zones) #OK

    area_to_zone[8] = findfirst(==("NO3"), price_zones) #OK
    area_to_zone[9] = findfirst(==("NO4"), price_zones) #OK

    area_to_zone[10] = findfirst(==("NO4"), price_zones) #OK
    area_to_zone[11] = findfirst(==("NO4"), price_zones) #OK

    area_to_zone[7] = findfirst(==("NO2"), price_zones) #OK Var NO5

    area_to_zone[12] = findfirst(==("Others"), price_zones) #OK
    area_to_zone[13] = findfirst(==("Others"), price_zones) #OK
    area_to_zone[37] = findfirst(==("Others"), price_zones) #OK

    area_to_zone[14] = findfirst(==("Others"), price_zones) #OK
    area_to_zone[15] = findfirst(==("Others"), price_zones) #OK

    area_to_zone[16] = findfirst(==("Others"), price_zones) #OK
    area_to_zone[38] = findfirst(==("Others"), price_zones) #OK

    area_to_zone[17] = findfirst(==("Others"), price_zones) #OK

    area_to_zone[20] = findfirst(==("Others"), price_zones) #OK
    area_to_zone[42] = findfirst(==("Others"), price_zones) #OK

    area_to_zone[19] = findfirst(==("Others"), price_zones) #OK
    area_to_zone[41] = findfirst(==("Others"), price_zones) #OK

    area_to_zone[18] = findfirst(==("Others"), price_zones) #OK
    area_to_zone[40] = findfirst(==("Others"), price_zones) #OK

    areas_in_zone = [Int[] for _ in 1:NZ]
    for a in 1:NArea
        push!(areas_in_zone[area_to_zone[a]], a)
    end

    hydrosys_to_area = fill(0, NHSys)
    for a in 1:NArea
        for j in 1:NAreaSys[a]
            hydrosys_to_area[AreaSys[a,j]] = a
        end
    end
    
    h2_to_area = Int[]
    if !isnothing(H2Data) && hasproperty(H2Data, :NArea) && hasproperty(H2Data, :Ind) && H2Data.NArea > 0
        h2_to_area = fill(0, H2Data.NArea)  # h2_index -> model area
        for a in 1:NArea
            iH2 = H2Data.Ind[a]             # 0 hvis ikke H2 i området
            if iH2 > 0
                h2_to_area[iH2] = a
            end
        end
    end

    #Mark reserves
    lc(s) = lowercase(String(s))

    excluded(navn::String) = begin
        n = lc(navn)
        occursin("nucl", n) || occursin("nuclear", n) ||
        occursin("el-import", n) || occursin("el-export", n) ||
        occursin("h2-import", n) || occursin("waste", n) ||
        occursin("a/s union", n)
    end

    realistic_pos(navn::String) = !excluded(navn) && (
        occursin("pa. kjop dellast", lc(navn)) ||
        (occursin("hydro", lc(navn)) && !occursin("ps_hydro_con", lc(navn)))
    )

    realistic_neg(navn::String) = !excluded(navn) && (
        occursin("ps_hydro_con", lc(navn)) ||
        occursin("kraftintensiv", lc(navn)) ||
        occursin("kjelkraft", lc(navn)) ||
        occursin("pa. salg dellast", lc(navn))
    )

    pos_by_area = Dict{Int, Set{Int}}()
    neg_by_area = Dict{Int, Set{Int}}()

    for a in 1:NArea
        for iMark in 1:AMData[a].NMStep
            navn   = AMData[a].MSData[iMark].Name
            maxcap = maximum(AMData[a].MSData[iMark].Capacity)
            mincap = minimum(AMData[a].MSData[iMark].Capacity)

            if maxcap > 0 && realistic_pos(navn)
                push!(get!(pos_by_area, a, Set{Int}()), iMark)
            end
            if mincap < 0 && realistic_neg(navn)
                push!(get!(neg_by_area, a, Set{Int}()), iMark)
            end
        end
    end
    #=
    tidsserie_path = joinpath(dataset, "TidsserieData.h5")

    clean_str(x) = strip(replace(String(x), '\0' => ' ')) |>
                   s -> replace(s, r"\s+" => " ") |>
                   strip |>
                   s -> replace(s, r",\s*$" => "")

    lc(s) = lowercase(String(s))

    function get_attr_str(obj, key; default="")
        a = attrs(obj)
        haskey(a, key) || return default
        v = a[key]
        v = (v isa AbstractVector && length(v) > 0) ? v[1] : v
        return clean_str(v)
    end

    function get_attr_int(obj, key)
        s = get_attr_str(obj, key; default="")
        s2 = replace(s, r"[^0-9\\-]" => "")
        isempty(s2) ? missing : try parse(Int, s2) catch; missing end
    end

    excluded(navn::String) = begin
        n = lc(navn)
        occursin("nucl", n) || occursin("nuclear", n) ||
        occursin("el-import", n) || occursin("el-export", n) ||
        occursin("h2-import", n) || occursin("waste", n) ||
        occursin("a/s union", n)
    end

    realistic_pos(navn::String) = !excluded(navn) && (
        occursin("pa. kjop dellast", lc(navn)) ||
        (occursin("hydro", lc(navn)) && !occursin("ps_hydro_con", lc(navn)))
    )

    realistic_neg(navn::String) = !excluded(navn) && (
        occursin("ps_hydro_con", lc(navn)) ||
        occursin("kraftintensiv", lc(navn)) ||
        occursin("kjelkraft", lc(navn)) ||
        occursin("pa. salg dellast", lc(navn))
    )

    # map fra områdenavn -> områdeindeks
    area_index = Dict(strip(name) => i for (i, name) in enumerate(AreaName))

    pos_by_area = Dict{Int, Set{Int}}()
    neg_by_area = Dict{Int, Set{Int}}()

    h5open(tidsserie_path, "r") do f
        for area in keys(f)
            area_name = strip(String(area))
            haskey(area_index, area_name) || continue
            a = area_index[area_name]

            g_area = f[area]
            haskey(g_area, "TRINN") || continue
            g_trinn = g_area["TRINN"]

            for trinn_key in keys(g_trinn)
                startswith(String(trinn_key), "Trinn_") || continue
                g_step = g_trinn[trinn_key]
                haskey(g_step, "Mengde") || continue

                navn = get_attr_str(g_step, "Navn"; default="")
                isempty(navn) && continue
                idv = get_attr_int(g_step, "Id")
                idv === missing && continue

                m = read(g_step["Mengde"])
                mmin = minimum(m)
                mmax = maximum(m)

                if mmax > 0 && realistic_pos(navn)
                    push!(get!(pos_by_area, a, Set{Int}()), idv)
                end
                if mmin < 0 && realistic_neg(navn)
                    push!(get!(neg_by_area, a, Set{Int}()), idv)
                end
            end
        end
    end
    =#


    println("Read ORData.csv")

    return OperatingReserves(NZ,NZ-1,price_zones,zone_reqs,area_to_zone,areas_in_zone,hydrosys_to_area,h2_to_area,LH2Reserves,LMarkReserves,pos_by_area,neg_by_area)
end