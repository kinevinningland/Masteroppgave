module HydProb
   using JuMP 
   using MathOptInterface
   const MOI = MathOptInterface

   #Treat as function f(y,z), where y = (v,e,b,c) and z = (I,V0) 
   function Build(iWeek,iSys,iArea,MyHSys,NMod,MyAHData,MyUSMod,DMData,NK,InfUB,DT,CNS,M3S2MM3,GWH2MAGEFF,MAGEFF2GWH,MW2GWHWEEK,WeekFrac,MyRegFrac,StartRamp,EndRamp,optimizer)
      M = Model(optimizer)
      
      @variable(M,DMData.MaMin[iArea,MyHSys.ModNo[iMod],iWeek] <= res[iMod=1:NMod,k=1:NK] <= MyAHData.MData[MyHSys.ModNo[iMod]].MaxRes, base_name="res")                      #Mm3
      @variable(M,0.0 <= res0d[iMod=1:NMod] <= MyAHData.MData[MyHSys.ModNo[iMod]].MaxRes, base_name="res0d")                                                                  #Mm3
      @variable(M,0.0 <= disSeg[iMod=1:NMod,iSeg=1:MyAHData.PQData[MyHSys.ModNo[iMod]].NSeg,k=1:NK] <= MyAHData.PQData[MyHSys.ModNo[iMod]].DMax[iSeg], base_name="disSeg")    #m3s
      @variable(M,0.0 <= dis[iMod=1:NMod,k=1:NK] <= DMData.QMax[iArea,MyHSys.ModNo[iMod],iWeek], base_name="dis")                     #m3s                                
      #@variable(M,DMData.QMin[iArea,MyHSys.ModNo[iMod],iWeek] <= dis[iMod=1:NMod,k=1:NK] <= DMData.QMax[iArea,MyHSys.ModNo[iMod],iWeek], base_name="dis")                     #m3s                                
      @variable(M,0.0 <= spi[iMod=1:NMod,k=1:NK] <= InfUB, base_name="spi")                                                                                                   #m3s
      #@variable(M,DMData.QfoMin[iArea,MyHSys.ModNo[iMod],iWeek] <= byp[iMod=1:NMod,k=1:NK] <= InfUB, base_name="byp")                                                        #m3s
      @variable(M,0.0 <= byp[iMod=1:NMod,k=1:NK] <= InfUB, base_name="byp")                                                                                                   #m3s
      @variable(M,0.0 <= rel[iMod=1:NMod,k=1:NK] <= InfUB, base_name="rel")                                                                                                   #m3s
      @variable(M,0.0 <= ghy[iMod=1:NMod,k=1:NK] <= InfUB, base_name="ghy")                                                                                                   #GWh/step


      #Feasibility variables
      @variable(M,0.0 <= resPen <= InfUB, base_name="resPen")                                                                 #GWh
      @variable(M,0.0 <= enPen <= InfUB, base_name="enPen")                                                                   #GWh
      @variable(M,0.0 <= rampPen <= InfUB, base_name="rampPen")                                                               #GWh/step/stage
      @variable(M,0.0 <= capPen <= InfUB, base_name="capPen")                                                                 #GWh/step/stage
      
      @variable(M,0.0 <= ramp <= InfUB, base_name="ramp")                                                                     #GWh/step/stage
      @variable(M,0.0 <= resEnd <= InfUB, base_name="resEnd")                                                                 #GWh
      @variable(M,0.0 <= resCap <= InfUB, base_name="resCap")                                                                 #GWh/step/stage

      #State variables
      @variable(M,0.0 <= res0 <= MyHSys.MaxRes, base_name="res0")                                                             #GWh
      @variable(M,0.0 <= inf <= InfUB, base_name="inf")                                                                       #GWh/stage

      #OBJECTIVE [10E3 EUR]
      @objective(M,MathOptInterface.MIN_SENSE, CNS.CResPen*resPen+CNS.CEnPen*enPen+CNS.CRampPen*rampPen+CNS.CCapPen*capPen)

      #DEFINE DISCHARGE [m3s]
      @constraint(M,discharge[iMod=1:NMod,k=1:NK], 
                  dis[iMod,k]-sum(disSeg[iMod,iSeg,k] for iSeg=1:MyAHData.PQData[MyHSys.ModNo[iMod]].NSeg) >= 0.0)

      #RAMPING [m3s]
      #@constraint(M,ramping[iMod=1:NMod,k=2:NK], 
                  #-MyAHData.RampFrac[iMod]*MyAHData.PQData[iMod].DMax[1] <= dis[iMod,k]-dis[iMod,k-1] <= MyAHData.RampFrac[iMod]*MyAHData.PQData[iMod].DMax[1])
      
      
      #INITIAL RESERVOIR BALANCE [MM3/DT]
      @constraint(M,resbalReg0[iMod=1:NMod,k=[1]], res[iMod,k]
                +DT*M3S2MM3*(rel[iMod,k]+spi[iMod,k]) 
                -DT*M3S2MM3*(sum(dis[MyUSMod.USModA[iMod].Dis[iDisUS],k] for iDisUS=1:MyUSMod.USModA[iMod].NDis))
                -DT*M3S2MM3*(sum(byp[MyUSMod.USModA[iMod].Byp[iBypUS],k] for iBypUS=1:MyUSMod.USModA[iMod].NByp))
                -DT*M3S2MM3*(sum(spi[MyUSMod.USModA[iMod].Spi[iSpiUS],k] for iSpiUS=1:MyUSMod.USModA[iMod].NSpi))
                == res0d[iMod] + MyRegFrac*WeekFrac*MyAHData.RegShare[MyHSys.ModNo[iMod]]*(GWH2MAGEFF/MyAHData.EffSea[MyHSys.ModNo[iMod]])*inf)

      #RESERVOIR BALANCE [MM3/DT]
      @constraint(M,resbalReg[iMod=1:NMod,k=2:NK], res[iMod,k]-res[iMod,k-1]
                +DT*M3S2MM3*(rel[iMod,k]+spi[iMod,k]) 
                -DT*M3S2MM3*(sum(dis[MyUSMod.USModA[iMod].Dis[iDisUS],k] for iDisUS=1:MyUSMod.USModA[iMod].NDis))
                -DT*M3S2MM3*(sum(byp[MyUSMod.USModA[iMod].Byp[iBypUS],k] for iBypUS=1:MyUSMod.USModA[iMod].NByp))
                -DT*M3S2MM3*(sum(spi[MyUSMod.USModA[iMod].Spi[iSpiUS],k] for iSpiUS=1:MyUSMod.USModA[iMod].NSpi))
                == MyRegFrac*WeekFrac*MyAHData.RegShare[MyHSys.ModNo[iMod]]*(GWH2MAGEFF/MyAHData.EffSea[MyHSys.ModNo[iMod]])*inf)

      #BALANCE UNREG INFLOW [MM3/DT]
      @constraint(M,balUnReg[iMod=1:NMod,k=1:NK], DT*M3S2MM3*(dis[iMod,k]+byp[iMod,k]-rel[iMod,k]) 
                  == (1.0-MyRegFrac)*WeekFrac*MyAHData.URegShare[MyHSys.ModNo[iMod]]*(GWH2MAGEFF/MyAHData.EffSea[MyHSys.ModNo[iMod]])*inf)

      #HYDROPOWER GENERATION [GWh/Step]
      @constraint(M,prodfunc[iMod=1:NMod,k=1:NK],ghy[iMod,k]
                  -WeekFrac*MW2GWHWEEK*sum(MyAHData.PQData[MyHSys.ModNo[iMod]].Eff[iSeg]*disSeg[iMod,iSeg,k] for iSeg=1:MyAHData.PQData[MyHSys.ModNo[iMod]].NSeg) == 0.0)

      #PROFILE [GWh/Step]
      #   -Ramping according to pre-defined schedule
      if StartRamp > 0
          @constraint(M,updev[k=StartRamp:Int(NK/7):NK],sum(ghy[iMod,k]-ghy[iMod,k-1] for iMod=1:NMod) -ramp >= 0.0)   #Morning up-ramp
          @constraint(M,dndev[k=EndRamp:Int(NK/7):NK],sum(ghy[iMod,k-1]-ghy[iMod,k] for iMod=1:NMod) -ramp >= 0.0)   #Evening down-ramp
      end

      #CAPACITY CONSTRAINTS [GWh/step] - Assume that all stations can contribute
      @constraint(M,upcap[k=1:NK],sum(ghy[iMod,k] for iMod=1:NMod) + resCap <= WeekFrac*MyHSys.MaxProd)
      @constraint(M,dncap[k=1:NK],sum(ghy[iMod,k] for iMod=1:NMod) - resCap >= WeekFrac*MyHSys.MinProd[iWeek])


      #STORAGE AGGREGATION
      @constraint(M,sumres0,sum(MyAHData.EffSea[MyHSys.ModNo[iMod]]*MAGEFF2GWH*res0d[iMod] for iMod=1:NMod) - res0 == 0.0)
      @constraint(M,sumres,sum(MyAHData.EffSea[MyHSys.ModNo[iMod]]*MAGEFF2GWH*res[iMod,NK] for iMod=1:NMod) - resEnd == 0.0)

      #STORAGE DISAGGREGATION
      #   -Scale with regulation degree as a proxy for damping. Avoid hard constraint on end reservoir volumes (!)
      #@constraint(M,disaggr0[iMod=1:NMod],res0d[iMod]-MyAHData.RegDegComp[MyHSys.ModNo[iMod]]*MyAHData.MagShare[MyHSys.ModNo[iMod]]*res0 <= 0.0)
      #@constraint(M,disaggr0[iMod=1:NMod],res0d[iMod]-MyAHData.MagShare[MyHSys.ModNo[iMod]]*res0 <= 0.0)
      #@constraint(M,disaggr[iMod=1:NMod],res[iMod,NK]-MyAHData.RegDegComp[MyHSys.ModNo[iMod]]*MyAHData.MagShare[MyHSys.ModNo[iMod]]*resEnd <= 0.0)

      #RESERVOIR REQUIREMENT [GWh] (GAMMA-V)
      @constraint(M,resreq,resEnd + resPen >= 0.0)

      #ENERGY REQUIREMENT [GWh/Stage] (GAMMA-E)
      @constraint(M,enreq,sum(ghy[iMod,k] for iMod=1:NMod for k=1:NK) + enPen >=  0.0)

      #RAMP REQUIREMENT [GWh/Step/stage] (GAMMA-R)
      @constraint(M,rampreq,ramp + rampPen >=  0.0)

      #RESERVE CAPACITY (GWh/Step/stage) (GAMMA-C)
      @constraint(M,capreq,resCap + capPen >=  0.0)

      #INITIAL RESERVOIR [GWh] (KAPPA-V)
      @constraint(M,resinit,res0 <= 0.0)

      #PHYSICAL INFLOW [GWh] (KAPPA-I)
      @constraint(M,inflow,inf <= 0.0)

      return M
   end
end
