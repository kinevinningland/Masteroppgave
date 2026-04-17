#Aggregated stage problem 
module StageProb
   using JuMP 
   using MathOptInterface
   const MOI = MathOptInterface

   function Build(t,iWeek,IM,NHSys,HSys,AreaSys,AMData,EV,CCR,CCH,CCI,CFP,CFR,CFC,CFHD,CFHS,CFRHS,CNS,NCut,NCostCut,NArea,CTI,LEndVal,
           LFeasSpace,NFeasCut,FCC,H2Data,optimizer)

      M = Model(optimizer)

      WeekFrac = CTI.WeekFrac
      NH2Area = H2Data.NArea

      @variable(M,0.0 <= ner[iSys=1:NHSys] <= CNS.Big,base_name="ner")                                                     # GWh
      @variable(M,0.0 <= mir[iSys=1:NHSys] <= HSys[iSys].MinRes, base_name="mir")                                          # GWh
      @variable(M,0.0 <= res[iSys=1:NHSys] <= HSys[iSys].MaxRes, base_name="res")                                          # GWh
      @variable(M,0.0 <= spi[iSys=1:NHSys] <= CNS.Big, base_name="spi")                                                    # GWh/stage
      @variable(M, HSys[iSys].MinProd[iWeek] <= prod[iSys=1:NHSys] <=  HSys[iSys].MaxProd,base_name="prod")                # GWh/stage
      @variable(M, 0.0 <= cap[iSys=1:NHSys] <=  0.50*WeekFrac*(HSys[iSys].MaxProd-HSys[iSys].MinProd[iWeek]),base_name="cap")   # GWh/step/stage
      @variable(M, 0.0 <= ramp[iSys=1:NHSys] <=  WeekFrac*(HSys[iSys].MaxProd-HSys[iSys].MinProd[iWeek]),base_name="ramp") # GWh/step
      @variable(M,0.0 <= fsl[iSys=1:NHSys,c=1:NFeasCut[iSys]] <= CNS.Big, base_name="fsl")                                 # GWh
      @variable(M,zinit[iSys=1:NHSys],base_name="zinit")                                                                   # []   
      @variable(M,inf[iSys=1:NHSys],base_name="inf")                                                                       # GWh
      @variable(M,rinit[iSys=1:NHSys],base_name="rinit")                                                                   # GWh
      @variable(M,-CNS.AlphaMax <= alpha <= CNS.AlphaMax,base_name="alp")                                                  # 10E3 EUR
      @variable(M,-CNS.AlphaMax <= beta <= CNS.AlphaMax,base_name="bet")                                                   # 10E3 EUR
      @variable(M,0.0 <= h2res[iArea=1:NH2Area] <= H2Data.Areas[iArea].MaxRes, base_name="h2res")                          # GWh
      @variable(M,h2init[iArea=1:NH2Area],base_name="h2init")                                                              # GWh
      @variable(M, -CTI.NK*H2Data.Areas[iArea].MaxDis <= h2dis[iArea=1:NH2Area] <= CTI.NK*H2Data.Areas[iArea].MaxDis, base_name="h2dis") # GWh/stage

      #Min Cost [10E3 EUR]
      @objective(M,MathOptInterface.MIN_SENSE,
                 +sum(CNS.CNegRes*ner[iSys] for iSys=1:NHSys)
                 +sum(CNS.CMinRes*mir[iSys] for iSys=1:NHSys)
                 +sum(CNS.CSpi*spi[iSys] for iSys=1:NHSys)
                 +sum(CNS.CFeas*fsl[iSys,c] for iSys=1:NHSys for c=1:NFeasCut[iSys])
                 +alpha+beta)
     
      #Inflow = Phi*Sigma*zinit+Sigma*Eps(t) + MU. First term on LHS. DevCorrMat is the coefficient Phi(NHSys*NHSys).*Sigma(NHSys*1).
      DevCorrMat = IM.InflowSDev[1:NHSys,iWeek].*IM.CorrMat #Dim(NSer*NSer)
      MuCorrVec = zeros(Float64,NHSys,NCut)
      for c = 1:NCut
         MuCorrVec[1:NHSys,c] = IM.CorrMat*CCI[1:NHSys,t,c]
      end

      #Reservoir balance [GWh] 
      @constraint(M,resbal0[iSys=1:NHSys],res[iSys]+prod[iSys]+spi[iSys] -ner[iSys] -rinit[iSys] -inf[iSys] == 0.0) 

      #Minimum reservoir limit soft [GWh] 
      @constraint(M,resminSoft[iSys=1:NHSys],res[iSys]+mir[iSys] >= HSys[iSys].MinRes) 

      #Reservoir state: V_{t-1}
      @constraint(M,rstate[iSys=1:NHSys],rinit[iSys] == 0.0)

      #Normalized inflow: Z_{t-1} 
      @constraint(M,zstate[iSys=1:NHSys],zinit[iSys] == 0.0)

      #Physical inflow [GWh]
      @constraint(M,inflow[iSys=1:NHSys],inf[iSys] - sum(DevCorrMat[iSys,jSys]*zinit[jSys] for jSys=1:NHSys) <= 0.0)

      #Capacity constraints  [GWh/stage] 
      @constraint(M,upcap[iSys=1:NHSys],prod[iSys] + (1.0/WeekFrac)*cap[iSys] <= HSys[iSys].MaxProd)
      @constraint(M,dncap[iSys=1:NHSys],prod[iSys] - (1.0/WeekFrac)*cap[iSys] >= HSys[iSys].MinProd[iWeek])

      #H2 storage [GWh] 
      @constraint(M,h2storage[iArea=1:NH2Area],h2res[iArea]+h2dis[iArea]-h2init[iArea] == 0.0)

      #H2 storage state: s_{t-1}
      @constraint(M,h2state[iArea=1:NH2Area],h2init[iArea] == 0.0)

      #Cost cuts [10E3 EUR]
      @constraint(M,costcut[c=1:NCostCut],beta+sum(CFP[iSys,t,c]*prod[iSys] for iSys=1:NHSys)+sum(CFR[iSys,t,c]*ramp[iSys] for iSys=1:NHSys)+sum(CFC[iSys,t,c]*cap[iSys] for iSys=1:NHSys) 
                  +sum(CFHD[iArea,t,c]*h2dis[iArea] for iArea=1:NH2Area)+sum(CFHS[iArea,t,c]*h2init[iArea] for iArea=1:NH2Area) >= CFRHS[t,c])

      #Feasibility constraints []
      if LFeasSpace
         @constraint(M,feas[iSys=1:NHSys,c=1:NFeasCut[iSys]], FCC[iSys,c,1]*res[iSys] + FCC[iSys,c,2]*prod[iSys] 
                     +FCC[iSys,c,3]*ramp[iSys] +FCC[iSys,c,4]*cap[iSys] + FCC[iSys,c,5]*rinit[iSys] + FCC[iSys,c,6]*inf[iSys] -fsl[iSys,c] <= FCC[iSys,c,7])
      end

      if LEndVal
         if EV.NEndCut > 0
             @constraint(M,endc[c=1:EV.NEndCut],alpha+sum(EV.EndCutCoef[c,iSys]*res[iSys] for iSys=1:NHSys) >= EV.EndCutRHS[c])
         elseif EV.NEndSeg > 0
             SegFrac = 1.0/Float64(EV.NEndSeg)
             @variable(M,0.0 <= res_seg[iSys=1:NHSys,iSeg=1:EV.NEndSeg] <= SegFrac*HSys[iSys].MaxRes, base_name="res_seg") 
             @constraint(M,endsegval,alpha+sum(EV.EndSegCoef[iSeg]*res_seg[iSys,iSeg] for iSys=1:NHSys for iSeg=1:EV.NEndSeg) >= 0.0)
             @constraint(M,endsegbal[iSys=1:NHSys],sum(res_seg[iSys,iSeg] for iSeg=1:EV.NEndSeg) - res[iSys] == 0.0)
         else
            @constraint(M,endset,alpha == 0)
         end
      else
          @constraint(M,cut[c=1:NCut],alpha-sum(CCR[iSys,t,c]*res[iSys] for iSys=1:NHSys)-sum(CCH[iArea,t,c]*h2res[iArea] for iArea=1:NH2Area)-sum(MuCorrVec[iSys,c]*zinit[iSys] for iSys=1:NHSys) >= 0.0)
      end


      return M
   end
end

