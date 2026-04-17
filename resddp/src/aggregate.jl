using StatsBase
using Printf

#Depth-first upstream and downstream searches along all waterways
function SearchCascade(AHData,USMod,Visited,StartNodes,StartNodesVisited,iArea,iMod)
    if (!(iMod in Visited) && iMod > 0)
        push!(Visited,iMod)
        if (iMod in StartNodes)
            push!(StartNodesVisited,iMod)
        end
        #Upstream
        for iDisUS = 1:USMod[iArea].USModA[iMod].NDis
            SearchCascade(AHData,USMod,Visited,StartNodes,StartNodesVisited,iArea,USMod[iArea].USModA[iMod].Dis[iDisUS])
        end
        for iBypUS = 1:USMod[iArea].USModA[iMod].NByp
            SearchCascade(AHData,USMod,Visited,StartNodes,StartNodesVisited,iArea,USMod[iArea].USModA[iMod].Byp[iBypUS])
        end
        for iSpiUS = 1:USMod[iArea].USModA[iMod].NSpi
            SearchCascade(AHData,USMod,Visited,StartNodes,StartNodesVisited,iArea,USMod[iArea].USModA[iMod].Spi[iSpiUS])
        end
        #Downstream
        mDisTo = 0
        mBypTo = 0
        mSpiTo = 0
        for jMod =  1:AHData[iArea].NMod
            if AHData[iArea].MData[jMod].ModId == AHData[iArea].MData[iMod].DischargeTo
                mDisTo = jMod
            end
            if AHData[iArea].MData[jMod].ModId == AHData[iArea].MData[iMod].BypassTo
                mBypTo = jMod
            end
            if AHData[iArea].MData[jMod].ModId == AHData[iArea].MData[iMod].SpillTo
                mSpiTo = jMod
            end
        end
        SearchCascade(AHData,USMod,Visited,StartNodes,StartNodesVisited,iArea,mDisTo)
        SearchCascade(AHData,USMod,Visited,StartNodes,StartNodesVisited,iArea,mBypTo)
        SearchCascade(AHData,USMod,Visited,StartNodes,StartNodesVisited,iArea,mSpiTo)
    end
end

#Identify cascades in dataset
function IdentifyCascades(AHData,NArea,USMod,ModInfReg,ModInfUReg,CTR,CTI,CNS)
    HydroCascade = Vector{HydroSystemArea}()
    NCascadeTot = 0
    for iArea = 1:NArea
        #Count number of "root nodes", i.e., all waterways to sea. Require inflow.
        NCascade = 0
        StartNodes = Vector{Int}()
        sumInf = 0.0
        for iMod = 1:AHData[iArea].NMod
            MyCnt = AHData[iArea].MData[iMod].ModCnt
            sumInf = sumInf+sum(ModInfReg[MyCnt,1:CTI.NWeek,1:CTI.NInflowYear]+ModInfUReg[MyCnt,1:CTI.NWeek,1:CTI.NInflowYear])
            if sumInf > 0.1 && AHData[iArea].MData[iMod].DischargeTo == 0 && AHData[iArea].MData[iMod].BypassTo == 0 && AHData[iArea].MData[iMod].SpillTo == 0
                NCascade += 1
                push!(StartNodes,iMod)
                sumInf = 0.0
            end
        end
        #Identify modules in each cascade
        Cascades =  Vector{HydroSystemData}()
        StartNodesVisited = Set{Int}()
        for iCascade = 1:NCascade
            iMod = StartNodes[iCascade]
            if !(iMod in StartNodesVisited)
                Visited = Set{Int}()
                myModNos = Vector{Int}()
                myModIds = Vector{Int}()
                myModCnts = Vector{Int}()
                SearchCascade(AHData,USMod,Visited,StartNodes,StartNodesVisited,iArea,iMod)
                cMod = 0
                MinRes = 0.0
                MaxRes = 0.0
                MinProd = zeros(Float64,CTI.NWeek)
                MaxProd = 0.0
                AveInflow = 0.0
                for iMod = 1:AHData[iArea].NMod
                    if (iMod in Visited)
                        MyId = AHData[iArea].MData[iMod].ModId
                        MyCnt = AHData[iArea].MData[iMod].ModCnt
                        MaxRes += CTR.MaxResScale*AHData[iArea].MData[iMod].MaxRes*AHData[iArea].EffSea[iMod]*CNS.MAGEFF2GWH
                        MaxProd += AHData[iArea].MData[iMod].ProdCap*CNS.MW2GWHWEEK #GWh per week
                        for iWeek=1:CTI.NWeek
                            MinProd[iWeek] += minimum(ModInfUReg[MyCnt,iWeek,1:CTI.NInflowYear])*AHData[iArea].EffSea[iMod]*CNS.MAGEFF2GWH 
                            AveInflow += mean(ModInfReg[MyCnt,iWeek,1:CTI.NInflowYear]+ModInfUReg[MyCnt,iWeek,1:CTI.NInflowYear])*AHData[iArea].EffSea[iMod]*CNS.MAGEFF2GWH
                        end
                        push!(myModNos,iMod)
                        push!(myModIds,MyId)
                        push!(myModCnts,MyCnt)
                        cMod += 1
                    end
                end
                MinRes = CTR.ResMinFrac*MaxRes
                push!(Cascades,HydroSystemData(iArea,cMod,myModNos,myModIds,myModCnts,MinRes,MaxRes,MinProd,MaxProd,AveInflow))
            end
        end

        NCascade = length(Cascades)
        push!(HydroCascade,HydroSystemArea(NCascade,Cascades))
        NCascadeTot += NCascade
    end
    return NCascadeTot,HydroCascade
end


function AggrSystems(HCascade,AHData,NArea,ModInfReg,ModInfUReg,CTR,CTI,CNS,CAGR)
    Subsystems = []
    USModSys = []
    NAreaSys = zeros(Int,NArea)
    for iArea = 1:NArea
        NSys = HCascade[iArea].NSys 
        if NSys > 0
            KeepCascade = Vector{Int}()
            for iSys = 1:NSys
                MyCasc = HCascade[iArea].Systems[iSys]
                RegDeg = MyCasc.MaxRes/MyCasc.AveInflow
                Depletion = MyCasc.MaxRes/MyCasc.MaxProd # in weeks
                MaxProd = MyCasc.MaxProd*1000.0/CTI.NHoursWeek # in MW
                
                if (MaxProd > CAGR.ProdCutoff ||  MyCasc.MaxRes > CAGR.ResCutoff) && iArea <= 11 &&  MyCasc.NMod > CAGR.ModCutoff 
                    push!(KeepCascade,iSys)
                end
                #if MyCasc.NMod > CAGR.ModCutoff && MyCasc.MaxRes > CAGR.ResCutoff && MaxProd > CAGR.ProdCutoff && Depletion > CAGR.DeplCutoff
                    #push!(KeepCascade,iSys)
                #end
            end

            #Create an aggregated subsystem comprising cacscades that do not meet the criteria
            myModNos = Vector{Int}()
            myModIds = Vector{Int}()
            myModCnts = Vector{Int}()
            cMod = 0
            MinRes = 0.0
            MaxRes = 0.0
            MinProd = zeros(Float64,CTI.NWeek)
            MaxProd = 0.0
            AveInflow = 0.0
            for iSys = 1:NSys
                if !(iSys in KeepCascade)
                    MyCasc = HCascade[iArea].Systems[iSys]
                    for iMod = 1:MyCasc.NMod
                        MyMod = MyCasc.ModNo[iMod] 
                        MyId = AHData[iArea].MData[MyMod].ModId
                        MyCnt = AHData[iArea].MData[MyMod].ModCnt
                        MaxRes += CTR.MaxResScale*AHData[iArea].MData[MyMod].MaxRes*AHData[iArea].EffSea[MyMod]*CNS.MAGEFF2GWH
                        MaxProd += AHData[iArea].MData[MyMod].ProdCap*CNS.MW2GWHWEEK #GWh per week
                        for iWeek=1:CTI.NWeek
                            MinProd[iWeek] += minimum(ModInfUReg[MyCnt,iWeek,1:CTI.NInflowYear])*AHData[iArea].EffSea[MyMod]*CNS.MAGEFF2GWH 
                            AveInflow += mean(ModInfReg[MyCnt,iWeek,1:CTI.NInflowYear]+ModInfUReg[MyCnt,iWeek,1:CTI.NInflowYear])*AHData[iArea].EffSea[iMod]*CNS.MAGEFF2GWH
                        end
                        push!(myModNos,MyMod)
                        push!(myModIds,MyId)
                        push!(myModCnts,MyCnt)
                        cMod += 1
                    end
                end
            end
            MinRes = CTR.ResMinFrac*MaxRes

            iCnt = 0
            #First copy the unmodified cascades
            for iSys = 1:HCascade[iArea].NSys
                if iSys in KeepCascade
                    iCnt += 1
                    push!(Subsystems,HCascade[iArea].Systems[iSys]) #copy
                end
            end

            #Then add the aggragted cascades
            if cMod > 0
                push!(Subsystems,HydroSystemData(iArea,cMod,myModNos,myModIds,myModCnts,MinRes,MaxRes,MinProd,MaxProd,AveInflow))
                iCnt += 1
            end
            NAreaSys[iArea] = iCnt
        end

        #Adapt datastructure for upstream modules to the aggregated subsystems
        sCnt = length(Subsystems)-NAreaSys[iArea]
        for iSys = 1:NAreaSys[iArea]
        sCnt += 1
        USModC = []
        MySys = Subsystems[sCnt]
        for iMod = 1:MySys.NMod
            UpDis = []; UpByp = []; UpSpi = [];
            MyId = MySys.ModId[iMod]
            for upMod = 1:MySys.NMod
                upModAreaNo = MySys.ModNo[upMod]
                if AHData[iArea].MData[upModAreaNo].DischargeTo == MyId
                    push!(UpDis,upMod)
                end
                if AHData[iArea].MData[upModAreaNo].BypassTo == MyId
                    push!(UpByp,upMod)
                end
                if AHData[iArea].MData[upModAreaNo].SpillTo == MyId
                    push!(UpSpi,upMod)
                end
            end
            LDis = max(0,length(UpDis))
            LByp = max(0,length(UpByp))
            LSpi = max(0,length(UpSpi))
            push!(USModC,UpstreamMod(LDis,LByp,LSpi,UpDis,UpByp,UpSpi))
        end
        push!(USModSys,AreaUpstreamMod(USModC))
        end
    end
    NSysTot = length(Subsystems)
    NMax = maximum(NAreaSys)
    AreaSys = zeros(Int,NArea,NMax)
    iCnt = 1
    for iArea = 1:NArea
        for iSys = 1:NAreaSys[iArea]
        AreaSys[iArea,iSys] = iCnt
        iCnt += 1
        end
    end

    return NSysTot,Subsystems,NAreaSys,AreaSys,USModSys
end


function AggrInflow(NHSys,HSys,ModInfReg,ModInfUReg,AHData,CTR,CTI,CNS)

    #Find energy inflow. Add regulated and unregulated.
    InfReg = zeros(Float64,NHSys,CTI.NWeek,CTI.NInflowYear) 
    InfUReg = zeros(Float64,NHSys,CTI.NWeek,CTI.NInflowYear) 
    RegFrac = zeros(Float64,NHSys,CTI.NWeek) #Average fraction of regulated inflow. 
    for iSys = 1:NHSys
        iArea = HSys[iSys].AreaNo
        for iWeek = 1:CTI.NWeek
            for iYear = 1:CTI.NInflowYear
                for iMod = 1:HSys[iSys].NMod
                    MyCnt = HSys[iSys].ModCnt[iMod]
                    MyMod = HSys[iSys].ModNo[iMod]
                    InfReg[iSys,iWeek,iYear] += (ModInfReg[MyCnt,iWeek,iYear]+ModInfUReg[MyCnt,iWeek,iYear])*AHData[iArea].EffSea[MyMod]*CNS.MAGEFF2GWH 
                    InfUReg[iSys,iWeek,iYear] += ModInfUReg[MyCnt,iWeek,iYear]*AHData[iArea].EffSea[MyMod]*CNS.MAGEFF2GWH 
                end
            end
        RegFrac[iSys,iWeek] = (mean(InfReg[iSys,iWeek,1:CTI.NInflowYear])-mean(InfUReg[iSys,iWeek,1:CTI.NInflowYear]))/mean(InfReg[iSys,iWeek,1:CTI.NInflowYear])
        end
    end
    return InfReg,RegFrac
end


#Print statistics for cascades 
function PrintCascadeStats(dataset,NHSys,HSys,CAGR)
    out = open(string(dataset,"HSysStats.dat"),"w")
    @printf(out,"%s %s %s %s %s %s %s \n","   Area   ", " No   ","NMod   ","MaxRes GWh  ","MaxProd MW     ","Reg. deg.  ","Max. depl. in weeks")
    sumMod = 0
    sumGW = 0.0
    sumTWh = 0.0
    for iSys = 1:NHSys
        MySys = HSys[iSys]
        RegDeg = MySys.MaxRes/MySys.AveInflow
        Depletion = MySys.MaxRes/MySys.MaxProd # in weeks
        MaxProd = MySys.MaxProd*1000/168 # in MW
        sumMod += MySys.NMod
        sumGW += 0.001*MaxProd
        sumTWh += 0.001*MySys.MaxRes
        @printf(out,"%6.0f %6.0f %6.0f %12.2f % 12.2f %12.2f %12.2f \n",MySys.AreaNo,iSys,MySys.NMod,MySys.MaxRes,MaxProd,RegDeg,Depletion)
    end
    @printf(out,"%s \n","------------------------------------------------------------------------")
    @printf(out,"%s %s %s \n","NModTot ","MaxProd GW  ","MaxRes TWh    ")
    @printf(out,"%6.0f %12.2f %12.2f \n",sumMod,sumGW,sumTWh)
    close(out)
end

#Print inflow (for treatment in R-script)
function PrintInflow(dataset,InfReg,NHSys,NWeek,NBranch,NYear)
    open(string(dataset,"InfReg.dat"),"w") do io
        writedlm(io,InfReg)
    end
    out = open(string(dataset,"InfParams.dat"),"w")
    @printf(out,"%8.0f %8.0f %8.0f %8.0f \n",NHSys,NWeek,NBranch,NYear)
    close(out)
end
