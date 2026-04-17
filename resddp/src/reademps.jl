#Preparing data for samplan

using HDF5
using StatsBase
using Printf

function ReadHDF5(dataset,CNS,CTI,CTR)
    t_setup = time_ns()
    DT = CTI.DT
    DT_int = Int(DT)
    NK = CTI.NK
    NWeek = CTI.NWeek
    NInflowYear = CTI.NInflowYear
    NHoursWeek = CTI.NHoursWeek
    NHoursWeek_int = Int(NHoursWeek)

    #AREA DATA#
    println("Reading from dataset: ",dataset)
    filename = joinpath(dataset,"model.h5")
    fid = h5open(filename,"r")

    dset = fid["area_names"]
    Names = read(dset)
    NArea = length(Names)
    AreaName = String[]
    AreaNumber = zeros(Int,NArea)
    for iArea =1:NArea
        push!(AreaName,strip(Names[iArea].name))
        AreaNumber[iArea] = Names[iArea].number
    end

    dset = fid["area_connections/connections"]
    AreaConnections = read(dset)
    NCon = length(AreaConnections)
    ConnectionName = String[]
    ConnectionTopo = zeros(Int,NCon,2)
    for iCon =1:NCon
        push!(ConnectionName,AreaConnections[iCon].line_id)
        ConnectionTopo[iCon,1] = AreaConnections[iCon].from_area
        ConnectionTopo[iCon,2] = AreaConnections[iCon].to_area
    end

    #WIND POWER PARKS
    dset = fid["market_data/wind_power"]
    WPPark = read(dset)
    NWPPark = length(WPPark)
    WPParkNames = String[]
    WPParkArea = zeros(Int,NWPPark)
    WPParkPeriod = zeros(Int,NWPPark)
    for iPark = 1:NWPPark
        push!(WPParkNames,strip(WPPark[iPark].name))
        WPParkArea[iPark] = WPPark[iPark].area_id
        WPParkPeriod[iPark] = WPPark[iPark].periode
    end

    dset = fid["market_data/list_of_periodes/numb_in_list"]
    NWindPer = read(dset)[1]
    WPParkWeight = zeros(Float64,NWindPer)
    for iPer = 1:NWindPer
        tag = string(string("market_data/list_of_periodes/Periode_",iPer),"/weight")
        locDset = fid[tag]
        WPParkWeight[iPer] = Float64(read(locDset)[1])
    end


    #HYDRO MODULES
    dset = fid["hydro_data"]
    HydroAreas = read(dset)
    NAreaHyd = length(HydroAreas)-1
    NAreaTerm = NArea-NAreaHyd
    MyKeys = HydroAreas.keys

    for iKey = 1:length(MyKeys)
        if ! isassigned(MyKeys,iKey)
            MyKeys[iKey] = " "
        end
    end
        
    modCnt = 0
    MaxModArea = 0
    ModNr = []
    AHData = []
    WMData = []
    for iArea = 1:NArea
        MData = []
        if AreaName[iArea] in MyKeys
            helptag = string("hydro_data/",AreaName[iArea])

            #INFLOW WATERMARKS#
            tag = string(helptag,"/Watermark_data")
            locDset = fid[tag]
            WaterMarks = read(locDset)
            NWaterMark = length(WaterMarks)
            InfNames = String[]
            Scaling = zeros(Float64,NWaterMark)
            AveInflow = zeros(Float64,NWaterMark)
            for iMark = 1:NWaterMark
                push!(InfNames,strip(WaterMarks[iMark].infl_name))
                Scaling[iMark] = WaterMarks[iMark].scaling_factor
                AveInflow[iMark] = WaterMarks[iMark].average_inflow
            end
            push!(WMData,WaterMarkData(AreaName[iArea],NWaterMark,Scaling,AveInflow,InfNames))

            #HYDROPOWER DATA#
            tag = string(helptag,"/Module_data")
            locDset = fid[tag]
            ModData = read(locDset)
            NMod = length(ModData)
            if NMod > MaxModArea 
                MaxModArea=NMod 
            end

            for iMod = 1:NMod
                UInfName = strip(ModData[iMod].u_infl_name)
                RInfName = strip(ModData[iMod].r_infl_name)
                PlantName = strip(ModData[iMod].plant_name)
                ResName = strip(ModData[iMod].res_name)
                RInfVol = ModData[iMod].r_infl_rvol
                UInfVol = ModData[iMod].u_infl_rvol
                RegDeg = ModData[iMod].reg_level
                ResBottom = ModData[iMod].res_bottom
                MaxBypass = ModData[iMod].max_bypass
                MaxFlow = ModData[iMod].max_flow
                NomElevation = ModData[iMod].nom_elevation
                SubmersionLevel = ModData[iMod].Undervannstand
                MaxRes = ModData[iMod].max_res
                ProdCap = ModData[iMod].prod_cap
                ConvFac = ModData[iMod].conv_factor
                ModId = ModData[iMod].res_id
                InfLastYr = ModData[iMod].infl_l_year
                InfFirstYr = ModData[iMod].infl_f_year
                CoplNumber = ModData[iMod].copl_number
                SpillTo = ModData[iMod].spill_to
                BypassTo = ModData[iMod].bypass_to
                DischargeTo = ModData[iMod].flow_to

                if CoplNumber > 0
                    #println("Hydraulic coupling in Area ",iArea," and module ",iMod," is not treated!")
                end

                #Create a module object
                modCnt += 1
                push!(ModNr,modCnt)
                push!(MData,ModuleData(AreaName[iArea],UInfName,RInfName,PlantName,ResName,UInfVol,RInfVol,RegDeg,ResBottom,
                                        MaxBypass,MaxFlow,NomElevation,SubmersionLevel,MaxRes,ProdCap,
                                        ConvFac,ModId,modCnt,InfFirstYr,InfLastYr,CoplNumber,SpillTo,BypassTo,DischargeTo))
            end

            #EFFICIENCY CURVES#
            PQTol = 1.0E-5
            PQPointTol = 1.0E-4
            PQCData = Array{PQCurveData}(undef,NMod)
            tag = string(helptag,"/PQ_curve/")
            mytag = string(tag,"npoint")
            locDset = fid[mytag]
            NCurve = read(locDset)[1]+1

            for iCurve = 1:NCurve
                curveTag = string(string(tag,string("curve_",iCurve-1)))
                mytag = string(curveTag,"/nPkt")
                NSeg = read(fid[mytag])[1]
                mytag = string(curveTag,"/M3s")
                M3S = read(fid[mytag])
                mytag = string(curveTag,"/MW")
                MW = read(fid[mytag])
                mytag = string(curveTag,"/PMax")
                PMax = read(fid[mytag])[1]
                mytag = string(curveTag,"/QMax")
                QMax = read(fid[mytag])[1]
                mytag = string(curveTag,"/first_week")
                FirstWeek = read(fid[mytag])[1]
                mytag = string(curveTag,"/last_week")
                LastWeek = read(fid[mytag])[1]
                mytag = string(curveTag,"/id")
                MyId = read(fid[mytag])[1]

                #Find efficiency
                BestEff = 0.0
                Eff = zeros(Float64,NSeg)
                DMax = zeros(Float64,NSeg)
                for iSeg = 1:NSeg
                    #Ensure strictly increasing PQ-curve
                    if iSeg > 1
                        if MW[iSeg] <= MW[iSeg-1] 
                            MW[iSeg] = MW[iSeg-1]+PQTol
                        end
                        if M3S[iSeg] <= M3S[iSeg-1]
                            M3S[iSeg] = M3S[iSeg-1]+PQTol
                        end
                    end
                    if iSeg == 1
                        if M3S[iSeg] > PQTol
                            Eff[iSeg] = MW[iSeg]/M3S[iSeg]
                        end
                        DMax[iSeg] = max(M3S[iSeg],0.0)
                    else
                        Eff[iSeg] = (MW[iSeg]-MW[iSeg-1])/(M3S[iSeg]-M3S[iSeg-1])
                        DMax[iSeg] = M3S[iSeg]-M3S[iSeg-1]
                    end
                    EffZero = 0.0
                    if M3S[iSeg] > PQTol
                        EffZero = MW[iSeg]/M3S[iSeg]
                        if EffZero > BestEff
                            BestEff = EffZero
                        end
                    end
                end


                MaxDis = sum(DMax)
                if NSeg > 0
                    PMax = MW[NSeg]
                    PMin = MW[1]
                else
                    PMax = 0.0
                    PMin = 0.0
                end

                if PMax <= PMin
                    PMin = 0.0
                end
                
                for iMod = 1:NMod
                    #Sort PQ-curves according to module number
                    if MData[iMod].ModId == MyId
                        PQCData[iMod] = PQCurveData(NSeg,MyId,M3S,DMax,MW,Eff,BestEff,PMax,PMin,QMax,MaxDis,FirstWeek,LastWeek,iCurve-1)
                        continue
                    end
                end
            end
            #Insert dummy-stations for modules missing PQ-curve
            dNSeg=0;dMyId=0;dM3S=Array{Float64}(undef,1);dDMax=Array{Float64}(undef,1);dMW=Array{Float64}(undef,1);dEff=Array{Float64}(undef,1);
            dBestEff=0.0;dPMax=0.0;dPMin=0.0;dQMax=0.0;dMaxDis=0.0;dFirstWeek=0;dLastWeek=0;diCurve=0;
            for iMod = 1:NMod
                if !isassigned(PQCData,iMod)
                    PQCData[iMod] = PQCurveData(dNSeg,dMyId,dM3S,dDMax,dMW,dEff,dBestEff,dPMax,dPMin,dQMax,dMaxDis,dFirstWeek,dLastWeek,diCurve)
                end
            end

            #RESERVOIR CURVES#
            RCData = Array{ResCurveData}(undef,NMod)
            rctag = string(helptag,"/res_curves/")
            for iMod = 1:NMod
                curvetag = string(rctag,MData[iMod].ModId)
                mytag = string(curvetag,"/nPkt")
                NPkt = read(fid[mytag])[1]
                if NPkt > 0
                    tag = string(curvetag,"/Kote")
                    Elevation = read(fid[tag])
                    tag = string(curvetag,"/Vol")
                    Volume = read(fid[tag])
                    RCData[iMod] = ResCurveData(NPkt,Elevation,Volume)
                end
            end

            #Find energy efficiency to sea
            EffSea = 1.0E-3*ones(Float64,NMod) #A lower level introduced to avoid NaN
            for iMod = 1:NMod
                iTo = iMod
                while (iTo > 0)
                    EffSea[iMod] += PQCData[iTo].BestEff
                    NextDs = 0
                    for jMod =  1:NMod
                        if MData[jMod].ModId == MData[iTo].DischargeTo
                            NextDs = jMod
                            break
                        end
                    end
                    iTo = NextDs
                end
            end


            #MODULE SHARE (of energy in area)
            MagShare = zeros(Float64,NMod)
            RegShare = zeros(Float64,NMod)
            URegShare = zeros(Float64,NMod)
            SumMaxGWH = 0.0
            SumRegGWH = 0.0
            SumURegGWH = 0.0
            for iMod = 1:NMod
                SumMaxGWH += ModData[iMod].max_res*(EffSea[iMod]*CNS.MAGEFF2GWH)
                SumRegGWH += ModData[iMod].r_infl_rvol*(EffSea[iMod]*CNS.MAGEFF2GWH) 
                SumURegGWH += ModData[iMod].u_infl_rvol*(EffSea[iMod]*CNS.MAGEFF2GWH) 
            end
            for iMod = 1:NMod
                if EffSea[iMod] > 1.0E-3 && ModData[iMod].max_res > 1.0E-1 
                   MagShare[iMod] = ModData[iMod].max_res*(EffSea[iMod]*CNS.MAGEFF2GWH)/SumMaxGWH
                end
                if MagShare[iMod] < 1.0E-3
                   MagShare[iMod] = 0.0
                end
                if SumRegGWH > 1.0E-3
                   RegShare[iMod] =  ModData[iMod].r_infl_rvol*(EffSea[iMod]*CNS.MAGEFF2GWH)/SumRegGWH
                end
                if RegShare[iMod] < 1.0E-3
                   RegShare[iMod] = 0.0
                end
                if SumURegGWH > 1.0E-3
                   URegShare[iMod] =  ModData[iMod].u_infl_rvol*(EffSea[iMod]*CNS.MAGEFF2GWH)/SumURegGWH
                end
                if URegShare[iMod] < 1.0E-3
                   URegShare[iMod] = 0.0
                end
            end

            #Regulation degree 
            RegDegComp = zeros(Float64,NMod)
            for iMod = 1:NMod
                aInflow = MData[iMod].RInfVol
                mRes = MData[iMod].MaxRes
                
                #Consider immediate upstream modules only
                for uMod = 1:NMod
                    if MData[uMod].DischargeTo == MData[iMod].ModId
                        aInflow += MData[uMod].UInfVol
                        if MData[uMod].MaxRes < 1.0
                            aInflow += MData[uMod].RInfVol
                        end
                    end
                end

                if mRes > 10.0 && aInflow > 5.0
                    RegDegComp[iMod] = max(mRes/aInflow,0.10)
                end
            end

            RampFrac = ones(Float64,NMod)
            rampfile = joinpath(dataset,"RampFrac.dat")
            rampfid = open(rampfile,"r")
            NRamp = Int(ceil(length(readlines(rampfid))))
            close(rampfid)
            rampfid = open(rampfile,"r")
            for iRamp = 1:NRamp
                line = readline(rampfid); items = split(line,","); 
                if parse(Int,items[1]) == iArea
                    iMod = parse(Int, items[2])
                    tmp = parse(Float64, items[3])
                    RampFrac[iMod] = tmp
                end
            end
            close(rampfid)

            push!(AHData,AreaDataHydro(NMod,EffSea,RegDegComp,MagShare,RegShare,URegShare,RampFrac,MData,PQCData,RCData))
        
        else #No hydro, make dummy structures
            UInfName = " "
            RInfName = " "
            PlantName = " "
            ResName = " "
            UInfVol = 0.0
            RInfVol = 0.0
            RegDeg = 0.0
            ResBottom = 0.0
            MaxBypass = 0.0
            MaxFlow = 0.0
            NomElevation = 0.0
            SubmersionLevel = 0.0
            MaxRes = 0.0
            ProdCap = 0.0
            ConvFac = 0.0
            ModId = 0
            InfFirstYr = 0
            InfLastYr = 0
            CoplNumber = 0
            SpillTo = 0
            BypassTo = 0
            DischargeTo = 0
            push!(MData,ModuleData(AreaName[iArea],UInfName,RInfName,PlantName,ResName,UInfVol,RInfVol,RegDeg,ResBottom,
                MaxBypass,MaxFlow,NomElevation,SubmersionLevel,MaxRes,ProdCap,ConvFac,ModId,ModId,
                InfFirstYr,InfLastYr,CoplNumber,SpillTo,BypassTo,DischargeTo))

            NWaterMark = 0
            Scaling = Array{Float64}(undef,1)
            AveInflow = Array{Float64}(undef,1)
            InfNames = Array{String}(undef,1)
            push!(WMData,WaterMarkData(AreaName[iArea],NWaterMark,Scaling,AveInflow,InfNames))

            NMod = 0
            PQCData = Array{PQCurveData}(undef,1)
            dNSeg=0;dMyId=0;dM3S=Array{Float64}(undef,1);dDMax=Array{Float64}(undef,1);dMW=Array{Float64}(undef,1);dEff=Array{Float64}(undef,1);
            dBestEff=0.0;dPMax=0.0;dPMin=0.0;dQMax=0.0;dMaxDis=0.0;dFirstWeek=0;dLastWeek=0;diCurve=0;
            PQCData[1] = PQCurveData(dNSeg,dMyId,dM3S,dDMax,dMW,dEff,dBestEff,dPMax,dPMin,dQMax,dMaxDis,dFirstWeek,dLastWeek,diCurve)
            RCData = Array{ResCurveData}(undef,1)
            dNPkt=0;dEl=Array{Float64}(undef,1);dVol=Array{Float64}(undef,1);
            RCData[1] = ResCurveData(dNPkt,dEl,dVol)

            dA = zeros(Float64,NMod)
            push!(AHData,AreaDataHydro(NMod,dA,dA,dA,dA,dA,dA,MData,PQCData,RCData))
        end
    end

    for iArea = 1:NArea
        if AreaName[iArea] in MyKeys
            helptag = string("hydro_data/",AreaName[iArea])
            tag = string(helptag,"/Pumpe_data")
            if false && exists(fid,tag)
                locDset = fid[tag]
                PUData = read(locDset)
                NPumpArea = length(PUData)
                if (NPumpArea > 0) 
                    println(NPumpArea, " in area ", iArea, " are not considered")
                end
            end 
        end
    end


    NModTot = modCnt
    USMod = []
    for iArea = 1:NArea
        USModArea = []
        if AreaName[iArea] in MyKeys
            for iMod = 1:AHData[iArea].NMod
                UpDis = []; UpByp = []; UpSpi = [];
                MyId = AHData[iArea].MData[iMod].ModId
                for upMod = 1:AHData[iArea].NMod
                    if AHData[iArea].MData[upMod].DischargeTo == MyId
                        push!(UpDis,upMod)
                    end
                    if AHData[iArea].MData[upMod].BypassTo == MyId
                        push!(UpByp,upMod)
                    end
                    if AHData[iArea].MData[upMod].SpillTo == MyId
                        push!(UpSpi,upMod)
                    end
                end
                LDis = max(0,length(UpDis))
                LByp = max(0,length(UpByp))
                LSpi = max(0,length(UpSpi))
                push!(USModArea,UpstreamMod(LDis,LByp,LSpi,UpDis,UpByp,UpSpi))
            end
        end
        push!(USMod,AreaUpstreamMod(USModArea))
    end

    close(fid) #close model.h5
    println("Read model.h5")


    #INFLOW DATA#
    #Store module inflow in [Mm3/week]
    #Formula: scale*(ave_mod/ave_wm)*inf_wm
    filename = joinpath(dataset,"historical.h5")
    fid = h5open(filename,"r")
    ModInfReg = zeros(Float64,NModTot,NWeek,NInflowYear) 
    ModInfUReg = zeros(Float64,NModTot,NWeek,NInflowYear) 
    for iArea = 1:NArea
        if AreaName[iArea] in MyKeys
            for iMod = 1:AHData[iArea].NMod
                MyNr = AHData[iArea].MData[iMod].ModCnt
                #Regulated inflow
                MyName = AHData[iArea].MData[iMod].RInfName
                locDset = fid[string("/historical_series/",MyName)]
                Indx = findfirst(x->x==MyName,WMData[iArea].SeriesNames)
                WMScale = WMData[iArea].Scaling[Indx]
                WMAverage = WMData[iArea].AverageInflow[Indx]
                startTime = attrs(locDset)["start_time"]
                timeResol = parse(Int,attrs(locDset)["TS_object_type"][1])
                MyScale = WMScale*AHData[iArea].MData[iMod].RInfVol/WMAverage
                WMInflow = read(locDset)
                #println(iArea," ",iMod," ",MyNr," ",MyName," ",AHData[iArea].MData[iMod].RInfVol," ",MyScale," ",WMAverage," ")
                if timeResol == 124     #Daily values in m3/s, scaling_factor from model.h5 turns this into Mm3/day
                    iCnt = 1
                    for iYear = 1:NInflowYear
                        for iWeek = 1:NWeek
                            ModInfReg[MyNr,iWeek,iYear] = sum(WMInflow[iCnt:(iCnt+6)])*MyScale
                            iCnt += 7
                        end
                    end
                elseif timeResol == 123 #Week values in Mm3/week, scaling_factor from model.h5 turns this into Mm3/day
                    iCnt = 1
                    for iYear = 1:NInflowYear
                        for iWeek = 1:NWeek
                            ModInfReg[MyNr,iWeek,iYear] = WMInflow[iCnt]*MyScale
                            iCnt += 1
                        end
                    end
                else
                    error(println("Invalid TS_object_type ",timeResol))
                end

                #Unregulated inflow
                MyName = AHData[iArea].MData[iMod].UInfName
                locDset = fid[string("/historical_series/",MyName)]
                Indx = findfirst(x->x==MyName,WMData[iArea].SeriesNames)
                WMScale = WMData[iArea].Scaling[Indx]
                WMAverage = WMData[iArea].AverageInflow[Indx]
                startTime = attrs(locDset)["start_time"]
                timeResol = parse(Int,attrs(locDset)["TS_object_type"][1])
                MyScale = WMScale*AHData[iArea].MData[iMod].UInfVol/WMAverage
                WMInflow = read(locDset)
                if isnan(MyScale); MyScale = 0.0; end
                if timeResol == 124     #Daily values in m3/s, scaling_factor from model.h5 turns this into Mm3/day
                    iCnt = 1
                    for iYear = 1:NInflowYear
                        for iWeek = 1:NWeek
                            ModInfUReg[MyNr,iWeek,iYear] = sum(WMInflow[iCnt:(iCnt+6)])*MyScale
                            iCnt += 7
                        end
                    end
                elseif timeResol == 123 #Week values in Mm3/week, scaling_factor from model.h5 turns this into Mm3/day
                    iCnt = 1
                    for iYear = 1:NInflowYear
                        for iWeek = 1:NWeek
                            ModInfUReg[MyNr,iWeek,iYear] = WMInflow[iCnt]*MyScale
                            iCnt += 1
                        end
                    end
                else
                    error(println("Invalid TS_object_type ",timeResol))
                end
                WMInflow = read(locDset)

                if isnan(ModInfReg[MyNr]); ModInfReg[MyNr] = 0.0; end
                if isnan(ModInfUReg[MyNr]); ModInfUReg[MyNr] = 0.0; end
            end
        end
    end

    close(fid) #close historical.h5
    println("Read historical.h5")


    # SCENARIO DATA
    # Could/should use this file to read all inflow data. Currently used for reading StartYear.
    FirstYear = 0
    filename = joinpath(dataset,"ScenarioData.h5")
    fid = h5open(filename,"r")
    for iArea = 1:NArea
        if AreaName[iArea] in MyKeys
            for iMod = 1:AHData[iArea].NMod
                MyNr = AHData[iArea].MData[iMod].ModCnt
                MyName = AHData[iArea].MData[iMod].RInfName
                locDset = fid[string(string(AreaName[iArea],"/"),MyName)]
                tmp = parse(Int,attrs(locDset)["STAAR"][1])
                if tmp != FirstYear
                    println("FirstYear = ",tmp)
                    FirstYear = tmp
                end
            end
        end
    end

    close(fid) #close ScenarioData.h5
    println("Read ScenarioData.h5")

    #MARKET DATA#
    #frequency = {123:'W',124:'D',125:'H'}
    filename = joinpath(dataset,"TidsserieData.h5")
    fid = h5open(filename,"r")
    AMData = []

    for iArea = 1:NArea
        #Market steps
        trinnTag = string(AreaName[iArea],"/TRINN")
        AllKeys = read(fid[trinnTag]).keys
        for iKey = 1:length(AllKeys)
            if ! isassigned(AllKeys,iKey)
                AllKeys[iKey] = " "
            end
        end

        NMStep = read(fid[string(trinnTag,"/nTrinn")])[1]
        MSData = Array{MarketStepData}(undef,NMStep)
        stepCnt = 0
        for iStep = 1:NMStep
            testTag = string("Trinn_",iStep)
            if testTag in AllKeys
                iTag = string(trinnTag,string("/Trinn_",iStep))
                timeResol = parse(Int,attrs(fid[string(iTag,"/Mengde")])["TS_object_type"][1])
                typeNr = parse(Int,attrs(fid[string(iTag)])["Id"][1])
                marketStepName = attrs(fid[string(iTag)])["Navn"][1]
                if timeResol == 125 #hourly
                    Capacity = read(fid[string(iTag,"/Mengde")]) #GWh
                    Price = read(fid[string(iTag,"/Pris")]) # [10E3 MU/GWh] 
                    stepCnt += 1
                    WeekCapacity = zeros(Float64,NWeek)
                    WeekPrice = zeros(Float64,NWeek)
                    hourFrom = 1
                    hourTo = NHoursWeek_int
                    for iWeek = 1:NWeek
                        WeekCapacity[iWeek] = sum(Capacity[hourFrom:hourTo]) 
                        WeekPrice[iWeek] = mean(Price[hourFrom:hourTo]) 
                        hourFrom += NHoursWeek_int
                        hourTo += NHoursWeek_int
                    end
                    MSData[stepCnt] = MarketStepData(marketStepName, WeekCapacity[1:NWeek],WeekPrice[1:NWeek]) 
                else
                    error(println("This time resolution is not treated!"))
                end
            end
        end
        NMStep = stepCnt
        
        #LOAD
        #Data is originally in GWh/step
        loadTag = string(AreaName[iArea],"/DELLASTER")
        NLoad = read(fid[string(loadTag,"/nDellaster")])[1]
        MLData = Array{MarketLoadData}(undef,NLoad)

        for iLoad = 1:NLoad
            iTag = string(loadTag,string("/DELLAST_",iLoad))
            timeResol = parse(Int,attrs(fid[iTag])["TS_object_type"][1])
            marketLoadName = attrs(fid[iTag])["TS_name"][1]
            if timeResol == 125
                Load = read(fid[iTag])
                WeekLoad = zeros(Float64,NWeek,NK)
                hourFrom = 1
                hourTo = DT_int
                for iWeek = 1:NWeek
                    for k = 1:NK
                        WeekLoad[iWeek,k] = sum(Load[hourFrom:hourTo]) 
                        hourFrom += DT_int
                        hourTo += DT_int
                    end
                end
                MLData[iLoad] = MarketLoadData(marketLoadName, CTR.LoadScale*WeekLoad[1:NWeek,1:NK]) 
            else
                error(println("This time resolution is not treated!"))
            end
        end
        push!(AMData,AreaDataMarket(NMStep,NLoad,MSData,MLData))
    end
    close(fid) #close TidsserieData.h5
    println("Read TidsserieData.h5")


    #WIND POWER DATA
    #data are provided in GWh/hour. Sum to GWh/time step. 
    filename = joinpath(dataset,"wind_data.h5")
    fid = h5open(filename,"r")
    WPData = zeros(Float64,NArea,NInflowYear,NWeek,NK)

    for iPark = 1:NWPPark
        iArea = WPParkArea[iPark]
        iPeriod = WPParkPeriod[iPark]
        myWeight = WPParkWeight[iPeriod]
        iYear = 1
        for yearCnt = FirstYear:(FirstYear+NInflowYear-1)
            locDset = fid[string(uppercase(WPParkNames[iPark]),string("/",yearCnt))]
            WPD = read(locDset)
            timeResol = parse(Int,attrs(locDset)["TS_object_type"][1])
            if timeResol == 125 #hourly
                iHour = 1
                for iWeek = 1:NWeek
                    for k = 1:NK
                        WPData[iArea,iYear,iWeek,k] += sum(WPD[iHour:(iHour+DT_int-1)])*myWeight 
                        iHour += DT_int
                    end
                end
            elseif timeResol == 124 #daily
                error(println("Daily resolution in wind not implemented"))
            elseif timeResol == 125 #weekly
                error(println("Weekly resolution in wind not implemented"))
            end
            iYear += 1
        end
    end

    close(fid)
    println("Read wind_data.h5")

    #println(" ")
    #println("Reading time: ",(time_ns()-t_setup)*1.0E-9)
    #println(" ")

    return AHData,AMData,USMod,WPData,NArea,ModInfReg,ModInfUReg,AreaName,MyKeys,MaxModArea
end

#HYDROGEN MODEL
function ReadH2(dataset,AHData,NArea,AreaName,CNS,CTI,CTR)
    H2Areas = []
    NH2Area = 0
    H2Ind = zeros(Int,NArea)
    for iArea = 1:NArea
        if contains(AreaName[iArea],"H2")
            H2MaxDis =  AHData[iArea].MData[1].MaxFlow*CNS.M3S2MM3*CTI.DT  #GWh/step
            H2MaxRes = AHData[iArea].MData[1].MaxRes                       #GWh
            LStrategic = false
            if H2MaxRes > 2.0
                LStrategic = true
            end
            NH2Area += 1
            H2Ind[iArea] = NH2Area
            push!(H2Areas,H2Area(iArea,H2MaxDis,H2MaxRes,CTR.H2CompLoss,LStrategic))
        end
    end
    H2Data = H2System(NH2Area,H2Ind,H2Areas)
    return H2Data
end

function PrintMarketSummary(datapath,NArea,AMData,WPData,MCon,LineCap,CTI)
    msum = open(joinpath(datapath,"MarketSummary.dat"), "w")
    CutOffPrice = 100
    sWeek = 1

    Thermal = zeros(Float64,NArea)
    WindSolar = zeros(Float64,NArea)
    Exchange = zeros(Float64,NArea)
    DemandFlex = zeros(Float64,NArea)
    DemandFixed = zeros(Float64,NArea)

    @printf(msum,"%s %s %s %s %s %s \n","    ","     Thermal","   WindSolar","    Exchange","  DemandFlex"," DemandFixed")
    for iArea = 1:NArea
        AMD = AMData[iArea]
        for iMark = 1:AMD.NMStep 
            if AMD.MSData[iMark].Capacity[sWeek] > 0.0 && AMD.MSData[iMark].Price[sWeek] < CutOffPrice
                Thermal[iArea] += AMD.MSData[iMark].Capacity[sWeek]/CTI.NHoursWeek #GW
            elseif AMD.MSData[iMark].Capacity[sWeek] < 0.0
                DemandFlex[iArea] -= AMD.MSData[iMark].Capacity[sWeek]/CTI.NHoursWeek #GW
            end
        end
        WindSolar[iArea] = mean(WPData[iArea,:,sWeek,:])/CTI.DT #GW

        MCN = MCon[iArea]
        for iLine = 1:MCN.NCon
            Exchange[iArea] += 0.5*(LineCap[MCN.LIndxOut[iLine]]+LineCap[MCN.LIndxIn[iLine]])/CTI.NHoursWeek #GW
        end

        for iLoad = 1:AMD.NLoad
            DemandFixed[iArea] += mean(AMD.MLData[iLoad].Load[sWeek,:])/CTI.DT #GW
        end
        @printf(msum,"%4.0f %12.2f %12.2f %12.2f %12.2f %12.2f \n",iArea,Thermal[iArea],WindSolar[iArea],Exchange[iArea],DemandFlex[iArea],DemandFixed[iArea])
    end
    @printf(msum,"%s \n","----------------------------------------------------------------------")
    @printf(msum,"%s %12.2f %12.2f %12.2f %12.2f %12.2f \n"," SUM",sum(Thermal),sum(WindSolar),sum(Exchange),sum(DemandFlex),sum(DemandFixed))
    close(msum)

end

#Read grid data from file MASKENETT.DATA
function ReadMaske(dataset,NArea,CNS,CTR)
    filename = joinpath(dataset,"MASKENETT.DATA")
    fid = open(filename,"r")
    lineTmp = readline(fid) #skip header
    lineCnt = 0
    readMaske = true
    MaxLine = 1000
    Topo = zeros(Int,MaxLine,2)
    LineCap = zeros(Float64,0)
    LineLoss = zeros(Float64,0)
    while readMaske
        line = readline(fid); items = split(line,","); 
        if parse(Int,items[1]) == -1; break; end
        lineCnt += 1
        FromIndx = parse(Int, items[1]); FromArea = items[2];
        ToIndx = parse(Int, items[3]); ToArea = items[4];
        Topo[lineCnt,1] = FromIndx; Topo[lineCnt,2] = ToIndx;
        Topo[lineCnt+1,1] = ToIndx; Topo[lineCnt+1,2] = FromIndx;
        line = readline(fid); items = split(line,","); 
        Loss= parse(Float64,items[1])
        append!(LineLoss,Loss/100.0)
        append!(LineLoss,Loss/100.0) 
        line = readline(fid); items = split(line,","); 
        CapFrom = parse(Float64,items[2]); CapTo = parse(Float64,items[3])
        append!(LineCap,CTR.LineCapScale*CapFrom*CNS.MW2GWHWEEK) 
        append!(LineCap,CTR.LineCapScale*CapTo*CNS.MW2GWHWEEK)
        line = readline(fid); items = split(line,","); 
        NRev =  parse(Float64,items[1]);
        for iRev = 1:NRev
            line = readline(fid) 
        end
        lineCnt += 1
    end
    NLine = lineCnt

    MCon = []
    for iArea = 1:NArea
        LineIndxIn = []
        LineIndxOut = []
        for iLine = 1:NLine
            if Topo[iLine,1] == iArea
                push!(LineIndxOut,iLine)
            elseif Topo[iLine,2] == iArea
                push!(LineIndxIn,iLine)
            end
        end
        NL = length(LineIndxOut)
        push!(MCon,Connections(NL,LineIndxOut,LineIndxIn))
    end

    close(fid)
    println("Read MASKENETT.DATA")

    return MCon,LineCap,LineLoss,NLine
end

#Read time-dependent data for hydropower system
function ReadDynmod(dataset,AHData,NArea,AreaName,MyKeys,MaxModArea,NWeek)
    #Assumption: we need NWeek values
    MaMaxTD = zeros(Float64,NArea,MaxModArea,NWeek)
    MaMinTD = zeros(Float64,NArea,MaxModArea,NWeek)
    QMaxTD = zeros(Float64,NArea,MaxModArea,NWeek)
    QMinTD = zeros(Float64,NArea,MaxModArea,NWeek)
    QfoMaxTD = zeros(Float64,NArea,MaxModArea,NWeek)
    QfoMinTD = zeros(Float64,NArea,MaxModArea,NWeek)
    ModNrTD = zeros(Int,NArea,MaxModArea)
    for iArea = 1:NArea
        if AreaName[iArea] in MyKeys
            if iArea < 10
                fname = joinpath(dataset,"DYNMODELL.OM0$iArea")
            else
                fname = joinpath(dataset,"DYNMODELL.OM$iArea")
            end
            open(fname,"r") do fh
                
                I6 = zeros(Int32,6)
                for i=1:6
                    I6[i] = read(fh,Int32)
                end
                BLKL = I6[1]
                NModA = I6[3]
                JANT = I6[5]
                if JANT < NWeek 
                    error(println("JANT < NWeek"))
                end
                dp = 0
                MaMaxTmp = zeros(Float64,NModA,NWeek)
                MaMinTmp = zeros(Float64,NModA,NWeek)
                QMaxTmp = zeros(Float64,NModA,NWeek)
                QMinTmp = zeros(Float64,NModA,NWeek)
                QfoMaxTmp = zeros(Float64,NModA,NWeek)
                QfoMinTmp = zeros(Float64,NModA,NWeek)
                for juke = 1:NWeek
                    dp += BLKL
                    seek(fh,dp)
                    for iMod = 1:NModA
                        MaMaxTmp[iMod,juke] = Float64(read(fh,Float32))
                        MaMinTmp[iMod,juke] = Float64(read(fh,Float32))
                        tmp = Float64(read(fh,Float32))
                        tmp = Float64(read(fh,Float32))
                        QMaxTmp[iMod,juke] = Float64(read(fh,Float32))
                        QMinTmp[iMod,juke] = Float64(read(fh,Float32))
                        QfoMaxTmp[iMod,juke] = Float64(read(fh,Float32))
                        QfoMinTmp[iMod,juke] = Float64(read(fh,Float32))
                    end
                end
                seek(fh,BLKL*(JANT+1))
                for iMod = 1:NModA
                    ModNrTD[iArea,iMod] = Int(read(fh,Int32))
                end
                #They seem to be sorted similarly, could be skipped
                for iMod = 1:NModA
                    Indx = findfirst(x->x==AHData[iArea].MData[iMod].ModId,ModNrTD[iArea,1:NModA])
                    MaMaxTD[iArea,iMod,1:NWeek] = MaMaxTmp[Indx,1:NWeek]
                    MaMinTD[iArea,iMod,1:NWeek] = MaMinTmp[Indx,1:NWeek]
                    QMaxTD[iArea,iMod,1:NWeek] = QMaxTmp[Indx,1:NWeek]
                    QMinTD[iArea,iMod,1:NWeek] = QMinTmp[Indx,1:NWeek]
                    QfoMaxTD[iArea,iMod,1:NWeek] = QfoMaxTmp[Indx,1:NWeek]
                    QfoMinTD[iArea,iMod,1:NWeek] = QfoMinTmp[Indx,1:NWeek]
                end
            end # close DYNMODELL.OM
        end 
    end
    DMData = DynModData(MaMaxTD,MaMinTD,QMaxTD,QMinTD,QfoMaxTD,QfoMinTD)
    println("Read DYNMODELL.OM")
    
    return DMData
end
