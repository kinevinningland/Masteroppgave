using LinearAlgebra
using StatsBase
using Printf
using Random

function ReadInflow(dataset,NHSys,NWeek,NScen,NResid,LIgnoreCrossCorr)
    #Read inflow statistics all data in GWh
    filename = joinpath(dataset,"InflowModelLOG.dat")
    f = open(filename,"r")
    line = readline(f); items = split(line," ");
    NWeekRead = parse(Int,items[1]);
    if (NWeekRead != NWeek); error("NWeek != NWeek in inflow model"); end
    NResIdRead = parse(Int,items[2]);
    if (NResIdRead != NResid); error("NResIdRead != NResid in inflow model"); end
    NSysRead = parse(Int,items[3]);
    if (NSysRead != NHSys); error("NHSys != NHSys in inflow model --> Run inflow script!"); end

    InflowAv = zeros(Float64,NHSys)
    InflowMean = zeros(Float64,NHSys,NWeek)
    InflowSDev = zeros(Float64,NHSys,NWeek)
    Resid = zeros(Float64,NHSys,NWeek,NResid)
    line = readline(f)
    for iSys = 1:NHSys
        line = readline(f) #Skip comment 
        for iWeek = 1:NWeek
            line = readline(f); items = split(line," "); 
            InflowMean[iSys,iWeek] = parse(Float64,items[2]);
            InflowAv[iSys] += InflowMean[iSys,iWeek]
            InflowSDev[iSys,iWeek] = parse(Float64,items[3]);
        end
        for iWeek = 1:NWeek
            line = readline(f); items = split(line," "); 
            for iYear = 1:NResid
                Resid[iSys,iWeek,iYear] = parse(Float64,items[iYear])
            end
        end
    end
    
    CorrMat = zeros(Float64,NHSys,NHSys)
    line = readline(f)
    for iSys = 1:NHSys
        line = readline(f); items = split(line," "); 
        for jSer in 1:NHSys
            if LIgnoreCrossCorr
                if jSer == iSys
                    CorrMat[iSys,jSer] = parse(Float64,items[jSer])
                end
            else
                CorrMat[iSys,jSer] = parse(Float64,items[jSer])
            end
        end
    end
    close(f)
    NSer = NHSys
    println("Read inflow model parameters")

    PossibleSamples = collect(1:NResid)

    return InflowModel(NSer,NResid,InflowAv,InflowMean,InflowSDev,Resid,CorrMat)
end


function SampleScenario(NScen,NStage,NWeek,IM,LExtreme; fixed_seed = false)

    if fixed_seed
        Random.seed!(12345)
    end

    SScen = zeros(Int,NScen,NStage)
    Zscen = zeros(Float64,NScen,NStage,IM.NSer)
    PossibleSamples = collect(1:IM.NResid)
    for iScen = 1:NScen
        SScen[iScen,1:NStage] = sample(PossibleSamples,NStage,replace=true)
        SScen[iScen,1] = Int(ceil(IM.NResid/2)) #Assumption: fixed starting point
        for iStage = 1:NStage
            iWeek = mod1(iStage,NWeek)
            mySample = SScen[iScen,iStage]
            CorrZ = zeros(Float64,IM.NSer)
            if iStage > 1
                CorrZ[1:IM.NSer] = IM.CorrMat*Zscen[iScen,iStage-1,1:IM.NSer] #(NSerXNSer)*(NSerX1)
            end
            for iSer = 1:IM.NSer
                Zscen[iScen,iStage,iSer] = CorrZ[iSer]+IM.Resid[iSer,iWeek,mySample]
            end
        end 
    end
    if LExtreme && NScen > 2
        for iStage = 1:NStage
            SScen[1,iStage] = 1
            SScen[NScen,iStage] = IM.NResid
        end
    end

    LCheck = false
    if LCheck
        GWH = zeros(Float64,NStage,NScen,IM.NSer)
        SumGWH = zeros(Float64,NScen,IM.NSer)
        NSimYear = NStage/52.0
        for iScen = 1:NScen
           for iStage = 1:NStage
               iWeek = mod1(iStage,NWeek)
               for iSer = 1:IM.NSer
                   GWH[iStage,iScen,iSer] = Zscen[iScen,iStage,iSer]*IM.InflowSDev[iSer,iWeek]+IM.InflowMean[iSer,iWeek]
                   SumGWH[iScen,iSer] += GWH[iStage,iScen,iSer]
               end
           end
        end
        SumGWH = SumGWH/NSimYear
        for iSer = 1:IM.NSer
            MinS = minimum(GWH[1:NStage,1:NScen,iSer])
            MaxS = maximum(GWH[1:NStage,1:NScen,iSer])
            Min = minimum(SumGWH[1:NScen,iSer])
            Mean = mean(SumGWH[1:NScen,iSer])
            Max = maximum(SumGWH[1:NScen,iSer])
            @printf("%s %3.0f %s %12.0f %12.0f %12.0f %12.0f %12.0f \n","MM3 @",iSer,"=",MinS,Min,Mean,Max,MaxS)
        end
    end

    return SampleScen(SScen,Zscen)
end
