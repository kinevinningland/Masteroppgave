#Read extra data from file
using DelimitedFiles

function ReadDataLocation()
    f = open("../dataset.txt","r")
    dataset = readline(f) 
    return dataset
end


function ReadEndValue(dataset,NHSys)
    f = open(joinpath(dataset,"EndValueCuts.dat"),"r")
    line = readline(f)
    items = split(line," "); EVType = items[1];

    NEndCut = 0
    EndCutCoef = Array{Float64}(undef,1) 
    EndCutRHS = Array{Float64}(undef,1)
    NEndSeg = 0
    EndSegCoef = Array{Float64}(undef,1)

    if EVType == "CUT"
        line = readline(f)
        items = split(line," "); NEndCut = parse(Int, items[1]);
        println("Read end-value cuts")

        EndCutCoef = zeros(NEndCut,NHSys)
        EndCutRHS = zeros(NEndCut)
        for iCut = 1:NEndCut
            line = readline(f); items = split(line,";")
            for iArea = 1:NHSys
                EndCutCoef[iCut,iArea] = parse(Float64, items[iArea]);
            end
            EndCutRHS[iCut] = parse(Float64,items[NHSys+1]);
        end
    elseif EVType == "SEG"
        line = readline(f)
        items = split(line," "); NEndSeg = parse(Int, items[1]);
        println("Read end-value segments")

        EndSegCoef = zeros(NEndSeg)
        line = readline(f); items = split(line,";")
        for iSeg = 1:NEndSeg
            EndSegCoef[iSeg] = parse(Float64, items[iSeg]);
        end
    elseif EVType == "NONE"
        println("Zero end-valuation specified")
    else
        println("No end-valuation identified")
        exit()
    end

    close(f)

    return EndValuation(EVType,NEndCut,EndCutCoef,EndCutRHS,NEndSeg,EndSegCoef)
end


function ReadFeasCut(LFeasSpace,NHSys,dataset)
    CutDim = 7 #gammav, gammae, gammab, gammac, kappai, kappav, rhs
    NFeasCut = zeros(Int,NHSys)

    #Initial read to find dimensions
    if LFeasSpace
        for iSys = 1:NHSys
            try
                cutfilename = string(dataset,string(string("feas",iSys),".dat"))
                CutMatRaw = readdlm(cutfilename)
                NFeasCut[iSys] = size(CutMatRaw)[1]
                close(cutfilename)
            catch
                continue
            end
        end
    end
    
    MaxStatCut = maximum(NFeasCut)
    FCC = zeros(Float64,NHSys,MaxStatCut,CutDim)
    if LFeasSpace
        println("*Optimizing with Feasibility cuts*")
        for iSys = 1:NHSys
            try
                cutfilename = string(dataset,string(string("feas",iSys),".dat"))
                CutMatRaw = readdlm(cutfilename)
                dimCutMat = (NFeasCut[iSys],CutDim)
                CutMat = reshape(CutMatRaw,dimCutMat)
                for iCut = 1:NFeasCut[iSys] 
                    FCC[iSys,iCut,1:CutDim] = CutMat[iCut,1:CutDim]
                end
                println("...read ",NFeasCut[iSys]," feasibility cuts for area ",iSys)
                close(cutfilename)
            catch
                continue
            end
        end
    else
        println("*Optimizing without Feasibility cuts*")
    end
    return NFeasCut,FCC
end

function ReadDemandResponse(dataset,NArea,NWeek,AreaName,LDemandResponse)

    #Return dummy object if DR is not included
    if !LDemandResponse
        return DemandResponse(1, [0], Float64[], Float64[], [false], [0], 0)
    end

    #Hardcoded settings for "extra constraints" in load shifting DR formulation
    MaxLoadRecExtraConstr  = 2     #Max LoadRec where extra constraints are used
    ExtraConstrFilterScale = 2     #Filter size relative to LoadRec
    ExtraConstrSigma       = 0.5   #Sigma value for extra constraints

    f = open(string(dataset,"DRData.csv"),"r")

    #Reading parameter line
    line = readline(f)
    items = split(line,";")
    NAreaRead = parse(Int,items[1])
    if NAreaRead > NArea
        error(println("NArea in DRData.csv cannot be larger than NArea in params."))
    end
    NWeekRead = parse(Int,items[2])
    if NWeekRead > NWeek
        error(println("NWeek in DRData.csv cannot be larger than NWeek in params."))
    end
    NLoadRecStep = parse(Int,items[3])

    #Initializing arrays
    LoadRec    = zeros(Int, NLoadRecStep)
    MaxUpShift = zeros(Float64, NArea, NWeek, NLoadRecStep)
    MaxDnShift = zeros(Float64, NArea, NWeek, NLoadRecStep)
    LIncludeExtraConstr = fill(false, NLoadRecStep)
    ExtraConstrFilter   = zeros(Int, NLoadRecStep)

    for iStep = 1:NLoadRecStep
        #Reading LoadRec
        line = readline(f)
        items = split(line,";")
        LoadRec[iStep] = parse(Int,items[1])
        for iAreaInFile = 1:NAreaRead
            #Reading area name, finds the correct index in model
            line = readline(f)
            items = split(line,";")
            iArea = findfirst(x -> x == items[1], AreaName)
            if iArea === nothing
                error(println("Area name in DRData.csv: ", items[1], " does not exist in the model."))
            end
            #Reading MaxUpShift per week
            MaxUpShift[iArea, 1:NWeekRead, iStep] = [parse(Float64, val) for val in items[2:(NWeekRead+1)]]
            #Reading MaxDnShift per week
            line = readline(f)
            items = split(line,";")
            MaxDnShift[iArea, 1:NWeekRead, iStep] = [parse(Float64, val) for val in items[2:(NWeekRead+1)]]
        end
        #Setting "extra constraint" parameters for current LoadRec
        if LoadRec[iStep] <= MaxLoadRecExtraConstr
            LIncludeExtraConstr[iStep] = true
        end
        ExtraConstrFilter[iStep] = ExtraConstrFilterScale*LoadRec[iStep]
    end

    println("Read DRData.csv")

    return DemandResponse(NLoadRecStep, LoadRec, MaxUpShift, MaxDnShift, LIncludeExtraConstr, ExtraConstrFilter, ExtraConstrSigma)
end

function ReadCuts(NHSys,NStage,IM,dataset)
    f = open(string(dataset,"SDDPdims.txt"),"r")
    line = readline(f) 
    items = split(line," ")
    NScen = parse(Int,items[3])
    NCut = parse(Int,items[4])
    MaxIter = parse(Int,items[5])
    CutNoIter = parse(Int,items[6])
    close(f)

    NScenXMaxIter = NScen*MaxIter

    #Read cuts
    dimCCR = (NHSys,NStage,NScenXMaxIter)
    CCRfile = string(dataset,"CCR.dat")
    CCRraw = readdlm(CCRfile)
    CCR = reshape(CCRraw,dimCCR)
    #close(CCRfile)

    dimCCI = (IM.NSer,NStage,NScenXMaxIter) 
    CCIfile = string(dataset,"CCI.dat")
    CCIraw = readdlm(CCIfile)
    CCI = reshape(CCIraw,dimCCI)
    #close(CCIfile)

    dimCRHS = (NStage,NScenXMaxIter) 
    CRHSfile = string(dataset,"CRHS.dat")
    CRHSraw = readdlm(CRHSfile)
    CRHS = reshape(CRHSraw,dimCRHS)
    #close(CRHSfile)

    println("Completed reading ",NCut," cuts per ",NStage," stage")

    return NCut,CCR,CCI,CRHS
end
