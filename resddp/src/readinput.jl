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

function ReadOperatingReserves(NArea, NHSys, NAreaSys, AreaSys, H2Data, AMData,AreaName,LOperatingReserves,MCon)
    #Return dummy object if OR is not included
    if !LOperatingReserves
        return OperatingReserves(0,0,String[],ReserveZoneReq[],Int[],Vector{Vector{Int}}(),Int[],Int[],false,false,Dict{Int, Set{Int}}(),Dict{Int, Set{Int}}(),false,Set{Tuple{Int, Int}}())
    end

    LH2Reserves = false
    LMarkReserves = false
    LSharing = true

    price_zones = ["NO1", "NO2", "NO3", "NO4"]
    NZ = length(price_zones)
    zone_reqs = [
        ReserveZoneReq("NO1", 0.344, 0.172, 0.58, 0.56),
        ReserveZoneReq("NO2", 1.4, 1.4, 0.37, 0.37),
        ReserveZoneReq("NO3", 0.29, 0.145, 0.33, 0.38),
        ReserveZoneReq("NO4", 0.35, 0.175, 0.31, 0.33),
        #ReserveZoneReq("NO5", 1.4, 1.4, 0.37, 0.37),
    ]

    area_to_zone = fill(findfirst(==("Others"), price_zones),NArea) #Areas not explicitly listed default to "Others"
    area_to_zone[33] = findfirst(==("NO3"), price_zones) #OK
    area_to_zone[34] = findfirst(==("NO2"), price_zones) #OK Var NO5
    area_to_zone[35] = findfirst(==("NO2"), price_zones) #OK
    area_to_zone[36] = findfirst(==("NO2"), price_zones) #OK
    area_to_zone[66] = findfirst(==("NO3"), price_zones) # H2-M?
    area_to_zone[67] = findfirst(==("NO2"), price_zones) #H2-S?
    area_to_zone[68] = findfirst(==("NO4"), price_zones) #H2-N?
    area_to_zone[73] = findfirst(==("NO4"), price_zones) #OK
    area_to_zone[74] = findfirst(==("NO4"), price_zones) #OK
    area_to_zone[75] = findfirst(==("NO4"), price_zones) #OK
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

    neighboring_zones = Set{Tuple{Int, Int}}()#legge inn i ORData?
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

    println("Read ORData.csv")

    return OperatingReserves(NZ,NZ-1,price_zones,zone_reqs,area_to_zone,areas_in_zone,hydrosys_to_area,h2_to_area,LH2Reserves,LMarkReserves,pos_by_area,neg_by_area,LSharing,neighboring_zones)
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
