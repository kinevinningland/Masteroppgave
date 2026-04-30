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

#=
function ReadOperatingReserves(NArea, NHSys, NAreaSys, AreaSys, H2Data, AMData,AreaName,LOperatingReserves,MCon)
    #Return dummy object if OR is not included
    if !LOperatingReserves
        return OperatingReserves(0,String[],ReserveZoneReq[],Vector{Vector{Int}}(),Int[],false,0,0)
    end
    LMarkReserves = false
    LZoneReq = true

    if LZoneReq
        price_zones = ["NO1", "NO2", "NO3", "NO4"]
        NZ = length(price_zones)

        area_to_zone = fill(0,NArea) 
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
    else
        price_zones = ["NO"]
        NZ = length(price_zones)

        area_to_zone = fill(0,NArea) 
        area_to_zone[33] = findfirst(==("NO"), price_zones) 
        area_to_zone[34] = findfirst(==("NO"), price_zones) 
        area_to_zone[35] = findfirst(==("NO"), price_zones) 
        area_to_zone[36] = findfirst(==("NO"), price_zones) 
        area_to_zone[66] = findfirst(==("NO"), price_zones) 
        area_to_zone[67] = findfirst(==("NO"), price_zones) 
        area_to_zone[68] = findfirst(==("NO"), price_zones) 
        area_to_zone[73] = findfirst(==("NO"), price_zones) 
        area_to_zone[74] = findfirst(==("NO"), price_zones) 
        area_to_zone[75] = findfirst(==("NO"), price_zones) 
        area_to_zone[1] = findfirst(==("NO"), price_zones) 
        area_to_zone[2] = findfirst(==("NO"), price_zones) 
        area_to_zone[3] = findfirst(==("NO"), price_zones) 
        area_to_zone[4] = findfirst(==("NO"), price_zones) 
        area_to_zone[5] = findfirst(==("NO"), price_zones) 
        area_to_zone[6] = findfirst(==("NO"), price_zones) 
        area_to_zone[8] = findfirst(==("NO"), price_zones) 
        area_to_zone[9] = findfirst(==("NO"), price_zones) 
        area_to_zone[10] = findfirst(==("NO"), price_zones) 
        area_to_zone[11] = findfirst(==("NO"), price_zones) 
        area_to_zone[7] = findfirst(==("NO"), price_zones) 

    end

    areas_in_zone = [Int[] for _ in 1:NZ]
    for a in 1:NArea
        z = area_to_zone[a]
        if z > 0
            push!(areas_in_zone[z], a)
        end
    end

    hydrosys_to_area = fill(0, NHSys)
    for a in 1:NArea
        for j in 1:NAreaSys[a]
            hydrosys_to_area[AreaSys[a,j]] = a
        end
    end

    max_load_per_zone = zeros(Float64,NZ)
    for z in 1:NZ
        for a in areas_in_zone[z]
            for iLoad in 1:AMData[a].NLoad
                max_load_per_zone[z] += maximum(AMData[a].MLData[iLoad].Load)
            end
        end
    end
    
    owp_areas_in_zone = [Int[] for _ in 1:NZ]
    for z in 1:NZ
        for iArea in areas_in_zone[z]
            if endswith(AreaName[iArea], "-OWP")
                push!(owp_areas_in_zone[z], iArea)
            end
        end
    end
    
    if LZoneReq
        zone_reqs = [
            ReserveZoneReq("NO1", 0.344, 0.172, 0.58, 0.56, 0.48, 0.4,max_load_per_zone[1],owp_areas_in_zone[1]),
            ReserveZoneReq("NO2", 1.4, 1.4, 0.37, 0.37, 0.48, 0.4,max_load_per_zone[2],owp_areas_in_zone[2]),
            ReserveZoneReq("NO3", 0.29, 0.145, 0.33, 0.38, 0.48, 0.4,max_load_per_zone[3],owp_areas_in_zone[3]),
            ReserveZoneReq("NO4", 0.35, 0.175, 0.31, 0.33, 0.48, 0.4,max_load_per_zone[4],owp_areas_in_zone[4]),
            #ReserveZoneReq("NO5", 1.4, 1.4, 0.37, 0.37),
        ]
    else
        zone_reqs = [
            ReserveZoneReq("NO",1.4,1.4,0.17,0.17,0.48,0.4,max_load_per_zone[1]+max_load_per_zone[2]+max_load_per_zone[3]+max_load_per_zone[4],owp_areas_in_zone[1]),
        ]
    end
    
    a = 10
    b = 150

    println("Read ORData.csv")

    return OperatingReserves(NZ,price_zones,zone_reqs,areas_in_zone,hydrosys_to_area,LMarkReserves,a,b)
end
=#

function ReadOperatingReserves(dataset,NArea, NHSys, NAreaSys, AreaSys, AMData,AreaName,LOperatingReserves)
    #Return dummy object if OR is not included
    if !LOperatingReserves
        return OperatingReserves(0,String[],ReserveZoneReq[],Vector{Int}(),Vector{Vector{Int}}(),Int[],false,0,0)
    end
    LMarkReserves = false
    LZoneReq = true

    zone_reqs = Dict{String, Vector{Float64}}()
    price_zones = String[]

    # Read zone requirements
    filename = LZoneReq ? "ORData_zone_reqs.csv" : "ORData_zone_reqs_system.csv"
    f = open(joinpath(dataset, filename), "r")
    readline(f)
    for line in eachline(f)
        items = split(line, ",")
        zone = strip(items[1])
        push!(price_zones, zone)
        zone_reqs[zone] = parse.(Float64, items[2:end])
    end
    close(f)
    NZ = length(price_zones)

    # Read area-to-zone mapping
    area_to_zone = fill(0, NArea)
    filename = LZoneReq ? "ORData_area_zones.csv" : "ORData_area_system.csv"
    f = open(joinpath(dataset, filename), "r")
    readline(f)
    for line in eachline(f)
        items = split(line, ",")
        iArea = parse(Int, items[1])
        zone = strip(items[2])
        z = findfirst(==(zone), price_zones)
        if !isnothing(z)
            area_to_zone[iArea] = z
        end
    end
    
    
    areas_in_zone = [Int[] for _ in 1:NZ]
    for a in 1:NArea
        z = area_to_zone[a]
        if z > 0
            push!(areas_in_zone[z], a)
        end
    end

    hydrosys_to_area = fill(0, NHSys)
    for a in 1:NArea
        for j in 1:NAreaSys[a]
            hydrosys_to_area[AreaSys[a,j]] = a
        end
    end

    max_load_per_zone = zeros(Float64,NZ)
    for z in 1:NZ
        for a in areas_in_zone[z]
            for iLoad in 1:AMData[a].NLoad
                max_load_per_zone[z] += maximum(AMData[a].MLData[iLoad].Load)
            end
        end
    end
    
    owp_areas_in_zone = [Int[] for _ in 1:NZ]
    for z in 1:NZ
        for iArea in areas_in_zone[z]
            if endswith(AreaName[iArea], "-OWP")
                push!(owp_areas_in_zone[z], iArea)
            end
        end
    end
    
    zone_reqs = [
        ReserveZoneReq(
            price_zones[z],
            zone_reqs[price_zones[z]][1], # RI_up
            zone_reqs[price_zones[z]][2], # RI_down
            zone_reqs[price_zones[z]][3], # NI_up
            zone_reqs[price_zones[z]][4], # NI_down
            zone_reqs[price_zones[z]][5], # NI_up_OWP
            zone_reqs[price_zones[z]][6], # NI_down_OWP
            max_load_per_zone[z],
            owp_areas_in_zone[z]
        ) for z in 1:NZ
    ]

    a = 10
    b = 150

    println("Read ORData.csv")

    return OperatingReserves(NZ,price_zones,zone_reqs,area_to_zone,areas_in_zone,hydrosys_to_area,LMarkReserves,a,b)
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
