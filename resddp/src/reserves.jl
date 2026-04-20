module Reserves

function ReadOperatingReserves(NArea, NHsys, NAreaSys, H2Data, AMData)::RRData
    LH2Reserves = false
    LMarkReserves = false

    price_zones = ["NO1", "NO2", "NO3", "NO4", "NO5","Others"]
    NZ = length(price_zones)
    zone_reqs = [
        ReserveZoneReq("NO1", 0.344, 0.172, 0.58, 0.56),
        ReserveZoneReq("NO2", 1.4, 1.4, 0.37, 0.37),
        ReserveZoneReq("NO3", 0.29, 0.145, 0.33, 0.38),
        ReserveZoneReq("NO4", 0.35, 0.175, 0.31, 0.33),
        ReserveZoneReq("NO5", 1.4, 1.4, 0.37, 0.37),
        ReserveZoneReq("Others", 0.0, 0.0, 0.0, 0.0)
    ]

    area_to_zone = fill(findfirst(==("Others"), price_zones),NArea) #Areas not explicitly listed default to "Others"
    area_to_zone[33] = findfirst(==("NO3"), price_zones) #OK
    area_to_zone[34] = findfirst(==("NO5"), price_zones) #OK
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

    area_to_zone[7] = findfirst(==("NO5"), price_zones) #OK

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

    pos_by_area = Dict{Int, Set{Int}}()
    neg_by_area = Dict{Int, Set{Int}}()

    return RRData(NZ,NZ-1,price_zones,zone_reqs,area_to_zone,areas_in_zone,hydrosys_to_area,h2_to_area,LH2Reserves,LMarkReserves,pos_by_area,neg_by_area)

end

end # module
