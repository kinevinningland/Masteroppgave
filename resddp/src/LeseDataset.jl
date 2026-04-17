
module LeseData
export disable_linecaps_by_area!, disable_prodcap_by_area!, disable_load_by_area!, disable_wind_by_area!, disable_trinn_by_area!, disable_h2caps_by_area!

"""
function disable_linecaps_by_area!: 

Nullstiller model.LineCap for linjer koblet til områder i `areas_off`.

mode:
  :both -> nuller både inn+ut (default)
  :to   -> nuller bare INN (LIndxIn)
  :from -> nuller bare UT (LIndxOut)
"""
function disable_linecaps_by_area!(model, areas_off::Vector{Int}; mode::Symbol = :both)
    @assert hasproperty(model, :MCon)    "model har ikke MCon"
    @assert hasproperty(model, :LineCap) "model har ikke LineCap"

    areaset = Set(areas_off)

    for a in areaset
        con = model.MCon[a]

        if mode === :both || mode === :from
            for idx in con.LIndxOut
                model.LineCap[idx] = 0.0
            end
        end

        if mode === :both || mode === :to
            for idx in con.LIndxIn
                model.LineCap[idx] = 0.0
            end
        end
    end

    return model
end

"""
Kobler ut vannkraftproduksjon i valgte områder (area-indekser).

Gjør to ting:
1) Setter AHData[a].MData[i].ProdCap = 0.0 for alle moduler i området a
2) Setter HSys[iSys].MaxProd = 0.0 og HSys[iSys].MinProd .= 0.0 for alle systemer som er knyttet til området via AreaSys

Merk:
- Model, ModuleData og HydroSystemData er immutable, så vi erstatter elementer i arrayene.
- Returnerer model (samme objekt, men med oppdaterte arrays).
"""
function disable_generation_by_area!(model, areas_off::Vector{Int})
    areaset = Set(areas_off)

    # --- 1) Null ut ProdCap i AHData (modulnivå) ---
    if hasproperty(model, :AHData)
        for a in areaset
            area_h = model.AHData[a]
            # area_h.MData er Array{ModuleData}
            oldM = area_h.MData
            newM = similar(oldM)
            MT = eltype(oldM)
            m_fields = fieldnames(MT)

            for i in eachindex(oldM)
                m = oldM[i]
                vals = map(f -> (f == :ProdCap ? 0.0 : getfield(m, f)), m_fields)
                newM[i] = MT(vals...)
            end

            # bygg ny AreaDataHydro med MData = newM
            AT = typeof(area_h)
            a_fields = fieldnames(AT)
            avals = map(f -> (f == :MData ? newM : getfield(area_h, f)), a_fields)
            model.AHData[a] = AT(avals...)
        end
    end

    # --- 2) Null ut aggregert produksjonskapasitet i HSys (det som brukes i stageprob/hprob) ---
    @assert hasproperty(model, :HSys) "model har ikke HSys"
    @assert hasproperty(model, :NAreaSys) "model har ikke NAreaSys"
    @assert hasproperty(model, :AreaSys) "model har ikke AreaSys"

    HS = model.HSys
    for a in areaset
        ns = model.NAreaSys[a]
        for j = 1:ns
            iSys = model.AreaSys[a, j]

            sys = HS[iSys]

            # MinProd er en Array inni immutable struct -> kan muteres direkte
            sys.MinProd .= 0.0

            # MaxProd er Float -> må erstatte hele structen
            ST = typeof(sys)
            s_fields = fieldnames(ST)
            s_vals = map(f -> (f == :MaxProd ? 0.0 : getfield(sys, f)), s_fields)
            HS[iSys] = ST(s_vals...)
        end
    end

    return model
end

"""
Setter all last = 0.0 i områdene som sendes inn.

areas_off: vektor med area-indekser (1..NArea)
"""
function disable_load_by_area!(model, areas_off::Vector{Int})

    @assert hasproperty(model, :AMData) "model har ikke AMData"

    areaset = Set(areas_off)

    for a in areaset
        area_m = model.AMData[a]

        # Hvert område kan ha flere dellaster
        for iLoad in 1:area_m.NLoad
            # Load er en Array → kan muteres direkte
            area_m.MLData[iLoad].Load .= 0.0
        end
    end

    return model
end

"""
Setter vindproduksjon (WPData) = 0.0 i områdene som sendes inn.

areas_off: vektor med area-indekser (1..NArea)
"""
function disable_wind_by_area!(model, areas_off::Vector{Int})

    @assert hasproperty(model, :WPData) "model har ikke WPData"

    areaset = Set(areas_off)

    for a in areaset
        # WPData dimensjoner:
        # (area, inflowYear, week, k)
        model.WPData[a, :, :, :] .= 0.0
    end

    return model
end

"""
Setter alle TRINN (markedssteg) til null i områdene som sendes inn.

Dette fjerner:
- Termisk produksjon
- Importkapasitet
- Fleksibelt forbruk

areas_off: vektor med area-indekser (1..NArea)
"""
function disable_trinn_by_area!(model, areas_off::Vector{Int})

    @assert hasproperty(model, :AMData) "model har ikke AMData"

    areaset = Set(areas_off)

    for a in areaset
        area_m = model.AMData[a]

        # Antall markedssteg i området
        for iStep in 1:area_m.NMStep
            # Capacity er Array{Float64,1}
            area_m.MSData[iStep].Capacity .= 0.0
        end
    end

    return model
end

function disable_h2caps_by_area!(model, areas_off::Vector{Int})
    @assert hasproperty(model, :H2Data) "model har ikke H2Data"

    H2 = model.H2Data
    if H2.NArea == 0
        return model
    end

    areaset = Set(areas_off)

    # H2Area er immutable -> må erstatte struct i H2.Areas
    HT = eltype(H2.Areas)
    h_fields = fieldnames(HT)

    for a in areaset
        # Map fra områdeindex -> H2-area index (0 hvis ikke H2-område)
        iH2 = H2.Ind[a]
        if iH2 > 0
            h = H2.Areas[iH2]
            vals = map(f -> (f == :MaxDis ? 0.0 : getfield(h, f)), h_fields)
            H2.Areas[iH2] = HT(vals...)
        end
    end

    return model
end

end # module