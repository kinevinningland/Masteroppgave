using HDF5

function print_results_h5(dataset::String,RT::Result,model::Model,parameters::Parameters)

   file = h5open(string(dataset,"AggrSimResults_fixed_Mark.h5"),"w")
   
   attrs(file)["NArea"]  = model.NArea
   attrs(file)["NHSys"]  = model.NHSys
   attrs(file)["NScen"]  = parameters.Control.NScenSim
   attrs(file)["NStage"] = parameters.Control.NStageSim
   attrs(file)["NK"]     = parameters.Time.NK

   write(file, "ObjectiveValue", RT.ObjTable) #Added

   for iArea = 1:model.NArea
      areaGroup = create_group(file, model.AreaName[iArea])

      #Hydro results
      if iArea <= model.NHSys
         hydroGroup = create_group(areaGroup, "Hydro")

         write(hydroGroup, "Reservoir", RT.ReservoirTable[iArea,:,:,:])
         write(hydroGroup, "Production", RT.HProdTable[iArea,:,:,:])
         write(hydroGroup, "Ramping", repeat(RT.HRampTable[iArea,:,:], 1, 1, parameters.Time.NK))
         write(hydroGroup, "ReserveCapacity", repeat(RT.HCapTable[iArea,:,:,:], 1, 1, parameters.Time.NK))
         write(hydroGroup, "Spillage", RT.SpillTable[iArea,:,:,:])
         write(hydroGroup, "Inflow", parameters.Time.WeekFrac*repeat(RT.InflowTable[iArea,:,:,:], 1, 1, parameters.Time.NK))

         for dset in keys(hydroGroup)
            attrs(hydroGroup[dset])["Dim 1"] = "NScen"
            attrs(hydroGroup[dset])["Dim 2"] = "NStage"
            attrs(hydroGroup[dset])["Dim 3"] = "NK"
         end

         write(hydroGroup, "WaterValue", RT.WaterValueTable[iArea,:,:]) #Added
         attrs(hydroGroup["WaterValue"])["Dim 1"] = "NScen" #Added
         attrs(hydroGroup["WaterValue"])["Dim 2"] = "NStage" #Added

      end

      #Market results
      marketGroup = create_group(areaGroup, "Market")

      marketStepGroup = create_group(marketGroup, "Market_steps")
      #for iMark = 1:model.AMData[iArea].NMStep
      #   write(marketStepGroup, model.AMData[iArea].MSData[iMark].Name, RT.MarkTable[iArea,iMark,:,:,:])
      #end
      for iMark = 1:model.AMData[iArea].NMStep #byttet den over ut med denne

         base_name = model.AMData[iArea].MSData[iMark].Name
         name = base_name

         # Hvis navn finnes → legg til nummer
         if haskey(marketStepGroup, name)
            name = string(base_name, "_", iMark)
         end

         write(marketStepGroup, name, RT.MarkTable[iArea,iMark,:,:,:])
      end

      write(marketGroup, "Load", RT.LoadTable[iArea,:,:,:])
      write(marketGroup, "Wind", RT.WindTable[iArea,:,:,:])
      write(marketGroup, "Price", RT.PriceTable[iArea,:,:,:])
      write(marketGroup, "Rationing", RT.RationingTable[iArea,:,:,:])
      write(marketGroup, "DemandUpShift", RT.DemandUpTable[iArea,:,:,:])
      write(marketGroup, "DemandDownShift", RT.DemandDnTable[iArea,:,:,:])

      for dset in keys(marketStepGroup)
         attrs(marketStepGroup[dset])["Dim 1"] = "NScen"
         attrs(marketStepGroup[dset])["Dim 2"] = "NStage"
         attrs(marketStepGroup[dset])["Dim 3"] = "NK"
      end

      for dset in keys(marketGroup)
         attrs(marketGroup[dset])["Dim 1"] = "NScen"
         attrs(marketGroup[dset])["Dim 2"] = "NStage"
         attrs(marketGroup[dset])["Dim 3"] = "NK"
      end

      #H2 results
      if model.H2Data.Ind[iArea] > 0
          H2Group = create_group(areaGroup, "H2")
          write(H2Group, "Storage", RT.H2StoreTable[iArea,:,:,:])
          write(H2Group, "Discharge", RT.H2DisTable[iArea,:,:,:])

          for dset in keys(H2Group)
              attrs(H2Group[dset])["Dim 1"] = "NScen"
              attrs(H2Group[dset])["Dim 2"] = "NStage"
              attrs(H2Group[dset])["Dim 3"] = "NK"
          end
      end
      
   end

   #Flow results
   flowGroup = create_group(file, "Transmission")
   for iLine=1:model.NLine
      write(flowGroup, "Line "*string(iLine), RT.FlowTable[iLine,:,:,:])

      #Find "from" and "to" for line
      iFrom = -1
      iTo   = -1
      for iArea=1:model.NArea
         if iLine in model.MCon[iArea].LIndxOut
            iFrom = iArea
         elseif iLine in model.MCon[iArea].LIndxIn
            iTo = iArea
         end
      end

      attrs(flowGroup["Line "*string(iLine)])["From"]        = (iFrom != -1 ? model.AreaName[iFrom] : "Unknown")
      attrs(flowGroup["Line "*string(iLine)])["To"]          = (iTo   != -1 ? model.AreaName[iTo]   : "Unknown")
      attrs(flowGroup["Line "*string(iLine)])["Loss factor"] = string(model.LineLoss[iLine])
      attrs(flowGroup["Line "*string(iLine)])["Dim 1"]       = "NScen"
      attrs(flowGroup["Line "*string(iLine)])["Dim 2"]       = "NStage"
      attrs(flowGroup["Line "*string(iLine)])["Dim 3"]       = "NK"
   end

   if parameters.Control.LOperatingReserves #ADDED, mangler mark

      NZ = model.ORData.NZ
      price_zones = model.ORData.price_zones
      area_to_zone = model.ORData.area_to_zone
      areas_in_zone = model.ORData.areas_in_zone

      pzGroup = create_group(file, "XNordic_Reserve_req") #hovedmappe

      write(pzGroup, "Names", price_zones)
      write(pzGroup, "AreaToZone", area_to_zone)      
      write(pzGroup, "price_zones", price_zones)

      byZA = create_group(pzGroup, "ByZoneArea")

      for z in 1:NZ
         zname = price_zones[z]
         zGroup = create_group(byZA, zname) #Per sone
         if hasproperty(RT, :CapZoneUpTable)
            write(zGroup, "ReserveUp",   RT.CapZoneUpTable[z,:,:,:])
            write(zGroup, "ReserveDown", RT.CapZoneDownTable[z,:,:,:])
            write(zGroup, "ReserveDualUp",   RT.CapDualUpTable[z,:,:,:])
            write(zGroup, "ReserveDualDown", RT.CapDualDownTable[z,:,:,:])

            for dset in ["ReserveUp","ReserveDown"]
               attrs(zGroup[dset])["Dim 1"] = "NScen"
               attrs(zGroup[dset])["Dim 2"] = "NStage"
               attrs(zGroup[dset])["Dim 3"] = "NK"
            end
         end
         
         areasGroup = create_group(zGroup, "Areas")
         write(zGroup, "AreaIndices", areas_in_zone[z])
         attrs(zGroup["AreaIndices"])["Dim 1"] = "NAreaInZone"
         for a in areas_in_zone[z]
            aname = model.AreaName[a]
            aGroup = create_group(areasGroup, aname)

            for iSys in 1:model.NHSys
               if iSys == a
                  # dette hydrosystemet tilhører område a
                  if hasproperty(RT, :HydroCapUpTable)
                     write(aGroup, "HydroCapUp", RT.HydroCapUpTable[a,:,:,:])
                     write(aGroup, "HydroCapDown", RT.HydroCapDownTable[a,:,:,:])

                     for dset in ["HydroCapUp", "HydroCapDown"]
                        attrs(aGroup[dset])["Dim 1"] = "NScen"
                        attrs(aGroup[dset])["Dim 2"] = "NStage"
                        attrs(aGroup[dset])["Dim 3"] = "NK"
                     end
                  end
               end
            end
            if model.H2Data.Ind[a] > 0
               # dette H2-systemet tilhører område a
               write(aGroup, "H2CapUpDis", RT.H2CapUpDisTable[a,:,:,:])
               write(aGroup, "H2CapDownDis", RT.H2CapDownDisTable[a,:,:,:])
               write(aGroup, "H2CapUpChg", RT.H2CapUpChgTable[a,:,:,:])
               write(aGroup, "H2CapDownChg", RT.H2CapDownChgTable[a,:,:,:])

               for dset in ["H2CapUpDis", "H2CapDownDis", "H2CapUpChg", "H2CapDownChg"]
                  attrs(aGroup[dset])["Dim 1"] = "NScen"
                  attrs(aGroup[dset])["Dim 2"] = "NStage"
                  attrs(aGroup[dset])["Dim 3"] = "NK"
               end
               
            end

            if model.ORData.LMarkReserves
               write(aGroup, "MarketUpAreaPos",   RT.MarkCapUpTablePos[a,:,:,:])
               write(aGroup, "MarketDownAreaPos", RT.MarkCapDownTablePos[a,:,:,:])
               write(aGroup, "MarketUpAreaNeg",   RT.MarkCapUpTableNeg[a,:,:,:])
               write(aGroup, "MarketDownAreaNeg", RT.MarkCapDownTableNeg[a,:,:,:])

               for dset in ["MarketUpAreaPos", "MarketDownAreaPos", "MarketUpAreaNeg", "MarketDownAreaNeg"]
                     attrs(aGroup[dset])["Dim 1"] = "NScen"
                     attrs(aGroup[dset])["Dim 2"] = "NStage"
                     attrs(aGroup[dset])["Dim 3"] = "NK"
               end
            end

            # Bidrag per område (finnes allerede i RT)
            if hasproperty(RT, :WindCapDownTable)
               write(aGroup, "WindDownArea",   RT.WindCapDownTable[a,:,:,:])
               write(aGroup, "SpotPrice", RT.PriceTable[a,:,:,:])

               for dset in ["WindDownArea", "SpotPrice"]
                     attrs(aGroup[dset])["Dim 1"] = "NScen"
                     attrs(aGroup[dset])["Dim 2"] = "NStage"
                     attrs(aGroup[dset])["Dim 3"] = "NK"
               end
            end
         end
      end 

      

   end

   close(file)
end

function print_detailed_results_h5(dataset::String,DRT::DetailedResult,model::Model,parameters::Parameters)

   file = h5open(string(dataset,"DetSimResults_4Zone_wRes_wSharing2_100%.h5"),"w")
   
   attrs(file)["NArea"]  = model.NArea
   attrs(file)["NHSys"]  = model.NHSys
   attrs(file)["NScen"]  = parameters.Control.NScenSim
   attrs(file)["NStage"] = parameters.Control.NStageSim
   attrs(file)["NK"]     = parameters.Time.NK

   write(file, "ObjectiveValue", DRT.ObjTable) #Added

   for iArea = 1:model.NArea
      areaGroup = create_group(file, model.AreaName[iArea])
      if iArea <= model.NHSys
         #Hydro results
         hydroGroup = create_group(areaGroup, "Hydro")

         attrs(hydroGroup)["NMod"] = model.AHData[iArea].NMod

         for iMod=1:model.AHData[iArea].NMod
            moduleGroup = create_group(hydroGroup, "Module "*string(iMod))

            write(moduleGroup, "Reservoir", DRT.ReservoirTable[iArea,iMod,:,:,:])
            write(moduleGroup, "Production", DRT.HProdTable[iArea,iMod,:,:,:])
            write(moduleGroup, "Discharge", DRT.DischargeTable[iArea,iMod,:,:,:])
            write(moduleGroup, "Spillage", DRT.SpillTable[iArea,iMod,:,:,:])
            write(moduleGroup, "Bypass", DRT.BypassTable[iArea,iMod,:,:,:])
   
            for dset in keys(moduleGroup)
               attrs(moduleGroup[dset])["Dim 1"] = "NScen"
               attrs(moduleGroup[dset])["Dim 2"] = "NStage"
               attrs(moduleGroup[dset])["Dim 3"] = "NK"
            end
         end
         write(hydroGroup, "WaterValue", DRT.WaterValueTable[iArea,:,:]) #Added
         attrs(hydroGroup["WaterValue"])["Dim 1"] = "NScen" #Added
         attrs(hydroGroup["WaterValue"])["Dim 2"] = "NStage" #Added
      end

      #Market results
      marketGroup = create_group(areaGroup, "Market")

      marketStepGroup = create_group(marketGroup, "Market_steps")
      #for iMark = 1:model.AMData[iArea].NMStep #tatt bort
      #   write(marketStepGroup, model.AMData[iArea].MSData[iMark].Name, DRT.MarkTable[iArea,iMark,:,:,:])
      #end
      for iMark = 1:model.AMData[iArea].NMStep #ADDED

         base_name = model.AMData[iArea].MSData[iMark].Name
         name = base_name

         # Hvis navn finnes → legg til nummer
         if haskey(marketStepGroup, name)
            name = string(base_name, "_", iMark)
         end

         write(marketStepGroup, name, DRT.MarkTable[iArea,iMark,:,:,:])
      end

      write(marketGroup, "Load", DRT.LoadTable[iArea,:,:,:])
      write(marketGroup, "Wind", DRT.WindTable[iArea,:,:,:])
      write(marketGroup, "Price", DRT.PriceTable[iArea,:,:,:])
      write(marketGroup, "Rationing", DRT.RationingTable[iArea,:,:,:])
      write(marketGroup, "DemandUpShift", DRT.DemandUpTable[iArea,:,:,:])
      write(marketGroup, "DemandDownShift", DRT.DemandDnTable[iArea,:,:,:])

      for dset in keys(marketStepGroup)
         attrs(marketStepGroup[dset])["Dim 1"] = "NScen"
         attrs(marketStepGroup[dset])["Dim 2"] = "NStage"
         attrs(marketStepGroup[dset])["Dim 3"] = "NK"
      end

      for dset in keys(marketGroup)
         attrs(marketGroup[dset])["Dim 1"] = "NScen"
         attrs(marketGroup[dset])["Dim 2"] = "NStage"
         attrs(marketGroup[dset])["Dim 3"] = "NK"
      end
   end

   #Flow results
   flowGroup = create_group(file, "Transmission")
   for iLine=1:model.NLine
      write(flowGroup, "Line "*string(iLine), DRT.FlowTable[iLine,:,:,:])

      #Find "from" and "to" for line
      iFrom = -1
      iTo   = -1
      for iArea=1:model.NArea
         if iLine in model.MCon[iArea].LIndxOut
            iFrom = iArea
         elseif iLine in model.MCon[iArea].LIndxIn
            iTo = iArea
         end
      end

      attrs(flowGroup["Line "*string(iLine)])["From"]        = (iFrom != -1 ? model.AreaName[iFrom] : "Unknown")
      attrs(flowGroup["Line "*string(iLine)])["To"]          = (iTo   != -1 ? model.AreaName[iTo]   : "Unknown")
      attrs(flowGroup["Line "*string(iLine)])["Loss factor"] = string(model.LineLoss[iLine])
      attrs(flowGroup["Line "*string(iLine)])["Dim 1"]       = "NScen"
      attrs(flowGroup["Line "*string(iLine)])["Dim 2"]       = "NStage"
      attrs(flowGroup["Line "*string(iLine)])["Dim 3"]       = "NK"
   end

   if parameters.Control.LOperatingReserves #ADDED, mangler mark

      NZ = model.ORData.NZ
      price_zones = model.ORData.price_zones
      area_to_zone = model.ORData.area_to_zone
      areas_in_zone = model.ORData.areas_in_zone

      pzGroup = create_group(file, "XNordic_Reserve_req") #hovedmappe

      write(pzGroup, "Names", price_zones)
      write(pzGroup, "AreaToZone", area_to_zone)      
      write(pzGroup, "price_zones", price_zones)

      byZA = create_group(pzGroup, "ByZoneArea")

      for z in 1:NZ
         zname = price_zones[z]
         zGroup = create_group(byZA, zname) #Per sone
         if hasproperty(DRT, :CapZoneUpTable)
            write(zGroup, "ReserveUp",   DRT.CapZoneUpTable[z,:,:,:])
            write(zGroup, "ReserveDown", DRT.CapZoneDownTable[z,:,:,:])
            write(zGroup, "ReserveDualUp",   DRT.CapDualUpTable[z,:,:,:])
            write(zGroup, "ReserveDualDown", DRT.CapDualDownTable[z,:,:,:])

            for dset in ["ReserveUp","ReserveDown"]
               attrs(zGroup[dset])["Dim 1"] = "NScen"
               attrs(zGroup[dset])["Dim 2"] = "NStage"
               attrs(zGroup[dset])["Dim 3"] = "NK"
            end
         end
         write(zGroup, "SlackUp", DRT.SlackUpTable[z,:,:,:])
         write(zGroup, "SlackDown", DRT.SlackDownTable[z,:,:,:])
         
         areasGroup = create_group(zGroup, "Areas")
         write(zGroup, "AreaIndices", areas_in_zone[z])
         attrs(zGroup["AreaIndices"])["Dim 1"] = "NAreaInZone"
         for a in areas_in_zone[z]
            aname = model.AreaName[a]
            aGroup = create_group(areasGroup, aname)

            if a <= model.NHSys
               for iMod=1:model.AHData[a].NMod
                  moduleGroupRes = create_group(aGroup, "Module "*string(iMod))
                  write(moduleGroupRes, "HydroCapUp", DRT.HydroCapUpTable[a,iMod,:,:,:])
                  write(moduleGroupRes, "HydroCapDown", DRT.HydroCapDownTable[a,iMod,:,:,:])

                  for dset in keys(moduleGroupRes)
                     attrs(moduleGroupRes[dset])["Dim 1"] = "NScen"
                     attrs(moduleGroupRes[dset])["Dim 2"] = "NStage"
                     attrs(moduleGroupRes[dset])["Dim 3"] = "NK"
                  end
               end
                  #=
                  # dette hydrosystemet tilhører område a
                  if hasproperty(DRT, :HydroCapUpTable)
                     write(aGroup, "HydroCapUp", DRT.HydroCapUpTable[a,:,:,:])
                     write(aGroup, "HydroCapDown", DRT.HydroCapDownTable[a,:,:,:])

                     for dset in ["HydroCapUp", "HydroCapDown"]
                        attrs(aGroup[dset])["Dim 1"] = "NScen"
                        attrs(aGroup[dset])["Dim 2"] = "NStage"
                        attrs(aGroup[dset])["Dim 3"] = "NK"
                     end
                  end
                  =#
               
            end

            # Bidrag per område (finnes allerede i RT)
            if hasproperty(DRT, :WindCapDownTable)
               write(aGroup, "WindDownArea",   DRT.WindCapDownTable[a,:,:,:])
               write(aGroup, "SpotPrice", DRT.PriceTable[a,:,:,:])

               for dset in ["WindDownArea", "SpotPrice"]
                     attrs(aGroup[dset])["Dim 1"] = "NScen"
                     attrs(aGroup[dset])["Dim 2"] = "NStage"
                     attrs(aGroup[dset])["Dim 3"] = "NK"
               end
            end
            if model.ORData.LMarkReserves
               write(aGroup, "MarketUpAreaPos",   DRT.MarkCapUpTablePos[a,:,:,:])
               write(aGroup, "MarketDownAreaPos", DRT.MarkCapDownTablePos[a,:,:,:])
               write(aGroup, "MarketUpAreaNeg",   DRT.MarkCapUpTableNeg[a,:,:,:])
               write(aGroup, "MarketDownAreaNeg", DRT.MarkCapDownTableNeg[a,:,:,:])

               for dset in ["MarketUpAreaPos", "MarketDownAreaPos", "MarketUpAreaNeg", "MarketDownAreaNeg"]
                     attrs(aGroup[dset])["Dim 1"] = "NScen"
                     attrs(aGroup[dset])["Dim 2"] = "NStage"
                     attrs(aGroup[dset])["Dim 3"] = "NK"
               end
            end
         end
      end 
      if hasproperty(DRT, :SharingUpTable)
         sharingGroup = create_group(pzGroup, "Sharing")
         for z1 in 1:NZ
            for z2 in 1:NZ
               if z1 != z2
                  dname = string(price_zones[z1], "_to_", price_zones[z2])
                  write(sharingGroup, dname, DRT.SharingUpTable[z1,z2,:,:,:])
                  attrs(sharingGroup[dname])["From"]  = price_zones[z1]
                  attrs(sharingGroup[dname])["To"]    = price_zones[z2]
                  attrs(sharingGroup[dname])["Dim 1"] = "NScen"
                  attrs(sharingGroup[dname])["Dim 2"] = "NStage"
                  attrs(sharingGroup[dname])["Dim 3"] = "NK"
               end
            end
         end
      end
   end


   close(file)
end

function print_results(dataset::String,RT::Result,model::Model,parameters::Parameters)

   out = open(string(dataset,"hres.dat"),"w")
   @printf(out,"%6.0f %6.0f %6.0f %6.0f %8.2f \n",model.NArea,model.NHSys,model.NLine,parameters.Control.NStageSim,parameters.Time.DT)
   #dim (NMod,NScen,NStage), column-major order
   for iSys = 1:model.NHSys
      for iScen = 1:parameters.Control.NScenSim
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #Reservoir
               @printf(out,"%16.6f ",RT.ReservoirTable[iSys,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #hydro production
               @printf(out,"%16.6f ",RT.HProdTable[iSys,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #hydro ramping
               @printf(out,"%16.6f ",RT.HRampTable[iSys,iScen,iStage])
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #hydro reserve capacity
               @printf(out,"%16.6f ",RT.HCapTable[iSys,iScen,iStage])
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #Spillage
               @printf(out,"%16.6f ",RT.SpillTable[iSys,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #Inflow 
               @printf(out,"%16.6f ",parameters.Time.WeekFrac*RT.InflowTable[iSys,iScen,iStage])
            end
         end
         @printf(out,"%s \n","")
      end
   end
   @printf(out,"%s \n","")
   close(out)

   out = open(string(dataset,"mres.dat"),"w")
   @printf(out,"%6.0f %6.0f %6.0f %8.2f \n",model.NArea,model.NLine,parameters.Control.NStageSim,parameters.Time.DT)
   for iArea = 1:model.NArea
      for iScen = 1:parameters.Control.NScenSim
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #sum market
               val = sum(RT.MarkTable[iArea,iMark,iScen,iStage,k] for iMark=1:model.AMData[iArea].NMStep; init=0.0)#Added
               @printf(out, "%16.6f ", val) #Added
               #@printf(out,"%16.6f ",sum(RT.MarkTable[iArea,iMark,iScen,iStage,k] for iMark=1:model.AMData[iArea].NMStep))
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #Load
               @printf(out,"%16.6f ",RT.LoadTable[iArea,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #Wind
               @printf(out,"%16.6f ",RT.WindTable[iArea,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #Price 
               @printf(out,"%16.6f ",RT.PriceTable[iArea,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #Rationing
               @printf(out,"%16.6f ",RT.RationingTable[iArea,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #Demand Up 
               @printf(out,"%16.6f ",RT.DemandUpTable[iArea,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #Demand Down
               @printf(out,"%16.6f ",RT.DemandDnTable[iArea,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
      end 
   end
   @printf(out,"%s \n","")
   close(out)

   out = open(string(dataset,"flow.dat"),"w")
   for iLine = 1:model.NLine
      for iScen = 1:parameters.Control.NScenSim
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               @printf(out,"%16.6f ",RT.FlowTable[iLine,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
      end
   end
   @printf(out,"%s \n","")
   close(out)

end

function print_detailed_results(dataset::String,DRT::DetailedResult,model::Model,parameters::Parameters)

   out = open(string(dataset,"det_hres.dat"),"w")
   @printf(out,"%6.0f %6.0f %6.0f %6.0f %8.2f \n",model.NArea,model.NHSys,model.NLine,parameters.Control.NStageSim,parameters.Time.DT)
   for iSys = 1:model.NHSys
      @printf(out,"%6.0f ",model.AHData[iSys].NMod)
   end
   @printf(out,"%s \n","")
   #dim (NMod,NScen,NStage), column-major order
   for iSys = 1:model.NHSys
      for iMod = 1:model.AHData[iSys].NMod
         for iScen = 1:parameters.Control.NScenSim
            for iStage = 1:parameters.Control.NStageSim
               for k = 1:parameters.Time.NK
                  #Reservoir
                  @printf(out,"%16.6f ",DRT.ReservoirTable[iSys,iMod,iScen,iStage,k])
               end
            end
            @printf(out,"%s \n","")
            for iStage = 1:parameters.Control.NStageSim
               for k = 1:parameters.Time.NK
                  #hydro production
                  @printf(out,"%16.6f ",DRT.HProdTable[iSys,iMod,iScen,iStage,k])
               end
            end
            @printf(out,"%s \n","")
            for iStage = 1:parameters.Control.NStageSim
               for k = 1:parameters.Time.NK
                  #Discharge
                  @printf(out,"%16.6f ",DRT.DischargeTable[iSys,iMod,iScen,iStage,k])
               end
            end
            @printf(out,"%s \n","")
            for iStage = 1:parameters.Control.NStageSim
               for k = 1:parameters.Time.NK
                  #Spillage
                  @printf(out,"%16.6f ",DRT.SpillTable[iSys,iMod,iScen,iStage,k])
               end
            end
            @printf(out,"%s \n","")
            for iStage = 1:parameters.Control.NStageSim
               for k = 1:parameters.Time.NK
                  #Bypass
                  @printf(out,"%16.6f ",DRT.BypassTable[iSys,iMod,iScen,iStage,k])
               end
            end
            @printf(out,"%s \n","")
         end
      end
   end
   @printf(out,"%s \n","")
   close(out)

   out = open(string(dataset,"det_mres.dat"),"w")
   @printf(out,"%6.0f %6.0f %6.0f %8.2f \n",model.NArea,model.NLine,parameters.Control.NStageSim,parameters.Time.DT)
   for iArea = 1:model.NArea
      for iScen = 1:parameters.Control.NScenSim
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #sum market
               #@printf(out,"%16.6f ",sum(DRT.MarkTable[iArea,iMark,iScen,iStage,k] for iMark=1:model.AMData[iArea].NMStep))
               val = sum((DRT.MarkTable[iArea,iMark,iScen,iStage,k] for iMark=1:model.AMData[iArea].NMStep); init=0.0)
               @printf(out,"%16.6f ", val)
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #Load
               @printf(out,"%16.6f ",DRT.LoadTable[iArea,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #Wind
               @printf(out,"%12.2f ",DRT.WindTable[iArea,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #Price 
               @printf(out,"%16.6f ",DRT.PriceTable[iArea,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #Rationing
               @printf(out,"%16.6f ",DRT.RationingTable[iArea,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #Demand Up 
               @printf(out,"%16.6f ",DRT.DemandUpTable[iArea,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               #Demand Dn
               @printf(out,"%16.6f ",DRT.DemandDnTable[iArea,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
      end 
   end
   @printf(out,"%s \n","")
   close(out)

   out = open(string(dataset,"det_flow.dat"),"w")
   for iLine = 1:model.NLine
      for iScen = 1:parameters.Control.NScenSim
         for iStage = 1:parameters.Control.NStageSim
            for k = 1:parameters.Time.NK
               @printf(out,"%16.6f ",DRT.FlowTable[iLine,iScen,iStage,k])
            end
         end
         @printf(out,"%s \n","")
      end
   end
   @printf(out,"%s \n","")
   close(out)
end


function print_dims(datapath::String,NHSys::Int,NH2Area::Int,NStage::Int,NScen::Int,NCut::Int,MaxIter::Int,CCMaxIter::Int)
   out = open(string(datapath,"SDDPdims.txt"),"w")
   @printf(out,"%.0f ",NHSys)
   @printf(out,"%.0f ",NH2Area)
   @printf(out,"%.0f ",NStage)
   @printf(out,"%.0f ",NScen)
   @printf(out,"%.0f ",NCut)
   @printf(out,"%.0f ",MaxIter)
   @printf(out,"%.0f ",CCMaxIter)
   @printf(out,"%s \n"," ")
   close(out)
end

function print_strategy(datapath::String,strategy::Strategy,LCostApprox::Bool)
   
   open(string(datapath,"CCR.dat"),"w") do io
      writedlm(io,strategy.CCR)
   end
   open(string(datapath,"CCH.dat"),"w") do io
      writedlm(io,strategy.CCH)
   end
   open(string(datapath,"CCI.dat"),"w") do io
      writedlm(io,strategy.CCI)
   end
   open(string(datapath,"CRHS.dat"),"w") do io
      writedlm(io,strategy.CRHS)
   end

   if LCostApprox
      open(string(datapath,"CFP.dat"),"w") do io
         writedlm(io,strategy.CFP)
      end
      open(string(datapath,"CFR.dat"),"w") do io
         writedlm(io,strategy.CFR)
      end
      open(string(datapath,"CFC.dat"),"w") do io
         writedlm(io,strategy.CFC)
      end
      open(string(datapath,"CFHD.dat"),"w") do io
         writedlm(io,strategy.CFHD)
      end
      open(string(datapath,"CFHS.dat"),"w") do io
         writedlm(io,strategy.CFHS)
      end
      open(string(datapath,"CFRHS.dat"),"w") do io
         writedlm(io,strategy.CFRHS)
      end
   end
end

function print_feas(datapath::String,feas_space::FeasibilitySpace,NHSys::Int)
   for iSys = 1:NHSys
      fcout = open(string(datapath,string(string("feas",iSys),".dat")),"w")
      for c = 1:feas_space.NFeasCut[iSys]
         @printf(fcout,"%12.4f %s %12.4f %s %12.4f %s %12.4f %s %12.4f %s %12.4f %s %18.4f %s %4.1f \n",feas_space.FCC[iSys,c,1],",",feas_space.FCC[iSys,c,2],",",
                 feas_space.FCC[iSys,c,3],",",feas_space.FCC[iSys,c,4],",",feas_space.FCC[iSys,c,5],",",feas_space.FCC[iSys,c,6],",",feas_space.FCC[iSys,c,7],",",c)
      end
      close(fcout)
   end
end

