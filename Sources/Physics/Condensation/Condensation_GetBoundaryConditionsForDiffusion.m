function [gasTransportPores,inletVaporPressure,outletVaporPressure,vaporInletLinks,...
          vaporOutletLinks] = Condensation_GetBoundaryConditionsForDiffusion(network,...
                              options,invadedPore,equilibriumVaporPressure)
    
                          
    gasTransportPores = intersect(options.VaporTransportPores, find(not(invadedPore)));
    
    equilibriumVaporPressureLinks = network.InterpolatePoreDataToLink(equilibriumVaporPressure);  
    
    vaporInletLinks = options.VaporInletLinks;
    inletVaporPressure = options.RelativeHumidityInlet*equilibriumVaporPressureLinks(vaporInletLinks); 
    
    %     Attention, fonction utilisée dans 2 contextes : nucléation et
    %     diffusionCondensation
    
    
    %Impose RH=1 on the invaded pores boundary links
    vaporOutletLinks = options.VaporOutletLinks;
    outletPores = network.GetPoresOfLink(vaporOutletLinks);
    outletLinksToKeep = not(invadedPore(outletPores(:,2))); % remove outlet pores wich are invaded with water
    vaporOutletLinks = vaporOutletLinks(outletLinksToKeep);
    outletVaporPressure = options.RelativeHumidityOutlet*equilibriumVaporPressureLinks(vaporOutletLinks(outletLinksToKeep));
    
    [boundaryLinks,~] = network.GetPoreRegionBoundaryLinks(find(invadedPore));
    boundaryLinks = setdiff(boundaryLinks,network.GetLinksFrontiere(1:network.GetNumberOfBoundaries));
    RH = 1;
    boundaryVaporPressure = RH*equilibriumVaporPressureLinks(boundaryLinks);
    
    vaporOutletLinks = [vaporOutletLinks,boundaryLinks];
    outletVaporPressure = [outletVaporPressure;boundaryVaporPressure];
    
    
end