function [gasTransportPores,inletVaporPressure,outletVaporPressure,vaporInletLinks,...
          vaporOutletLinks] = Condensation_GetBoundaryConditionsForDiffusion(network,...
                              options,invadedPore,equilibriumVaporPressure)
    
    gasTransportPores = find(not(invadedPore));
    
    inletVaporPressure = options.RelativeHumidityInlet*options.AirPressure; %TODO : Check this formula    
    vaporInletLinks=options.VaporInletLinks;
    
    %     Attention, fonction utilisée dans 2 contextes : nucléation et
    %     diffusionCondensation
    
    
    %Impose RH=1 on the invaded pores boundary links
    vaporOutletLinks=options.VaporOutletLinks;
    vaporOutletPores = network.GetPoresOfLink(vaporOutletLinks);
    vaporOutletPores = vaporOutletPores(:,2); % Take the neighboor pore inside the domain
    assert(length(vaporOutletPores)==length(vaporOutletLinks));
    outletVaporPressure = options.RelativeHumidityOutlet*options.AirPressure*equilibriumVaporPressure(vaporOutletPores);
    
    [boundaryLinks,~] = network.GetPoreRegionBoundaryLinks(find(invadedPore));
    boundaryPores = network.GetPoresOfLink(boundaryLinks);
    boundaryPores = boundaryPores(:,2); % Take one of the 2 neighboor pores
    RH=1;
    boundaryVaporPressure = RH*options.AirPressure*equilibriumVaporPressure(boundaryPores);
    
    vaporOutletLinks = [vaporOutletLinks,boundaryLinks];
    outletVaporPressure = [outletVaporPressure,boundaryVaporPressure];
    
    
    
end