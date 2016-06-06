function partialVaporPressure = Condensation_ComputePartialVaporPressure(network,gasTransportPores,diffusionConductances,inletVaporPressure,outletVaporPressure,vaporInletLinks,vaporOutletLinks,airPressure,temperature)
    
    
    %Convert inletVaporPressure,outletVaporPressure to concentration
    R = 8.314 ;
    pressureToConcentrationRatioPore = 1./(R*temperature);
    
    pressureToConcentrationRatioLink = network.InterpolatePoreDataToLink(pressureToConcentrationRatioPore);
    
    inletConcentration = inletVaporPressure.*pressureToConcentrationRatioLink(vaporInletLinks) ;   
    outletConcentration = outletVaporPressure.*pressureToConcentrationRatioLink(vaporOutletLinks);
    
    %Compute diffusion of vapor
    boundaryConditions=struct;
    boundaryConditions.inletLink = vaporInletLinks;
    boundaryConditions.outletLink = vaporOutletLinks;
    boundaryConditions.inletType = 'Dirichlet' ;
    boundaryConditions.outletType = 'Dirichlet' ;
    boundaryConditions.inletValue = inletConcentration;
    boundaryConditions.outletValue = outletConcentration;
    
    waterConcentration = ComputeLinearTransport(network,gasTransportPores,diffusionConductances,boundaryConditions);
    
    
    %Convert back concentrations to vapor pressure
    
    partialVaporPressure = waterConcentration./pressureToConcentrationRatioPore ;
    
end
