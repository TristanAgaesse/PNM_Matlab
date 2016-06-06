function partialVaporPressure = Condensation_ComputePartialVaporPressure(network,gasTransportPores,diffusionConductances,inletVaporPressure,outletVaporPressure,vaporInletLinks,vaporOutletLinks,airPressure,temperature)
    
    
    %Convert inletVaporPressure,outletVaporPressure to concentration
    R = 8.314 ;
    airConcentration = airPressure./(R*temperature);
    
    linkData = network.InterpolatePoreDataToLink(poreData)
    
    inletConcentration = inletVaporPressure*airConcentration(vaporInletLinks) ;       % Approximate air concentration on inner links
    outletConcentration = outletVaporPressure*airConcentration(vaporOutletLinks);
    
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
    
    partialVaporPressure = waterConcentration./airConcentration ;
    
end
