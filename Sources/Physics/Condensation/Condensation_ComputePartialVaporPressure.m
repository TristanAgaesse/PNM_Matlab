function partialVaporPressure = Condensation_ComputePartialVaporPressure(network,gasTransportPores,diffusionConductances,inletVaporPressure,outletVaporPressure,vaporInletLinks,vaporOutletLinks,airPressure,temperature)
    
    
    %Convert inletVaporPressure,outletVaporPressure to concentration
    R = 8.314 ;
    airConcentration = airPressure./(R*temperature);
    
    inletConcentration = inletVaporPressure*airConcentration(vaporInletLinks) ;
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
