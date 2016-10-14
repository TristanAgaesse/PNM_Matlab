function [partialVaporPressure,waterConcentration,diffusionFlux,effectiveDiffusion ] = Condensation_ComputePartialVaporPressure(network,...
                                    gasTransportPores,diffusionConductances,inletVaporPressure,outletVaporPressure,...
                                    vaporInletLinks,vaporOutletLinks,airPressure,temperature)
    
    
    %Convert inletVaporPressure,outletVaporPressure to concentration
    R = 8.314 ;
    pressureToConcentrationRatioPore = 1./(R*temperature);
    
    pressureToConcentrationRatioLink = network.InterpolatePoreDataToLinkInAPoreRegion(pressureToConcentrationRatioPore,gasTransportPores);
    
    inletConcentration = inletVaporPressure.*pressureToConcentrationRatioLink(vaporInletLinks) ;   
    outletConcentration = outletVaporPressure.*pressureToConcentrationRatioLink(vaporOutletLinks);
    
    % Attention si lien inlet bouché par l’eau : reporter le flux sur autres liens.
    
    %Compute diffusion of vapor
    boundaryConditions=struct;
    boundaryConditions.inletLink = vaporInletLinks;
    boundaryConditions.outletLink = vaporOutletLinks;
    boundaryConditions.inletType = 'Dirichlet' ;
    boundaryConditions.outletType = 'Dirichlet' ;
    boundaryConditions.inletValue = inletConcentration;
    boundaryConditions.outletValue = outletConcentration;
    boundaryConditions.solver='mldivide';
    
    [waterConcentration,diffusionFlux,effectiveDiffusion ] = ComputeLinearTransport(network,gasTransportPores,diffusionConductances,boundaryConditions);
    
    
    %Convert back concentrations to vapor pressure
    
    partialVaporPressure = waterConcentration./pressureToConcentrationRatioPore ;
    
end
