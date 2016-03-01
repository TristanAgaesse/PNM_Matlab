function  [temperature,heatTransferCoefficient] = Condensation_ComputeTemperatureField(network,temperatureInlet,temperatureOutlet,temperatureInletLinks,temperatureOutletLinks) 
    %Temperature field resulting from a temperature difference between 
    %temperature inlet and temperature outlet

    heatDiffusivity = 1; % TODO : check value
    conductancesHeat = LocalScaleComputeConductancesDiffusion(network,heatDiffusivity);
    
    boundaryConditions=struct;
    boundaryConditions.inletLink = temperatureInletLinks;
    boundaryConditions.outletLink = temperatureOutletLinks;
    boundaryConditions.inletType = 'Dirichlet' ;
    boundaryConditions.outletType = 'Dirichlet' ;
    boundaryConditions.inletValue = temperatureInlet*ones(1,length(boundaryConditions.inletLink));
    boundaryConditions.outletValue = temperatureOutlet*ones(1,length(boundaryConditions.outletLink));

    transportPores = 1:network.GetNumberOfPores ;

    [ temperature, ~,heatTransferCoefficient ]=ComputeLinearTransport(network,transportPores,conductancesHeat,boundaryConditions);


end
