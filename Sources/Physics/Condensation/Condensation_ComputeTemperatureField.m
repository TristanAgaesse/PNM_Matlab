function  [temperature,heatTransferCoefficient] = Condensation_ComputeTemperatureField(network,temperatureInlet,temperatureOutlet,temperatureInletLinks,temperatureOutletLinks,temperatureTransportPores,poreHeatConductivity) 
    %Temperature field resulting from a temperature difference between 
    %temperature inlet and temperature outlet

    
    
    
    parameters.GeometricModel.Pore = 'Cylinder' ;
    parameters.GeometricModel.Link = 'None' ;% 'SurfaceResistance_RealSurface'
    nPore = network.GetNumberOfPores;
    assert(length(poreHeatConductivity)==nPore);
    parameters.PoreBulkProp = poreHeatConductivity; %  heatDiffusivity = 1;  poreHeatDiffusivityheatDiffusivity*ones(nPore,1);
    parameters.LinkBulkProp =0; % scalar or array(nLink,1)
    
    
    
    conductancesHeat = LocalScaleComputeConductancesDiffusion(network,parameters);
    
    boundaryConditions=struct;
    boundaryConditions.inletLink = temperatureInletLinks;
    boundaryConditions.outletLink = temperatureOutletLinks;
    boundaryConditions.inletType = 'Dirichlet' ;
    boundaryConditions.outletType = 'Dirichlet' ;
    boundaryConditions.inletValue = temperatureInlet*ones(1,length(boundaryConditions.inletLink));
    boundaryConditions.outletValue = temperatureOutlet*ones(1,length(boundaryConditions.outletLink));
    boundaryConditions.solver='mldivide';

    [ temperature, ~,heatTransferCoefficient ]=ComputeLinearTransport(network,temperatureTransportPores,conductancesHeat,boundaryConditions);


end
