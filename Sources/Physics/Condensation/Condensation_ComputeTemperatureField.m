function  [temperature,heatTransferCoefficient] = Condensation_ComputeTemperatureField(network,heatOptions) 
    %Temperature field resulting from a temperature difference between 
    %temperature inlet and temperature outlet
    
    
    
    
    parameters.GeometricModel.Pore = 'Cylinder' ;
    parameters.GeometricModel.Link = 'None' ;% 'SurfaceResistance_RealSurface'
    nPore = network.GetNumberOfPores;
    assert(length(heatOptions.TemperaturePoreHeatConductivity)==nPore);
    parameters.PoreBulkProp = heatOptions.TemperaturePoreHeatConductivity; %  heatDiffusivity = 1;  poreHeatDiffusivityheatDiffusivity*ones(nPore,1);
    parameters.LinkBulkProp =0; % scalar or array(nLink,1)
    
    
    
    conductancesHeat = LocalScaleComputeConductancesDiffusion(network,parameters);
    
    boundaryConditions=struct;
    boundaryConditions.inletLink = heatOptions.TemperatureInletLinks;
    boundaryConditions.outletLink = heatOptions.TemperatureOutletLinks;
    boundaryConditions.inletType = heatOptions.TemperatureInletType;   %'Dirichlet' ;
    boundaryConditions.outletType = heatOptions.TemperatureOutletType;    %'Dirichlet' ;
    boundaryConditions.inletValue = heatOptions.TemperatureInlet;     %*ones(1,length(boundaryConditions.inletLink));
    boundaryConditions.outletValue = heatOptions.TemperatureOutlet ;  %*ones(1,length(boundaryConditions.outletLink));
    boundaryConditions.solver='mldivide';

    [ temperature, ~,heatTransferCoefficient ]=ComputeLinearTransport(network,heatOptions.TemperatureTransportPores,conductancesHeat,boundaryConditions);


end
