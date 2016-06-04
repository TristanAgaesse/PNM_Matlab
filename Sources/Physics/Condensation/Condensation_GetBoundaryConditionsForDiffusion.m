function [gasTransportPores,inletVaporPressure,outletVaporPressure,vaporInletLinks,...
                vaporOutletLinks] = GetBoundaryConditionsForDiffusion(network,options,invadedPore)

    gasTransportPores = 1:network.GetNumberOfPores;

    inletVaporPressure = options.RelativeHumidityInlet*options.AirPressure; %TODO : Check this formula
    outletVaporPressure = options.RelativeHumidityOutlet*options.AirPressure;
        
    vaporInletLinks=options.VaporInletLinks;
    vaporOutletLinks=options.VaporOutletLinks;    
    
    %TODO : change this when some pores are invaded
    
    % Attention, fonction utilisée dans 2 contextes : nuléation et
    % diffusionCondensation 
            
end

