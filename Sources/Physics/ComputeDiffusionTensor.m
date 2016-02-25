function diffusionTensor = ComputeDiffusionTensor( network, phaseConductances)
%COMPUTEDIFFUSIONTENSOR Assumes a rectangular domain. 
%   Input: - network
%          - phaseConductances(i) : conductance of phase 
%   
%   Output : diffusion tensor normalized by the network dimensions
%   
    
    
    %Computes conductances for diffusion
    
    parameters.GeometricModel.Pore = 'Cylinder';
    if isfield(network.GetPoreDataList,'Phase')
        phaseCode=network.GetPoreData('Phase');
    else
        phaseCode=0;
    end    
    assert(length(phaseConductances)>max(phaseCode),'conductance must be specified for all phases')
    parameters.PoreBulkProp = phaseConductances(phaseCode+1);
    parameters.LinkBulkProp = 0;
    parameters.GeometricModel.Link = 'None';
    
    conductancesDiffusion = LocalScaleComputeConductancesDiffusion(network,parameters);
    
    
    %Perform diffusion simulations
    
    transportPores = 1:network.GetNumberOfPores ;
    boundaryConditions.inletType = 'Dirichlet' ;
    boundaryConditions.outletType = 'Dirichlet' ;
    
    [boundaryXYZ,lengthXYZ,surfaceXYZ] = network.GetCuboidDomain;
    
    % X direction
    boundaryConditions.inletLink = network.GetLinksFrontiere(boundaryXYZ{1});
    boundaryConditions.outletLink = network.GetLinksFrontiere(boundaryXYZ{2});
    
    nInletLink = length(boundaryConditions.inletLink);
    nOutletLink = length(boundaryConditions.outletLink);
    boundaryConditions.inletValue = 1*ones(1,nInletLink);
    boundaryConditions.outletValue = 0.1*ones(1,nOutletLink);
    
    [ ~, ~, effConductanceX ] = ComputeLinearTransport( ...
                network,transportPores, ...
                conductancesDiffusion, ...
                boundaryConditions  ...
                );
    
    normalizedConductanceX = effConductanceX*(lengthXYZ(1)/surfaceXYZ(1));         
    
    
    
    % Y direction
    boundaryConditions.inletLink = network.GetLinksFrontiere(boundaryXYZ{3});
    boundaryConditions.outletLink = network.GetLinksFrontiere(boundaryXYZ{4});
    
    nInletLink = length(boundaryConditions.inletLink);
    nOutletLink = length(boundaryConditions.outletLink);
    boundaryConditions.inletValue = 1*ones(1,nInletLink);
    boundaryConditions.outletValue = 0.1*ones(1,nOutletLink);
    
    [ ~, ~, effConductanceY ] = ComputeLinearTransport( ...
                network,transportPores, ...
                conductancesDiffusion, ...
                boundaryConditions  ...
                );
    
    normalizedConductanceY = effConductanceY*(lengthXYZ(2)/surfaceXYZ(2));        
    
    % Z direction
    if network.GetDimension>2
        
        boundaryConditions.inletLink = network.GetLinksFrontiere(boundaryXYZ{5});
        boundaryConditions.outletLink = network.GetLinksFrontiere(boundaryXYZ{6});
        
        nInletLink = length(boundaryConditions.inletLink);
        nOutletLink = length(boundaryConditions.outletLink);
        boundaryConditions.inletValue = 1*ones(1,nInletLink);
        boundaryConditions.outletValue = 0.1*ones(1,nOutletLink);
        
        [ ~, ~, effConductanceZ ] = ComputeLinearTransport( ...
                    network,transportPores, ...
                    conductancesDiffusion, ...
                    boundaryConditions  ...
                    );
        
        normalizedConductanceZ = effConductanceZ*(lengthXYZ(3)/surfaceXYZ(3));
        
    end
    
    if network.GetDimension==2
        diffusionTensor=[normalizedConductanceX,normalizedConductanceY];
    elseif network.GetDimension==3
        diffusionTensor=[normalizedConductanceX,normalizedConductanceY,normalizedConductanceZ];
    end
	
end



