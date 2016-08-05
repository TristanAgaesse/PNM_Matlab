function diffusionTensor = ComputeDiffusionTensor( network, phaseConductances)
%COMPUTEDIFFUSIONTENSOR Assumes a rectangular domain. 
%   Input: - network
%          - phaseConductances(i) : conductance of phase 
%   
%   Output : diffusion tensor normalized by the network dimensions
%   
% Example :
%phaseConductances=zeros(1,256);
%phaseConductances(1)=1;
%phaseConductances(256)=10;
%diffusionTensor = ComputeDiffusionTensor( network, phaseConductances)
    
    %Computes conductances for diffusion
    
    parameters.GeometricModel.Pore = 'Cylinder';%'Cylinder' '2Cubes' 'PolynomialProfileDegree1' 'VolumeConservation_LinearProfile'  'VolumeConservation_ConstantProfile'
    if isfield(network.GetPoreDataList,'Phase')
        phaseCode=network.GetPoreData('Phase');
    else
        phaseCode=0;
    end    
    assert(length(phaseConductances)>max(phaseCode),'conductance must be specified for all phases')
    parameters.PoreBulkProp = phaseConductances(phaseCode+1);
    parameters.LinkBulkProp = 0;
    parameters.GeometricModel.Link = 'None';
    
    
    linkSurface_voxelUnit = double(network.GetLinkData('RawData_GeometricSurface'));
    linkSurface = linkSurface_voxelUnit*(network.GetVoxelEdgeLength)^2;
    network.AddNewLinkData(linkSurface,'Surface');
    
    
    conductancesDiffusion = LocalScaleComputeConductancesDiffusion(network,parameters);
    
    boundaryConditions.solver='mldivide';
    
    %Perform diffusion simulations
    
    transportPores = 1:network.GetNumberOfPores ;
    boundaryConditions.inletType = 'Dirichlet' ;
    boundaryConditions.outletType = 'Dirichlet' ;
    %boundaryConditions.outletType = 'Neumann' ;
    
    [boundaryXYZ,lengthXYZ,surfaceXYZ] = network.GetCuboidDomain;
    
    % X direction
    boundaryConditions.inletLink = network.GetLinksFrontiere(boundaryXYZ{1});
    boundaryConditions.outletLink = network.GetLinksFrontiere(boundaryXYZ{2});
    
    nInletLink = length(boundaryConditions.inletLink);
    nOutletLink = length(boundaryConditions.outletLink);
    boundaryConditions.inletValue = 1*ones(1,nInletLink);
    boundaryConditions.outletValue = 0*ones(1,nOutletLink);
%     allLinkSurfaces = transpose(network.GetLinkData('Surface'));
%     surfacicFlux = 1;
%     boundaryConditions.outletValue = surfacicFlux*allLinkSurfaces(boundaryConditions.outletLink);
    
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
    boundaryConditions.outletValue = 0*ones(1,nOutletLink);
%     allLinkSurfaces = transpose(network.GetLinkData('Surface'));
%     surfacicFlux = 1;
%     boundaryConditions.outletValue = surfacicFlux*allLinkSurfaces(boundaryConditions.outletLink);
    
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
        boundaryConditions.outletValue = 0*ones(1,nOutletLink);
%         allLinkSurfaces = transpose(network.GetLinkData('Surface'));
%         surfacicFlux = 1;
%         boundaryConditions.outletValue = surfacicFlux*allLinkSurfaces(boundaryConditions.outletLink);
        
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
