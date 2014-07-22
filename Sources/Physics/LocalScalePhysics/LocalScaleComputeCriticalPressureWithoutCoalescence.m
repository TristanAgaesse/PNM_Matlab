function Pc = LocalScaleComputeCriticalPressureWithoutCoalescence(network,iLien)
    %Pressure on toroid

    sigmaWater = 60e-3;
    theta = network.GetLinkDataList.ContactAngle(iLien);
    linkDiameter = network.GetLinkDataList.Diameter(iLien);
    
    if isa(network,'PoreNetworkMeshFibrous')
        numEdges = network.FacesToEdges{iLien};
        fibreDiameter = mean(network.GetEdgeDataList.FiberDiameter(numEdges));
    else
        fibreDiameter=linkDiameter/100;
    end

    alpha = theta-pi+asin(sin(theta)/(1+linkDiameter/fibreDiameter));
    Pc = -(4*sigmaWater/linkDiameter)*cos(theta-alpha)/(1+fibreDiameter*(1-cos(alpha))/(linkDiameter));
end