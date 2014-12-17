function Pc = LocalScaleComputeCriticalPressureWithoutCoalescence(network,liens,clusterOptions)
    

    if isempty(liens)
        Pc=[];
        return
    end

    sigmaWater = 60e-3;
    theta = network.GetLinkDataList.ContactAngle(liens);
    linkDiameter = network.GetLinkDataList.Diameter(liens);
    
    switch clusterOptions.ThroatPressure
        
        case 'LaplaceCylinder'
            %Laplace law in cylinders
            Pc = -(4*sigmaWater*cos(theta)./linkDiameter);

        case 'PurcellToroid'

        %Pressure on toroid

            fibreDiameter=zeros(1,length(liens));
            for i=1:length(liens)
                if isa(network,'PoreNetworkMeshFibrous')
                    numEdges = network.FacesToEdges{liens(i)};
                    fibreDiameter(i) = mean(network.GetEdgeDataList.FiberDiameter(numEdges));
                else
                    fibreDiameter(i)=linkDiameter(i)/100;
                end
            end

            alpha = theta-pi+asin(sin(theta)./(1+linkDiameter./fibreDiameter));
            Pc = -(4*sigmaWater./linkDiameter).*cos(theta-alpha)./(1+fibreDiameter.*(1-cos(alpha))./(linkDiameter));
    end
end