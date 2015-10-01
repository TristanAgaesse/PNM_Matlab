function Pc = LocalScaleComputeCriticalPressureWithoutCoalescence(network,liens,clusterOptions)
    

    if isempty(liens)
        Pc=[];
        return
    end
    
    if isfield(clusterOptions,'SurfaceTension')
        sigma=clusterOptions.('SurfaceTension');
    else
        sigma = 60e-3;   %Surface tension water/air at 80°C
    end
    
    theta = network.GetLinkDataList.ContactAngle(liens);
    linkDiameter = network.GetLinkDataList.Diameter(liens);
    
    switch clusterOptions.CapillaryPressureLaw
        
        case 'LaplaceCylinder'
            %Laplace law in cylinders
            Pc = -(4*sigma*cos(theta)./linkDiameter);

        case 'PurcellToroid'

        %Pressure on toroid

            fibreDiameter=zeros(length(liens),1);
            for i=1:length(liens)
                if isa(network,'PoreNetworkMeshFibrous')
                    numEdges = network.FacesToEdges{liens(i)};
                    fibreDiameter(i) = sum(network.GetEdgeDataList.FiberDiameter(numEdges))/length(numEdges);
                elseif isfield(clusterOptions,'FiberDiameterForPurcell')
                    fibreDiameter(i)=clusterOptions.FiberDiameterForPurcell;
                else
                    fibreDiameter(i)=linkDiameter(i)/20;
                end
            end

            alpha = theta-pi+asin(sin(theta)./(1+linkDiameter./fibreDiameter));
            Pc = -(4*sigma./linkDiameter).*cos(theta-alpha)./(1+fibreDiameter.*(1-cos(alpha))./(linkDiameter));
    end
    
    Pc(linkDiameter==0)=Inf; %instead of -Inf computed for hydrophilic links with linkDiameter=0
    
    assert(isempty(find(isnan(Pc),1)),'Nan found for a capillary pressure');

    Pc=transpose(Pc);
end