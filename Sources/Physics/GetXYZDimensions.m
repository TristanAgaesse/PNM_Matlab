function [boundaryXYZ,lengthXYZ,surfaceXYZ] = GetXYZDimensions(network)
    % boundaryXYZ{1 (resp.2,3,4,5,6) } = boundary code for Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
    % lengthXYZ = [Lx,Ly,Lz]
    % surfaceXYZ = 
    
    
    %Compute domain dimensions
    
    centers = network.GetLinkCenter(1:network.GetNumberOfLinks);
    dimension=network.GetDimension;
    
    Xmin = min(centers(:,1));
    Xmax = max(centers(:,1));
    Ymin = min(centers(:,2));
    Ymax = max(centers(:,2));
    lengthXYZ = [Xmax-Xmin,Ymax-Ymin];
    surfaceXYZ = [lengthXYZ(2),lengthXYZ(1)];
    if dimension==3
        Zmin = min(centers(:,3));
        Zmax = max(centers(:,3));
        lengthXYZ = [Xmax-Xmin,Ymax-Ymin,Zmax-Zmin];
        surfaceXYZ = [lengthXYZ(2)*lengthXYZ(3),lengthXYZ(1)*lengthXYZ(3),lengthXYZ(1)*lengthXYZ(2)];
    end
    
    
    % Assign a domain face to each network boundary
    
    nBoundary = network.GetNumberOfBoundaries;
    boundaryMean = zeros(nBoundary,network.GetDimension);
    for iBoundary = 1:nBoundary
        boundaryMean(iBoundary,:) = mean(network.GetLinkCenter(network.GetLinksFrontiere(iBoundary)));
    end
    
    boundaryMean(:,1) = boundaryMean(:,1)-Xmin;
    boundaryMean(:,1) = boundaryMean(:,1)./lengthXYZ(1);
    boundaryMean(:,2) = boundaryMean(:,2)-Ymin;
    boundaryMean(:,2) = boundaryMean(:,2)./lengthXYZ(2);
    if dimension==3
        boundaryMean(:,3) = boundaryMean(:,3)-Zmin;
        boundaryMean(:,3) = boundaryMean(:,3)./lengthXYZ(3);
    end
    
    boundaryXYZ = cell(1,2*dimension);   
    
    for iBoundary = 1:nBoundary
        
        if boundaryMean(iBoundary,1) < 0.1
            boundaryXYZ{1}=[boundaryXYZ{1},iBoundary];
        elseif boundaryMean(iBoundary,1) > 0.9
            boundaryXYZ{2}=[boundaryXYZ{2},iBoundary];
        end
        
        if boundaryMean(iBoundary,2) < 0.1
            boundaryXYZ{3}=[boundaryXYZ{3},iBoundary];
        elseif boundaryMean(iBoundary,2) > 0.9
            boundaryXYZ{4}=[boundaryXYZ{4},iBoundary];
        end
        
        if dimension==3
            if boundaryMean(iBoundary,3) < 0.1
                boundaryXYZ{5}=[boundaryXYZ{5},iBoundary];
            elseif boundaryMean(iBoundary,3) > 0.9
                boundaryXYZ{6}=[boundaryXYZ{6},iBoundary];
            end
        end
        
    end   
    
    checkList=zeros(1,nBoundary);
    for i = 1:length(boundaryXYZ)
        checkList(boundaryXYZ{i})=checkList(boundaryXYZ{i})+1;
    end
    assert(~any(checkList-ones(1,nBoundary)),'some boundaries are unassigned to a domain face') 
end

