function network=CreateGDLNetwork(dimension,nPore,xLength,yLength,zLength)


    myGeometry=MacroscopicGeometry();

    
    myGeometry.SetNetworkType('PoreNetworkMeshFibrous') ;
    myGeometry.SetPavageType('RandomVoronoi');


    if dimension==3
        
        myGeometry.SetDimension(3);
        [vertices,~, faces] = createCube ;
        vertices(:,1) = xLength*vertices(:,1);
        vertices(:,2) = yLength*vertices(:,2);
        vertices(:,3) = zLength*vertices(:,3);
        
    elseif dimension==2
        
        myGeometry.SetDimension(2);
        vertices = [0 0 ; 1 0 ; 1 1 ; 0 1 ];
        faces=[1 2 ; 2 3 ; 3 4 ; 4 1];
        vertices(:,1) = xLength*vertices(:,1);
        vertices(:,2) = yLength*vertices(:,2);

    end
    
    myGeometry.AddVertices(vertices);

    nBoundary=size(faces,1);
    for iBoundary=1:nBoundary
        boundary.Face=faces(iBoundary,:);
        boundary.Rugosity='rough';
        boundary.Type='surface';
        boundary.Name='';
        myGeometry.AddBoundary(boundary);
    end

    block.VertexNumbers=1:size(vertices,1);
    block.Filling.PoreNumber=nPore;
    block.Filling.FiberThickness=0.01;
    myGeometry.AddBlock(block)
    
    


    networkBuilder=NetworkBuilder(myGeometry);

    network=networkBuilder.BuildNetwork();


end