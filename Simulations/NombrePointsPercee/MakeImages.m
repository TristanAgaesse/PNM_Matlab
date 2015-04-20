
function MakeImages()


    %Cas 2D

    Lx = 1;
    Ly = 0.2;
    nPore = floor(3000*Ly);
    randomSeed=0;
    
    network = BuildNetwork2D(nPore,Lx,Ly,randomSeed);

    viewer = Viewer(network.PrivateInternalOutputStruct);
    
    [clusterOptions,nCluster,clustersInletLink,clustersOutletLink] = SetIPoptions(network);

    clusters = ComputeInvasionPercolationSeveralClusters( ...
                                network,nCluster,...
                                clustersInletLink,clustersOutletLink,...
                                'currentWettability',clusterOptions );

    invadedPores=zeros(1,network.GetNumberOfPores);
    for i =1:length(clusters)
        invadedPores(clusters{i}.GetInvadedPores)=i;
    end
    network.AddNewPoreData(invadedPores,'ClustersInvadedPores')
    
    viewer.View('PoreField',invadedPores);
    
    %Cas 3D

%     Lx = 1;
%     Lz = 0.2;
%     nPore = floor(30000*Lz);
% 
%     randomSeed=0;
%     
%     network = BuildNetwork3D(nPore,Lx,Lz,randomSeed);
% 
% 
%     [clusterOptions,nCluster,clustersInletLink,clustersOutletLink] = SetIPoptions(network);
% 
%     clusters = ComputeInvasionPercolationSeveralClusters( ...
%                                 network,nCluster,...
%                                 clustersInletLink,clustersOutletLink,...
%                                 'currentWettability',clusterOptions );
% 
%     invadedPores=zeros(1,network.GetNumberOfPores);
%     for i =1:length(clusters)
%         invadedPores(clusters{i}.GetInvadedPores)=i;
%     end
% 
%     network.AddNewPoreData(invadedPores,'ClustersInvadedPores')
% 
%     network.ExportToParaview('3DclustersVisualization.vtk')

end



function network = BuildNetwork3D(nPore,Lx,Lz,randomSeed)

    myGeometry=MacroscopicGeometry();

    myGeometry.SetDimension(3);
    myGeometry.SetConvertToMeters(1);  %scale unit =  1m
    myGeometry.SetNetworkType('PoreNetworkMeshFibrous') ;
    myGeometry.SetPavageType('RandomVoronoi');
    myGeometry.SetRandomSeed(randomSeed);


    [V, ~, F] = createCube;

    V(:,1) = Lx*V(:,1);
    V(:,2) = Lx*V(:,2);
    V(:,3) = Lz*V(:,3);

    myGeometry.AddVertices(V);

    block=struct; 
    block.VertexNumbers=1:size(V,1);
    block.Filling.PoreNumber=nPore;
    block.Filling.FiberThickness=0.0001;
    myGeometry.AddBlock(block)

    for iBoundary=1:size(F,1)
        boundary=struct;
        boundary.Face=F(iBoundary,:);

        if iBoundary==1
            boundary.Rugosity='rough';
        else
            boundary.Rugosity='rough';
        end

        boundary.Type='surface';
        boundary.Name='';
        myGeometry.AddBoundary(boundary);
    end



    networkBuilder=NetworkBuilder(myGeometry);

    network=networkBuilder.BuildNetwork();

end



function network = BuildNetwork2D(nPore,Lx,Ly,randomSeed)

    myGeometry=MacroscopicGeometry();

    myGeometry.SetDimension(2);
    myGeometry.SetConvertToMeters(1);  %scale unit =  1m
    myGeometry.SetNetworkType('PoreNetworkMeshFibrous') ;
    myGeometry.SetPavageType('RandomVoronoi');
    myGeometry.SetRandomSeed(randomSeed);


    V = [ 0, 0 ; Lx,0 ; Lx, Ly ; 0, Ly ];

    F = [ 1 2 ; 3 4 ; 2 3 ; 4 1 ];

    myGeometry.AddVertices(V);

    block=struct; 
    block.VertexNumbers=1:size(V,1);
    block.Filling.PoreNumber=nPore;
    block.Filling.FiberThickness=0.0001;
    myGeometry.AddBlock(block)

    for iBoundary=1:size(F,1)
        boundary=struct;
        boundary.Face=F(iBoundary,:);

        if iBoundary==1
            boundary.Rugosity='rough';
        else
            boundary.Rugosity='rough';
        end

        boundary.Type='surface';
        boundary.Name='';
        myGeometry.AddBoundary(boundary);
    end



    networkBuilder=NetworkBuilder(myGeometry);

    network=networkBuilder.BuildNetwork();

end



function [clusterOptions,nCluster,clustersInletLink,clustersOutletLink]=SetIPoptions(network)

    %Definissons les lois d'invasion et la mouillabilite
    clusterOptions.Coalescence = 'None';
    clusterOptions.CapillaryPressureLaw = 'LaplaceCylinder';
    clusterOptions.SurfaceTension = 72e-3 ; 

    contactAngle=110;
    theta=contactAngle*pi/180*ones(1,network.GetNumberOfLinks);
    network.AddNewLinkData(theta,'ContactAngle');

    %Define the injection points : here we are choosing that 100% of inlet links  
    %are independant injection points
    inletLink = network.GetLinksFrontiere(1);
    nCluster = length(inletLink);

    clustersInletLink = cell(1,nCluster);
    clustersOutletLink = cell(1,nCluster);
    for iCluster=1:nCluster
        clustersInletLink{iCluster} = inletLink(iCluster);
        clustersOutletLink{iCluster} = network.GetLinksFrontiere(2);
    end
end
