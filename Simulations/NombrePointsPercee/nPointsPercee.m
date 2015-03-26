function [nBreakthroughList,nInletLinkList] = nPointsPercee(thicknessList,randomSeedList)


    %Calcule le nombre de points de percée lors d'une invasion avec 100% de
    %points d'injection.
    

    
    nIteration = length(thicknessList)*length(randomSeedList);
    
    nBreakthroughList = zeros(1,nIteration);
    nInletLinkList = zeros(1,nIteration); 
    
    
    parfor iteration=1:nIteration
             
        [a,b]=ind2sub([length(thicknessList),length(randomSeedList)],iteration);
        thickness = thicknessList(a);
        randomSeed = randomSeedList(b);

        Lx = 1;
        Lz = thickness/10;
        nPore = 4000*Lz;

        network = BuildNetwork(nPore,Lx,Lz,randomSeed);


        [clusterOptions,nCluster,clustersInletLink,clustersOutletLink] = SetIPoptions(network);

        clusters = ComputeInvasionPercolationSeveralClusters( ...
                                    network,nCluster,...
                                    clustersInletLink,clustersOutletLink,...
                                    'currentWettability',clusterOptions );



        %Compute number of breakthrough points

        nBreakthroughList(iteration) = length(clusters);


        nInletLinkList(iteration) = nCluster;

    end
    
    %plot results
    figure;
    scatter(nInletLinkList,nBreakthroughList./(nInletLinkList).^2);
    
end



function network = BuildNetwork(nPore,Lx,Lz,randomSeed)

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





