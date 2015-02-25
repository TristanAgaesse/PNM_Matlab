
[network,viewer]=CreateNetwork('GDL_2D');

inletLink = network.GetLinksFrontiere([1 2 3]) ; % GDL/MPL interface
outletLink = network.GetLinksFrontiere([4 6]);   % Channel


dropPressure=[0.5,0.9,1.4,2];


for i =1:3
    %Invasion initiale du reseau
    contactAngle=pi*(100+20*rand(1,network.GetNumberOfLinks))/180.*ones(1,network.GetNumberOfLinks); %MAKE IT VARY
    network.AddNewLinkData(contactAngle,'ContactAngle')
    
    clusterOptions.ThroatPressure = 'LaplaceCylinder' ;
    clusterOptions.SurfaceTension = 60e-3; %Water/air surface tension at 80°C

    clusterOptions.Coalescence = 'numberOfInvadedNeighbours';  %MAKE IT VARY
    
    capillaryPressureInitiale = ComputeCapillaryPressureCurve(network,network.GetLinksFrontiere([1 2 3]),'currentWettability',clusterOptions);
    
    [initialCluster,breakthroughPressure,invasionPressureList] = ComputeInvasionPercolation(network,inletLink,outletLink,'currentWettability',clusterOptions);

    %Degradation
    options.nIterations = 300 ;
    options.DeltaContactAngle=30/180*pi;
    options.ClusterGrowth = true  ;
    options.DegradationMechanism='uniform';
    options.InvasionTimeSearch='exactTimeForLaplaceLaw' ;

    clusterPressure=dropPressure(i)*invasionPressureList(end);   %MAKE IT VARY

    outputInformation=ComputeHydrophobicityLoss(network,initialCluster,clusterPressure,inletLink,outletLink,options);

    infoDegradation=postTraitementDegradation(network,outputInformation);

    %Post processing


    %Image degradation
    figure
    viewer.View('PoreField',infoDegradation.PourcentageDegradationMoyennePore{end})

    %Capillary Pressure après degradation

    capillaryPressureEnd = ComputeCapillaryPressureCurve(network,network.GetLinksFrontiere([1 2 3]),'currentWettability',clusterOptions);

    figure
    plot(capillaryPressureInitiale(:,1),capillaryPressureInitiale(:,2),'black');hold on 
    plot(capillaryPressureEnd(:,1),capillaryPressureEnd(:,2),'red');hold on 
end