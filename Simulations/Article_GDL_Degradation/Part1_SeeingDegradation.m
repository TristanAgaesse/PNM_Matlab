
[network,viewer]=CreateNetwork('GDL_2D');

%Assign contact angle
contactAngle=pi*(100+20*rand(1,network.GetNumberOfLinks))/180.*ones(1,network.GetNumberOfLinks);
network.AddNewLinkData(contactAngle,'ContactAngle')

%Capillary pressure initiale
clusterOptions.ThroatPressure = 'LaplaceCylinder' ;
clusterOptions.Coalescence = 'numberOfInvadedNeighbours';
clusterOptions.SurfaceTension = 60e-3; %Water/air surface tension at 80°C

capillaryPressureInitiale1 = ComputeCapillaryPressureCurve(network,network.GetLinksFrontiere([1 2 3]),'currentWettability',clusterOptions);
capillaryPressureInitiale2 = ComputeCapillaryPressureCurve(network,network.GetLinksFrontiere([4 5 6]),'currentWettability',clusterOptions);


%Invasion initiale du reseau
inletLink = network.GetLinksFrontiere([1 2 3]) ; % GDL/MPL interface
outletLink = network.GetLinksFrontiere([4 6]);   % Channel

[initialCluster,breakthroughPressure,invasionPressureList] = ComputeInvasionPercolation(network,inletLink,outletLink,'currentWettability',clusterOptions);



%Degradation
degradationMechanism={'uniform','uniformInWater','waterSpeed'};
for i =1:3
    options.nIterations = 300 ;
    options.DeltaContactAngle=30/180*pi;
    options.ClusterGrowth = false  ;
    options.DegradationMechanism=degradationMechanism{i};
    options.InvasionTimeSearch='exactTimeForLaplaceLaw' ;

    clusterPressure=0.9*invasionPressureList(end);

    outputInformation=ComputeHydrophobicityLoss(network,initialCluster,clusterPressure,inletLink,outletLink,options);

    infoDegradation=postTraitementDegradation(network,outputInformation);
    %save('DegradationUniformeEau','infosUniformeEau')


    %Post processing


    %Image degradation
    figure
    viewer.View('PoreField',infoDegradation.PourcentageDegradationMoyennePore{end})

    %Capillary Pressure après degradation

    capillaryPressureEnd1 = ComputeCapillaryPressureCurve(network,network.GetLinksFrontiere([1 2 3]),'currentWettability',clusterOptions);
    capillaryPressureEnd2 = ComputeCapillaryPressureCurve(network,network.GetLinksFrontiere([4 5 6]),'currentWettability',clusterOptions);

    figure
    plot(capillaryPressureInitiale1(:,1),capillaryPressureInitiale1(:,2),'black');hold on 
    plot(capillaryPressureInitiale2(:,1),capillaryPressureInitiale2(:,2),'red');hold on 
    plot(capillaryPressureEnd1(:,1),capillaryPressureEnd1(:,2),'green');hold on 
    plot(capillaryPressureEnd2(:,1),capillaryPressureEnd2(:,2),'blue');hold on 

end

