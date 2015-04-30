
%% Initialisation

% [network,viewer] = CreateNetwork('GDL_2D');
% 
% inletLink = network.GetLinksFrontiere([1 2 3]) ; % GDL/MPL interface
% outletLink = network.GetLinksFrontiere([4 6]);   % Channel

dimension=3;
nPore=8000;
xLength=1;
yLength=1;
zLength=0.2;
network=CreateGDLNetwork(dimension,floor(nPore*zLength),xLength,yLength,zLength);

inletLink = network.GetLinksFrontiere(1) ; % GDL/MPL interface
outletLink = network.GetLinksFrontiere(2);   % Channel

clusterOptions.ThroatPressure = 'LaplaceCylinder' ;
clusterOptions.SurfaceTension = 60e-3; %Water/air surface tension at 80°C


%Degradation parameters
options.DegradationMechanism = 'uniform';

options.nIterations = 1000 ;
options.DeltaContactAngle = 30/180*pi;
options.ClusterGrowth = true  ;
options.InvasionTimeSearch = 'exactTimeForLaplaceLaw' ;


%Varying parameters
%dropletPressure=[0.2,0.5,0.9,1.3];
dropletPressure = 0.2;
%contactAngles
coalescence = {'numberOfInvadedNeighbours','none'};


nLink = network.GetNumberOfLinks;

contactAngle = {pi*(100+20*rand(1,nLink))/180.*ones(1,nLink),pi*(80+20*rand(1,nLink))/180.*ones(1,nLink)}; %MAKE IT VARY

simuTitle = {'1','2','3','4'};

for i =1
    %% Computing degradation
    
    %Varying initial hydrophobicity parameter
    
    network.AddNewLinkData(contactAngle{1},'ContactAngle')
    
    clusterOptions.Coalescence = coalescence{1};  %MAKE IT VARY
    
        
    capillaryPressureInitiale = ComputeCapillaryPressureCurve(network,...
                                         inletLink,'currentWettability',clusterOptions);
    
    
    [initialCluster,breakthroughPressure,invasionPressureList] = ComputeInvasionPercolation(...
                                         network,inletLink,outletLink,...
                                        'currentWettability',clusterOptions);
    
    %Varying cluster boundary condition
    clusterPressure=dropletPressure*invasionPressureList(end);   %MAKE IT VARY
                                    
                                    
    outputInformation=ComputeHydrophobicityLoss(network,initialCluster,clusterPressure,...
                                                inletLink,outletLink,options);
    
    
    %% Post processing                                        
    infoDegradation=postTraitementDegradation(network,outputInformation);


        %Image degradation
%     figure
%     viewer.View('PoreField',infoDegradation.PourcentageDegradationMoyennePore{end})

       
        %Capillary Pressure après degradation

    capillaryPressureEnd = ComputeCapillaryPressureCurve(network,inletLink,...
                                            'currentWettability',clusterOptions);

    figure
    plot(capillaryPressureInitiale(:,1),capillaryPressureInitiale(:,2),'black');hold on 
    plot(capillaryPressureEnd(:,1),capillaryPressureEnd(:,2),'red');hold on 
    title(strcat('Capillary Pressure ',simuTitle{i}))
    
    figure
    plot(infoDegradation.Times)
    title(strcat('Times ',simuTitle{i}))
    
    figure
    plot(diff(infoDegradation.Times))
    title(strcat('DiffInvasionTime ',simuTitle{i}))
    
    figure
    plot(infoDegradation.DiffusionCoefficient)
    title(strcat('DiffusionCoefficient ',simuTitle{i}))
    
    figure
    plot(infoDegradation.Pressure)
    title(strcat('Pressure ',simuTitle{i}))
    
    figure
    plot(infoDegradation.nInvasionPerStep)
    title(strcat('nInvasionPerBurst ',simuTitle{i}))
    
    figure
    plot(infoDegradation.Saturation)
    title(strcat('Saturation ',simuTitle{i}))
    
%     figure
%     plot(infoDegradation.BurstVolume)
%     title(strcat('BurstVolume ',simuTitle{i}))
%     
%     figure
%     plot(infoDegradation.InvadedPoreVolume)
%     title(strcat('InvadedPoreVolume ',simuTitle{i}))
        
%     figure
%     data=infoDegradation.InvasionTimeDistribution{1};
%     [a,b]=hist(data,30);
%     plot(b,a/length(data),'black');hold on 
%     data=infoDegradation.InvasionTimeDistribution{outputInformation.nIteration};
%     [a,b]=hist(data,30);
%     plot(b,a/length(data),'red');hold on
%     data=infoDegradation.InvasionTimeDistribution{floor(outputInformation.nIteration/3)};
%     [a,b]=hist(data,30);
%     plot(b,a/length(data),'green');hold on
%     data=infoDegradation.InvasionTimeDistribution{floor(outputInformation.nIteration*2/3)};
%     [a,b]=hist(data,30);
%     plot(b,a/length(data),'yellow');hold on
%     title(strcat('InvasionTimeDistribution ',simuTitle{i}))
% 
%     figure
%     plot(infoDegradation.MeanInvasionTime)
%     title(strcat('MeanInvasionTime ',simuTitle{i}))
    
    figure
    data=infoDegradation.ContactAngle{1};
    [a,b]=hist(data*180/pi,6);
    plot(b,a/length(data),'black');hold on 
    data=infoDegradation.ContactAngle{outputInformation.nIteration};
    [a,b]=hist(data*180/pi,6);
    plot(b,a/length(data),'red');hold on
    data=infoDegradation.ContactAngle{floor(outputInformation.nIteration/3)};
    [a,b]=hist(data*180/pi,6);
    plot(b,a/length(data),'green');hold on
    data=infoDegradation.ContactAngle{floor(outputInformation.nIteration*2/3)};
    [a,b]=hist(data*180/pi,6);
    plot(b,a/length(data),'yellow');hold on
    title(strcat('ContactAngle ',simuTitle{i}))
    
    figure
    plot(infoDegradation.MeanContactAngle)
    title(strcat('MeanContactAngle ',simuTitle{i}))
    
    
    figure
    data=infoDegradation.InterfaceLinkContactAngle{1};
    [a,b]=hist(data*180/pi,6);
    plot(b,a/length(data),'black');hold on 
    data=infoDegradation.InterfaceLinkContactAngle{outputInformation.nIteration};
    [a,b]=hist(data*180/pi,6);
    plot(b,a/length(data),'red');hold on
    data=infoDegradation.InterfaceLinkContactAngle{floor(outputInformation.nIteration/3)};
    [a,b]=hist(data*180/pi,6);
    plot(b,a/length(data),'green');hold on
    data=infoDegradation.InterfaceLinkContactAngle{floor(outputInformation.nIteration*2/3)};
    [a,b]=hist(data*180/pi,6);
    plot(b,a/length(data),'yellow');hold on
    title(strcat('InterfaceLinkContactAngle ',simuTitle{i}))
    
    figure
    plot(infoDegradation.InterfaceLinkMeanContactAngle)
    title(strcat('InterfaceLinkMeanContactAngle ',simuTitle{i}))
    
    
    figure
    data=infoDegradation.InterfaceLinkDiameter{1};
    [a,b]=hist(data,6);
    plot(b,a/length(data),'black');hold on 
    data=infoDegradation.InterfaceLinkDiameter{outputInformation.nIteration};
    [a,b]=hist(data,6);
    plot(b,a/length(data),'red');hold on
    data=infoDegradation.InterfaceLinkDiameter{floor(outputInformation.nIteration/3)};
    [a,b]=hist(data,6);
    plot(b,a/length(data),'green');hold on
    data=infoDegradation.InterfaceLinkDiameter{floor(outputInformation.nIteration*2/3)};
    [a,b]=hist(data,6);
    plot(b,a/length(data),'yellow');hold on
    title(strcat('InterfaceLinkDiameter ',simuTitle{i}))
    
    figure
    plot(infoDegradation.InterfaceLinkMeanDiameter)
    title(strcat('InterfaceLinkMeanDiameter ',simuTitle{i}))
    
    figure
    [a,b]=hist(network.GetLinkData('Diameter'),6);
    plot(b,a/network.GetNumberOfLinks);
    title(strcat('LinkDiameter ',simuTitle{i}))
end




