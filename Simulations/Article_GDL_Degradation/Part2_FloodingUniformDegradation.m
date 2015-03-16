
%% Initialisation

[network,viewer]=CreateNetwork('GDL_2D');

inletLink = network.GetLinksFrontiere([1 2 3]) ; % GDL/MPL interface
outletLink = network.GetLinksFrontiere([4 6]);   % Channel

clusterOptions.ThroatPressure = 'LaplaceCylinder' ;
clusterOptions.SurfaceTension = 60e-3; %Water/air surface tension at 80°C


%Degradation parameters
options.DegradationMechanism='uniform';

options.nIterations = 300 ;
options.DeltaContactAngle=30/180*pi;
options.ClusterGrowth = true  ;
options.InvasionTimeSearch='exactTimeForLaplaceLaw' ;


%Varying parameters
dropletPressure=[0.2,0.5,0.9,1.3];
%contactAngles
%coalescence


nLink = network.GetNumberOfLinks;

contactAngle=pi*(100+20*rand(1,nLink))/180.*ones(1,nLink); %MAKE IT VARY

for i =1:3
    %% Computing degradation
    
    %Varying initial hydrophobicity parameter
    
    network.AddNewLinkData(contactAngle,'ContactAngle')
    
    clusterOptions.Coalescence = 'numberOfInvadedNeighbours';  %MAKE IT VARY
    
    
    
        
    capillaryPressureInitiale = ComputeCapillaryPressureCurve(network,...
                                         inletLink,'currentWettability',clusterOptions);
    
    
    [initialCluster,breakthroughPressure,invasionPressureList] = ComputeInvasionPercolation(...
                                         network,inletLink,outletLink,...
                                        'currentWettability',clusterOptions);
    
    %Varying cluster boundary condition
    clusterPressure=dropletPressure(i)*invasionPressureList(end);   %MAKE IT VARY
                                    
                                    
    outputInformation=ComputeHydrophobicityLoss(network,initialCluster,clusterPressure,...
                                                inletLink,outletLink,options);
    
    
    %% Post processing                                        
    infoDegradation=postTraitementDegradation(network,outputInformation);


    

        %Image degradation
    figure
    viewer.View('PoreField',infoDegradation.PourcentageDegradationMoyennePore{end})

        %Capillary Pressure après degradation

    capillaryPressureEnd = ComputeCapillaryPressureCurve(network,inletLink,...
                                            'currentWettability',clusterOptions);

    figure
    plot(capillaryPressureInitiale(:,1),capillaryPressureInitiale(:,2),'black');hold on 
    plot(capillaryPressureEnd(:,1),capillaryPressureEnd(:,2),'red');hold on 
end




