function [condensationInfos,condensationClusters]=ComputeCondensation(network, options)
%COMPUTECONDENSATION Summary of this function goes here
%   Pore network condensation algorithm, as described in Benjamin
%   Straubhaar's PhD thesis.
%   
%Algorithme : 
%      - Demarrage de la condensation par une etape de nucleation :
%   les pore ayant une pression de vapeur superieure à la pression de 
%   vapeur saturante sont envahis par l'eau. Cela forme un ou plusieurs
%   clusters de condensation independants. Plus de details sur le calcul
%   de la pression de vapeur sont disponibles dans la fonction
%   Nucleation, avec des subtilites dans le cas ou il y a plusieurs
%   points de nucleations simultanes.
%      - Lors de la seconde etape, les clusters grossissent car de la   
%   vapeur d'eau arrive par diffusion sur les bords des clusters. Le
%   rythme de croissance des clusters est controle par le flux de vapeur
%   qui arrive par diffusion. On procede donc a un envahissement selon un 
%   algorithme temporel tenant compte d'un bilan de matiere d'eau. La 
%   direction de croissance des clusters est controle par les forces
%   capillaires : les liens ayant la pression capillaire critique la plus
%   faible sont envahis en premier. Plus de détails dans la fonction 
%   ClusterGrowth.
%   
%Input : network, options
%     network : pore network on which the algorithm is run
%     options.TemperatureInletLinks : links inlet for heat
%     options.TemperatureOutletLinks : links outlet for heat
%     options.TemperatureInlet : temperature inlet (in K)
%     options.TemperatureOutlet : temperature outlet (in K)
%     options.TemperaturePoreHeatConductivity 
%     options.LiquidWaterOutletLinks : links outlet for water
%     options.AirPressure : uniform pressure imposed to air
%     options.VaporInletLinks
%     options.VaporOutletLinks
%     options.VaporPoreDiffusivity
%     options.RelativeHumidityInlet
%     options.RelativeHumidityOutlet
%     options.ClusterOptions : options for the beheavior of liquid phase (see ClusterMonophasique)
%     
%Output : [cluster,outputInformation]
%     cluster : liquid water cluster resulting from condensation
%     outputInformation : information to analyse the degradation process 
%   
%Example :
%     [condensationInfos,condensationClusters]=ComputeCondensation(network, options)
%---------------------------------------------------------------------------------------------    

    %% Initialisation of the algorithm
    disp('Running condensation');
    tic;
    Condensation_CheckInputs(network, options)

    nPore = network.GetNumberOfPores;

    % Parametrisation of diffusion in the network 
    
    diffusionParams=struct;
    diffusionParams.GeometricModel.Pore = 'Cylinder' ;
    diffusionParams.GeometricModel.Link = 'None' ;% 'SurfaceResistance_RealSurface'
    diffusionParams.PoreBulkProp = options.VaporPoreDiffusivity ; %diffusivity = 2e-5; % 02 in N2 at ambiant conditions;diffusivity*ones(nPore,1);
    diffusionParams.LinkBulkProp =0; % scalar or array(nLink,1)
    diffusionConductances = LocalScaleComputeConductancesDiffusion(network,diffusionParams);  
        % change this to multicomponents diffusion conductance ?


    %% Compute equilibrum vapor pressure field

    %Temperature field
    [temperature,heatTransferCoefficient] = Condensation_ComputeTemperatureField(network,...
            options.TemperatureInlet,options.TemperatureOutlet,...
            options.TemperatureInletLinks,options.TemperatureOutletLinks,options.TemperatureTransportPores,options.TemperaturePoreHeatConductivity) ;
    condensationInfos.TemperatureField = temperature;
    condensationInfos.HeatTransferCoefficient = heatTransferCoefficient;

    %Equilibrum vapor pressure field
    equilibriumVaporPressure = Condensation_ComputeEquilibriumVaporPressure(temperature);
    condensationInfos.EquilibriumVaporPressure = equilibriumVaporPressure;


    %% Nucleation step
    disp('Nucleation')
    [ nucleationClusters, nucleationInfos ] = Condensation_Nucleation( network, ...
                                equilibriumVaporPressure,options,diffusionConductances,temperature );
                            
    condensationInfos.NucleationInfos = nucleationInfos;


    %% DiffusionControledCondensation
    disp('DiffusionControledCondensation')
    [condensationClusters, clusterGrowthInfos] = Condensation_ClusterGrowth(network,...
                         nucleationClusters, options, diffusionConductances,equilibriumVaporPressure,temperature);
                     
    condensationInfos.ClusterGrowthInfos = clusterGrowthInfos;

    %% Conclusion de l'algorithme                

    %TODO : update outputInformation avec condensationInfos

    duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
    fprintf('Calcul de condensation termine. Duree : %d minutes %f s.',minutes,secondes); 

end
