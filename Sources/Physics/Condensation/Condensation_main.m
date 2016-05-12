%COMPUTECONDENSATION Summary of this function goes here
%   Detailed explanation goes here
%
%Algorithm : - Demarrage de la condensation dans le pore ayant l'humidité
%                   relative la plus élevée
%            - parmis les pores condensables proches du cluster de condensation, 
%                   envahir celui qui peut être envahi avec la pressions 
%                   capillaire la plus faible  
%
%Input : network, options
%       network : pore network on which the algorithm is run
%       options.TemperatureInletLinks : links inlet for heat
%       options.TemperatureOutletLinks : links outlet for heat
%       options.TemperatureInlet : temperature inlet (in K)
%       options.TemperatureOutlet : temperature outlet (in K)
%       options.LiquidWaterOutletLinks : links outlet for water
%       options.AirPressure : uniform pressure imposed to air
%       options.VaporInletLinks
%       options.VaporOutletLinks
%       options.RelativeHumidityInlet
%       options.RelativeHumidityOutlet
%       options.ClusterOptions : options for the beheavior of liquid phase (see ClusterMonophasique)
%                                  
%Output : [cluster,outputInformation]
%       cluster : liquid water cluster resulting from condensation
%       outputInformation : information to analyse the degradation process 

%---------------------------------------------------------------------------------------------    


disp('Running condensation');
tic;
Condensation_CheckInputs(network, options)

%Compute temperature field

[temperature,heatTransferCoefficient] = Condensation_ComputeTemperatureField(network,...
        options.TemperatureInlet,options.TemperatureOutlet,...
        options.TemperatureInletLinks,options.TemperatureOutletLinks) ;

outputInformation.TemperatureField = temperature;
outputInformation.HeatTransferCoefficient = heatTransferCoefficient;



%Compute equilibrum vapor pressure field
equilibriumVaporPressure = Condensation_ComputeEquilibriumVaporPressure(temperature);
outputInformation.EquilibriumVaporPressure = equilibriumVaporPressure;


%Nucleation step
nPore = network.GetNumberOfPores;

diffusivity = 2e-5; % 02 in N2 at ambiant conditions TODO : check value for water vapor
diffusionParams=struct;
diffusionParams.GeometricModel.Pore = 'Cylinder' ;
diffusionParams.GeometricModel.Link = 'None' ;% 'SurfaceResistance_RealSurface'
diffusionParams.PoreBulkProp = diffusivity*ones(nPore,1);
diffusionParams.LinkBulkProp =0; % scalar or array(nLink,1)

diffusionConductances = LocalScaleComputeConductancesDiffusion(network,diffusionParams);  
%TODO : change this to multicomponents diffusion conductance


[ nucleationClusters, nucleationInfos ] = Condensation_Nucleation( network, ...
                            equilibriumVaporPressure,options,diffusionConductances );
outputInformation.PartialVaporPressure{1} = nucleationInfos.PartialVaporPressure{1};
outputInformation.InvadedPore{1}=firstInvadedPore;



%Begin invasion loop
%nPoreAccessible=FindNumberOfAccessiblePores(network,1:network.);
outlet_reached = false;
outletPores = network.GetPoresFrontiere(options.LiquidWaterOutletLinks);
iteration = 0;

while not(outlet_reached) %&& iteration<nPoreAccessible
    iteration = iteration+1;

    %Update partial pressure field

    partialVaporPressure = Condensation_ComputePartialVaporPressure(network,cluster.GetInvadedPoresComplementary,...
        diffusionConductances,inletVaporPressure,outletVaporPressure,...
        options.VaporInletLinks,options.VaporOutletLinks,options.AirPressure,temperature);

    outputInformation.PartialVaporPressure{end+1} = partialVaporPressure;


    %new invasion : des pores condensables proches d'une zone envahie, 
    %choisir celui qui peut être envahi par IP (Pc la plus faible)  
    [~,indexMinPressureLink] = Condensation_FindNextInvadedLink(cluster,partialVaporPressure,equilibriumVaporPressure);

    if indexMinPressureLink>0
        invadedPore = cluster.GetOutwardPore(indexMinPressureLink);
        outputInformation.InvadedPore{end+1} = invadedPore;

        %TODO : Gérer l'évolution du cluster
        
            %TODO : envahir les pores actifs des autres clusters en fonction du temps et des débits
        
        interfaceChangeInformation=cluster.InvadeNewPore(indexMinPressureLink);
        cluster.UpdateCriticalPressure(interfaceChangeInformation,[],options.LiquidWaterOutletLinks);


        %verifier si outlet_reached
        if ismember(invadedPore,outletPores)
            outlet_reached = true;
            breakthroughLinks = intersect(cluster.Network.GetLinksOfPore(invadedPore),options.LiquidWaterOutletLinks);
            cluster.InvadeOutletLink(breakthroughLinks);
        end
    else
        %no new pore condensable near cluster
        return
    end
end

duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
fprintf('Calcul de condensation termine. Duree : %d minutes %f s.',minutes,secondes); 


    
    
    