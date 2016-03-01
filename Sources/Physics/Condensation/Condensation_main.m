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
CheckInputs(network, options)

%Compute temperature field

[temperature,heatTransferCoefficient] = ComputeTemperatureField(network,options.TemperatureInlet,options.TemperatureOutlet,options.TemperatureInletLinks,options.TemperatureOutletLinks) ;

outputInformation.TemperatureField = temperature;
outputInformation.HeatTransferCoefficient = heatTransferCoefficient;



%Compute equilibrum vapor pressure field

equilibriumVaporPressure = ComputeEquilibriumVaporPressure(temperature);

outputInformation.EquilibriumVaporPressure = equilibriumVaporPressure;



%Compute partial pressure field
gasTransportPores = 1:network.GetNumberOfPores;

inletVaporPressure = options.RelativeHumidityInlet*options.AirPressure;
outletVaporPressure = options.RelativeHumidityOutlet*options.AirPressure;

diffusivity = 2e-5; % 02 in N2 at ambiant conditions TODO : check value for water vapor
diffusionConductances = LocalScaleComputeConductancesDiffusion(network,diffusivity);  %TODO : change this to multicomponents diffusion conductance

partialVaporPressure = ComputePartialVaporPressure(network,gasTransportPores,diffusionConductances,inletVaporPressure,outletVaporPressure,options.VaporInletLinks,options.VaporOutletLinks,options.AirPressure,temperature);

outputInformation.PartialVaporPressure{1} = partialVaporPressure;



%invade the pore which has the max partial pressure if > equilibrium
%partial pressure

condensationRatio = partialVaporPressure ./ equilibriumVaporPressure ;
[maxRatio,indexMaxRatio] = max(condensationRatio);

if maxRatio>1
    firstInvadedPore=indexMaxRatio;
else
    return
end

cluster = ClusterMonophasique.InitialiseCondensationCluster(network,options.ClusterOptions,firstInvadedPore,options.LiquidWaterOutletLinks);

outputInformation.InvadedPore{1}=firstInvadedPore;



%Begin invasion loop
%nPoreAccessible=FindNumberOfAccessiblePores(network,1:network.);
outlet_reached = false;
outletPores = network.GetPoresFrontiere(options.LiquidWaterOutletLinks);
iteration = 0;

while not(outlet_reached) %&& iteration<nPoreAccessible
    iteration = iteration+1;

    %Update partial pressure field

    partialVaporPressure = ComputePartialVaporPressure(network,cluster.GetInvadedPoresComplementary,diffusionConductances,inletVaporPressure,outletVaporPressure,options.VaporInletLinks,options.VaporOutletLinks,options.AirPressure,temperature);

    outputInformation.PartialVaporPressure{end+1}= partialVaporPressure;


    %new invasion : des pores condensables proches d'une zone envahie, 
    %choisir celui qui peut être envahi par IP (Pc la plus faible)  
    [~,indexMinPressureLink] = FindNextInvadedLink(cluster,partialVaporPressure,equilibriumVaporPressure);

    if indexMinPressureLink>0
        invadedPore = cluster.GetOutwardPore(indexMinPressureLink);
        outputInformation.InvadedPore{end+1} = invadedPore;

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


    
    
    