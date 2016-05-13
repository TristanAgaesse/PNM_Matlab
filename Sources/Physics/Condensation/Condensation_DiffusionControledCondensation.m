function [condensationClusters, condensationInfos] = Condensation_DiffusionControledCondensation(...
                network,nucleationClusters, options, diffusionConductances)
%CONDENSATION_DIFFUSIONCONTROLEDCONDENSATION Deuxieme partie de l'algorithme 
% de condensation de Benjamen Straubhaar et Marc Prat. 
%   
%   Lors de cette seconde etape, la condensation se poursuit sur les 
%   bords des clusters d'eau nuclees car de la vapeur d'eau arrive par 
%   diffusion. Le  rythme de croissance des clusters est controle par le 
%   flux de vapeur qui arrive par diffusion. On procede donc a un 
%   envahissement selon un algorithme temporel tenant compte d'un bilan de 
%   matiere d'eau. La direction de croissance des clusters est controle par
%   les forces capillaires : les liens ayant la pression capillaire 
%   critique la plus faible sont envahis en premier. Plus de détails dans 
%   la fonction DiffusionControledCondensation.


    %Begin invasion loop
    %nPoreAccessible=FindNumberOfAccessiblePores(network,1:network.);
    outlet_reached = false;
    outletPores = network.GetPoresFrontiere(options.LiquidWaterOutletLinks);
    iteration = 0;

    while not(outlet_reached) %&& iteration<nPoreAccessible
        iteration = iteration+1;

        %Update partial pressure field

        partialVaporPressure = Condensation_ComputePartialVaporPressure(network,...
            cluster.GetInvadedPoresComplementary,...
            diffusionConductances,inletVaporPressure,outletVaporPressure,...
            options.VaporInletLinks,options.VaporOutletLinks,options.AirPressure,temperature);

        outputInformation.PartialVaporPressure{end+1} = partialVaporPressure;


        %new invasion : des pores condensables proches d'une zone envahie, 
        %choisir celui qui peut être envahi par IP (Pc la plus faible)  
        [~,indexMinPressureLink] = Condensation_FindNextInvadedLink(cluster,...
                            partialVaporPressure,equilibriumVaporPressure);

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
                breakthroughLinks = intersect(cluster.Network.GetLinksOfPore(invadedPore),...
                                              options.LiquidWaterOutletLinks);
                cluster.InvadeOutletLink(breakthroughLinks);
            end
        else
            %no new pore condensable near cluster
            return
        end
    end

end

