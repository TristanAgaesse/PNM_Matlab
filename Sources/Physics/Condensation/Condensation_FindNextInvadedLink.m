% function [minPressure,indexMinPressureLink] = Condensation_FindNextInvadedLink(cluster,partialVaporPressure,equilibriumVaporPressure)
%     %new invasion : des pores condensables proches d'une zone envahie, 
%     %choisir celui qui peut être envahi par IP (Pc la plus faible)  
%     
%     %TODO : add invasion time
%     
%     %List of condensable pores next to cluster
%     poresCondensables = partialVaporPressure>equilibriumVaporPressure ;
%     
%     outwardPores = cluster.GetOutwardPore(1:length(cluster.GetInterfaceLinks));
%     hasVoisinCondensable = poresCondensables(outwardPores(outwardPores>0));
%     clusterLinkCondensable = find(hasVoisinCondensable);
% 
%     %Chose condensable pore accessible with min capillary pressure
%     if not(isempty(clusterLinkCondensable))
%         criticalPressures  =  cluster.GetCriticalPressures;
%         boundaryCriticalPressures = criticalPressures(cluster.GetInterfaceLinks);
%         [minPressure,indexMinPressureLink]  =  min(boundaryCriticalPressures(clusterLinkCondensable));
% 
%         indicesPoresCondensables = find(poresCondensables);
%         indexMinPressureLink=indicesPoresCondensables(indexMinPressureLink);
%     else
%         indexMinPressureLink=-1;
%         minPressure=-1;
%     end
% end