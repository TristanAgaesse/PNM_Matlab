function [ network,viewer ] = CreateNetwork( geometryFileName )
%CREATENETWORK Cree un reseau de pores et un viewer a partir d'un fichier de
%geometrie macroscopique. Cette fonction est l'equivalent en ligne de
%commande de UserInterface.
%   Input : geometryFileName : nom du fichier de geometrie macroscopique
%   Output : [ network,viewer ]
%               - network : objet reseau de pores correspondant au fichier
%               geometrie macroscopique
%               - viewer : objet viewer correspondant au network


    myGeometry=MacroscopicGeometry();
    myGeometry.LoadGeometry(geometryFileName);
    networkBuilder=NetworkBuilder(myGeometry);
    
    network=networkBuilder.BuildNetwork;

    viewer=Viewer(network.PrivateInternalOutputStruct);
end

