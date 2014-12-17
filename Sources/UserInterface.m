%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%       Matlab Pore Network Modeling      %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%           Outil de maillage             %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description : interface interactive pour construire un reseau 
%de pores. Prend un fichier de geometrie macroscopique et construit un
%fichier reseau de pore adape a cette geometrie. Ce script est un wrapper
%pour la fonction createNetwork
%



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
disp('Bonjour, c''est le Mailleur qui vous parle.');

%demande a l'utilisateur  le nom du fichier de geometrie macroscopique
file_name = input('Entrez un nom de fichier pour la geometrie macroscopique s''il vous plait : \n');

%construire le reseau de pores a partir du fichier de geometrie macroscopique
[ network,viewer ] = CreateNetwork( file_name ) ;









