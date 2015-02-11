%Install_PNM_Matlab

%Script to execute the first time you use PNM_Matlab on your computer
disp('Compiling some functions')
%Compile the mex file (C++ function) graph_conn_comp_mex.cpp
cd Sources/ExternalLibrairies/graph_conn_comp
mex -O -largeArrayDims graph_conn_comp_mex.cpp
cd ../../..
%If this compilation doesn't work, change the option 'graphconncomp' to 
%'HomeMadeTarjan' in the function PoreNetwork.FindComposantesConnexes
