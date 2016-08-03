%Tutorial : using Image based pore networks

%Create the pore network object from the mat file which contains the 
% results from pore network extraction

voxelEdgeLength=1e-6; % if voxel size is 1 micrometer
network = CreateNetworkImageFromPythonData(extractionFileName);

phaseConductances=ones(1,256);
phaseConductances(1)=1;
phaseConductances(101)=1;
phaseConductances(201)=1;
diffusionTensor = ComputeDiffusionTensor( network,phaseConductances);
