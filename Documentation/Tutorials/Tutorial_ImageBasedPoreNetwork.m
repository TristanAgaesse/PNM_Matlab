%Tutorial : using Image based pore networks

%First load the mat file which contains the results from python pore
%network extraction
%load(extractionResultFile)

%Create the pore network object
InputScriptImageBasedPNM
network=CreateImageBasedPoreNetwork(inputContainerMap);


phaseConductances=ones(1,256);
phaseConductances(1)=1;
phaseConductances(101)=1;
phaseConductances(201)=1;
diffusionTensor = ComputeDiffusionTensor( network,phaseConductances);
