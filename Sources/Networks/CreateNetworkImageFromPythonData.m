function network=CreateNetworkImageFromPythonData(fileName,voxelEdgeLength)
%CreateNetworkImageFromPythonData Build the matlab pore network object from the
%data extracted from an image by python VirtualMaterial
%Input : fileName,voxelEdgeLength
%Output : network


    load(fileName)
    
    InputScriptImageBasedPNM
    inputContainerMap('VoxelEdgeLength')=voxelEdgeLength;
    
    network = CreateImageBasedPoreNetwork(inputContainerMap);

end