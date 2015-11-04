function network=GetNetworkFromPythonData(fileName)
%GetNetworkFromPythonData Build the matlab pore network object from the
%data extracted from an image by python VirtualMaterial
%Input : fileName
%Output : network


    load(fileName)
    InputScriptImageBasedPNM
    network = CreateImageBasedPoreNetwork(inputContainerMap);

end