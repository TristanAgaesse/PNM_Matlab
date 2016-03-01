function nPoreAccessible= Condensation_FindNumberOfAccessiblePores(network,inletLink)

    fooCluster=network.CreateVoidCluster;
    totalFloodCluster=fooCluster.GetComplementaryCluster;
    percoPath=totalFloodCluster.FindPercolationPath(inletLink,1:network.GetNumberOfLinks);
    nPoreAccessible=0;
    for i=1:length(percoPath)
        nPoreAccessible=nPoreAccessible+length(percoPath{i}.GetInvadedPores);
    end

end