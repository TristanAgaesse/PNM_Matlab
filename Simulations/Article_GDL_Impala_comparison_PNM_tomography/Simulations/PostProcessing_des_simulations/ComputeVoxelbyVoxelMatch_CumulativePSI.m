function voxelByVoxelError=ComputeVoxelbyVoxelMatch_CumulativePSI(image1,image2,labels)


    nLabel = length(labels);
    voxelByVoxelError = zeros(nLabel,1);
    
    labels = sort(labels);

    [labelEnds1,orderLabels1,labelIndices1] = PoreNetworkImageBased.ParseLabeledImage(image1);
    [labelEnds2,orderLabels2,labelIndices2] = PoreNetworkImageBased.ParseLabeledImage(image2);
    
    voxels1=[];
    voxels2=[];
    
    for i = 1:nLabel
        
        iLabel = labels(i);
        incrementVoxels1 = PoreNetworkImageBased.GetVoxelsOfLabel(iLabel,labelEnds1,orderLabels1,labelIndices1);
        incrementVoxels2 = PoreNetworkImageBased.GetVoxelsOfLabel(iLabel,labelEnds2,orderLabels2,labelIndices2);
        
        voxels1=[voxels1,incrementVoxels1];
        voxels2=[voxels2,incrementVoxels2];
        
        voxelByVoxelError(i) = numel(intersect(voxels1,voxels2))/numel(voxels1);
        
    end
    
    
    
end

