function voxelByVoxelError=ComputeVoxelbyVoxelMatch(image1,image2,labels)

    nLabel = length(labels);
    voxelByVoxelError = zeros(nLabel,1);

    [labelEnds1,orderLabels1,labelIndices1] = PoreNetworkImageBased.ParseLabeledImage(image1);
    [labelEnds2,orderLabels2,labelIndices2] = PoreNetworkImageBased.ParseLabeledImage(image2);
    
    for i = 1:nLabel
        
        iLabel = labels(i);
        voxels1 = PoreNetworkImageBased.GetVoxelsOfLabel(iLabel,labelEnds1,orderLabels1,labelIndices1);
        voxels2 = PoreNetworkImageBased.GetVoxelsOfLabel(iLabel,labelEnds2,orderLabels2,labelIndices2);
        
        voxelByVoxelError(i) = numel(intersect(voxels1,voxels2))/numel(voxels1);
        
    end
    
    
    
end

