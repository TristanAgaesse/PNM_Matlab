function [ totalSaturation, saturationProfile ]=ComputeSaturationOnImage(image,nPointCurve,axe,codeForLiquid,codeForSolid)
%ComputeSaturationOnImage Computes the saturation profile on a image
% REMARK : Use PostprocessorVolumeAverage.PlotProfile for improved service.
%Input : image,nPointCurve,axe,codeForLiquid,codeForSolid
%       -image
%       -nPointCurve
%       -axe= '[0 0 1]', '[0 1 0]' or '[1 0 0]'
%       -codeForLiquid,codeForSolid
%Output : [ totalSaturation, saturationProfile ]
    

    %Computing total saturation
    totalLiquidVolume=nnz(image==codeForLiquid);
    totalSolidVolume=nnz(image==codeForSolid);
    totalVolume= numel(image) ;
    totalSaturation=totalLiquidVolume/(totalVolume-totalSolidVolume);

    %Computing saturation profile
    saturationProfile=zeros(nPointCurve,2);
    
    imThickness=size(image,find(axe));
    sliceThickness = floor(imThickness/nPointCurve);
    
    for iPointCurve=1:nPointCurve
        
        indices=sliceThickness*(iPointCurve-1)+1:sliceThickness*iPointCurve;
        
        if isequal(axe,[0 0 1])
            slice=image(:,:,indices);
            
        elseif isequal(axe,[0 1 0])
            slice=image(:,indices,:);
            
        elseif isequal(axe,[1 0 0])
            slice=image(indices,:,:);
        end

        liquidVolume = nnz(slice==codeForLiquid);
        solidVolume = nnz(slice==codeForSolid);
        sliceVolume = numel(slice) ;

        saturationProfile(iPointCurve,1)=iPointCurve/nPointCurve;
        saturationProfile(iPointCurve,2)=liquidVolume/(sliceVolume-solidVolume);
        
    end
    
end
