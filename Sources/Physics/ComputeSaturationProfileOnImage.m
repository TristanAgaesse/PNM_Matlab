function saturationProfile=ComputeSaturationProfileOnImage(image,nPointCurve,axe,codeForLiquid,codeForVoid)
%ComputeSaturationProfileOnImage Computes the saturation profile on a image
%Input : image,nPointCurve,axe,codeForLiquid,codeForVoid
%       -image
%       -nPointCurve
%       -axe= '[0 0 1]', '[0 1 0]' or '[1 0 0]'
%       -codeForLiquid,codeForVoid
%Output : saturationProfile
    
    saturationProfile=zeros(nPointCurve,2);
    
    imThickness=size(image,find(axe));
    sliceThickness = floor(imThickness/nPointCurve);
    
    for iPointCurve=1:nPointCurve
        
        indices=sliceThickness*(iPointCurve-1):sliceThickness*iPointCurve;
        indices=indices+1;
        
        if isequal(axe,[0 0 1])
            slice=image(:,:,indices);
            
        elseif isequal(axe,[0 1 0])
            slice=image(:,indices,:);
            
        elseif isequal(axe,[1 0 0])
            slice=image(indices,:,:);
        end

        liquidVolume=length(find(slice==codeForLiquid));
        voidVolume=length(find(slice==codeForVoid));

        saturationProfile(iPointCurve,1)=iPointCurve/nPointCurve;
        saturationProfile(iPointCurve,2)=liquidVolume/(liquidVolume+voidVolume);
        
    end
    
end
