function [figSatProf,figCapPres]=ProcessSaturationCurves(image,pressureStep,pressureCode,axe,codeForSolid,titleParameter)


    nPressure=length(pressureStep);

    color=colormap(hsv(nPressure));

    imageShape=size(image);
    nSlice = min(imageShape(axe>0));
    codeForLiquid = 100;
    
    saturation=zeros(1,nPressure);
    
    [labelEnds,orderLabels,labelIndices] = PoreNetworkImageBased.ParseLabeledImage(image);
    
    
    
    figure  
    for iPressure=1:nPressure
        linearIndices=PoreNetworkImageBased.GetVoxelsOfLabel(pressureCode(iPressure),labelEnds,orderLabels,labelIndices);
        image( linearIndices) = 100;
        
        [ totalSaturation, saturationProfile ] = ComputeSaturationOnImage(image,nSlice,axe,codeForLiquid,codeForSolid) ;
        saturation(iPressure) = totalSaturation;
        plot(saturationProfile(:,1),saturationProfile(:,2),'Color',color(iPressure,:));
        leg{iPressure} = sprintf('Pore network %d mbar', round(pressureStep(iPressure)/100));
        hold on;
    end
    figSatProf=gcf;
    title(strcat('Saturation profiles ',titleParameter));
    xlabel('Relative thickness');
    ylabel('Slice saturation');
    legend(leg)
    
    figure
    plot(pressureStep,saturation)
    legend('Pore network ');
    title(strcat('Capillary Pressure Curve ',titleParameter));
    xlabel('Pressure');
    ylabel('Saturation');
    figCapPres=gcf ;
end

