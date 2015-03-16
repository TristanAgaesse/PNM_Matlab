function infos=postTraitementDegradation(network,floodingStepInformation)
    

    infos=struct;
    nStep=length(floodingStepInformation.times);
    
    infos.SimuOptions=floodingStepInformation.options;
    infos.FixedPressure=floodingStepInformation.fixedPressure;
    
    inletLink=floodingStepInformation.inletLink;
    outletLink=floodingStepInformation.outletLink;
    
    disp('Post-Traitement : Times')
    infos.Times=floodingStepInformation.times(1:nStep);

    disp('Post-Traitement : Pressure')
    pressures=zeros(1,sum(cellfun(@length,floodingStepInformation.invasionPressures)));
    indice=1;
    for iStep=1:nStep
        pressures(indice:(indice+length(floodingStepInformation.invasionPressures{iStep})-1))=floodingStepInformation.invasionPressures{iStep};
        indice=indice+length(floodingStepInformation.invasionPressures{iStep});
    end
    infos.Pressure=pressures;
    
    disp('Post-Traitement : nInvasionPerStep')
    infos.nInvasionPerStep=cellfun('length',floodingStepInformation.invasionPressures);
    
    disp('Post-Traitement : Saturation')
    saturation=zeros(1,nStep);
    for iStep=1:nStep
        options.type ='totalSaturation';
        saturation(iStep)=ComputeSaturation(floodingStepInformation.clusters{iStep},network,options);
    end
    infos.Saturation=saturation;
    
    disp('Post-Traitement : NombrePointsDePercée')
    nPointsPercee=zeros(1,nStep);
    outletLink=floodingStepInformation.outletLink;
    for iStep=1:nStep
        foo=floodingStepInformation.clusters{iStep}.GetInvadedLinksBoolean;
        nPointsPercee(iStep)=nnz(foo(outletLink));
    end
    infos.nPointsPercee=nPointsPercee;
    
    
    disp('Post-Traitement : PourcentageDegradationMoyennePore')
    pourcentageDegradation=cell(1,nStep);
    for iStep=1:nStep
        invadedPores=floodingStepInformation.clusters{iStep}.GetInvadedPores;
        unprocessedDegradation=floodingStepInformation.degradationPercentage{iStep};
        pourcentageDegradation{iStep}=zeros(1,network.GetNumberOfPores);
        for iPore=invadedPores
            poreLinks=network.GetLinksOfPore(iPore);
            pourcentageDegradation{iStep}(iPore)=sum(unprocessedDegradation(poreLinks))/length(poreLinks);
        end
    end
    infos.PourcentageDegradationMoyennePore=pourcentageDegradation;
    
    disp('Post-Traitement : PoresEnvahis')
    poreEnvahis=cell(1,nStep);
    poreEnvahis{1}=floodingStepInformation.clusters{1}.GetInvadedPores;
    for iStep=2:nStep
        poreEnvahis{iStep}=setdiff(floodingStepInformation.clusters{iStep}.GetInvadedPores,poreEnvahis{iStep-1});
    end
    infos.PoreEnvahis=poreEnvahis;
    
    
    disp('Post-Traitement : BurstSize')
    
    % TODO  
    
    
    disp('Post-Traitement : InvasionTimeDistribution')
    
    % TODO  
    
    
    disp('Post-Traitement : CapillaryPressure')
    capillaryPressure=cell(1,nStep);
    for iStep=1:nStep
        [nelements,xcenters]=hist(floodingStepInformation.clusters{iStep}.GetCriticalPressures,100);
        capillaryPressure{iStep}{1}=nelements;
        capillaryPressure{iStep}{2}=xcenters;
    end
    infos.CapillaryPressure=capillaryPressure;
    
    disp('Post-Traitement : DiffusionCoefficient')
    diffusionCoefficient=zeros(1,nStep);
    for iStep=1:nStep
         try
            [ ~, ~, ~, diffCoeff ]=ComputeDiffusion(network, floodingStepInformation.clusters{iStep}.GetComplementaryCluster,outletLink, inletLink);
        catch
            disp('Erreur lors du calcul de la diffusion');
            diffCoeff=0;
        end
        diffusionCoefficient(iStep)=diffCoeff;
    end
    infos.DiffusionCoefficient=diffusionCoefficient;
    
    
end