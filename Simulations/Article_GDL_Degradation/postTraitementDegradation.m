function infos=postTraitementDegradation(network,outputInformation)
    

    infos=struct;
    nStep=outputInformation.nIteration;
    
    infos.SimuOptions=outputInformation.options;
    infos.FixedPressure=outputInformation.fixedPressure;
    
    inletLink=outputInformation.inletLink;
    outletLink=outputInformation.outletLink;
    
    disp('Post-Traitement : Times')
    infos.Times=outputInformation.times(1:nStep);

    disp('Post-Traitement : Pressure')
    pressures=zeros(1,sum(cellfun(@length,outputInformation.invasionPressures)));
    indice=1;
    for iStep=1:nStep
        pressures(indice:(indice+length(outputInformation.invasionPressures{iStep})-1))=outputInformation.invasionPressures{iStep};
        indice=indice+length(outputInformation.invasionPressures{iStep});
    end
    infos.Pressure=pressures;
    
    disp('Post-Traitement : nInvasionPerStep')
    infos.nInvasionPerStep=cellfun('length',outputInformation.invasionPressures);
    
    disp('Post-Traitement : Saturation')
    saturation=zeros(1,nStep);
    for iStep=1:nStep
        options.type ='totalSaturation';
        saturation(iStep)=ComputeSaturation(outputInformation.clusters{iStep},network,options);
    end
    infos.Saturation=saturation;
    
    disp('Post-Traitement : NombrePointsDePercée')
    nPointsPercee=zeros(1,nStep);
    outletLink=outputInformation.outletLink;
    for iStep=1:nStep
        foo=outputInformation.clusters{iStep}.GetInvadedLinksBoolean;
        nPointsPercee(iStep)=nnz(foo(outletLink));
    end
    infos.nPointsPercee=nPointsPercee;
    
    
    disp('Post-Traitement : PourcentageDegradationMoyennePore')
    pourcentageDegradation=cell(1,nStep);
    for iStep=1:nStep
        invadedPores=outputInformation.clusters{iStep}.GetInvadedPores;
        unprocessedDegradation=outputInformation.degradationPercentage{iStep};
        pourcentageDegradation{iStep}=zeros(1,network.GetNumberOfPores);
        for iPore=invadedPores
            poreLinks=network.GetLinksOfPore(iPore);
            pourcentageDegradation{iStep}(iPore)=sum(unprocessedDegradation(poreLinks))/length(poreLinks);
        end
    end
    infos.PourcentageDegradationMoyennePore=pourcentageDegradation;
    
    disp('Post-Traitement : PoresEnvahis')
    poreEnvahis=cell(1,nStep);
    poreEnvahis{1}=outputInformation.clusters{1}.GetInvadedPores;
    for iStep=2:nStep
        poreEnvahis{iStep}=setdiff(outputInformation.clusters{iStep}.GetInvadedPores,outputInformation.clusters{iStep-1}.GetInvadedPores);
    end
    infos.PoreEnvahis=poreEnvahis;
    
    
    disp('Post-Traitement : BurstVolume')
    
    burstVolume=zeros(1,nStep);  
    poreVolume=network.GetPoreData('Volume');
    for iStep=1:nStep
        burstVolume(iStep) = sum(poreVolume(infos.PoreEnvahis{iStep}));
    end
    
    infos.BurstVolume=burstVolume;
    
    disp('Post-Traitement : InvadedPoreVolume')
    
    nInvadedPore = length(outputInformation.clusters{outputInformation.nIteration}.GetInvadedPores);
    invadedPoreVolume=zeros(1,nInvadedPore);  
    poreVolume=network.GetPoreData('Volume');
    iInvadedPore=1;
    for iStep=1:nStep
        burstPores=infos.PoreEnvahis{iStep};
        for i=1:length(burstPores)
            invadedPoreVolume(iInvadedPore) = poreVolume(burstPores(i));
            iInvadedPore=iInvadedPore+1;
        end
    end
    assert(iInvadedPore==nInvadedPore+1)
    infos.InvadedPoreVolume=invadedPoreVolume;
    
    
    disp('Post-Traitement : InvasionTimeDistribution')
    infos.InvasionTimeDistribution=outputInformation.invasionTimeDistribution;
    
    
    disp('Post-Traitement : MeanInvasionTime')
    meanInvasionTime=zeros(1,nStep);
    for iStep=1:nStep
        meanInvasionTime(iStep)=mean(outputInformation.invasionTimeDistribution{iStep});
    end
    infos.MeanInvasionTime = meanInvasionTime;
    
    
    disp('Post-Traitement : ContactAngle')
    infos.ContactAngle=outputInformation.contactAngle;
    
    disp('Post-Traitement : MeanContactAngle')
    meanContactAngle=zeros(1,nStep);
    for iStep=1:nStep
        meanContactAngle(iStep)=mean(infos.ContactAngle{iStep} );
    end
    infos.MeanContactAngle=meanContactAngle;
    
    
    disp('Post-Traitement : InterfaceLinkContactAngle')
    interfaceLinkContactAngle=cell(1,nStep);
    for iStep=1:nStep
        contactAngles=outputInformation.contactAngle{iStep};
        interfaceLinkContactAngle{iStep}=contactAngles(outputInformation.clusters{iStep}.GetInterfaceLinks);
    end
    infos.InterfaceLinkContactAngle=interfaceLinkContactAngle;
    
    disp('Post-Traitement : InterfaceLinkMeanContactAngle')
    interfaceLinkMeanContactAngle=zeros(1,nStep);
    for iStep=1:nStep
        interfaceLinkMeanContactAngle(iStep)=mean(infos.InterfaceLinkContactAngle{iStep});
    end
    infos.InterfaceLinkMeanContactAngle=interfaceLinkMeanContactAngle;
    
    
    disp('Post-Traitement : InterfaceLinkDiameter')
    diameter=network.GetLinkData('Diameter');
    interfaceLinkDiameter=cell(1,nStep);
    for iStep=1:nStep
        interfaceLinkDiameter{iStep}=diameter(outputInformation.clusters{iStep}.GetInterfaceLinks);
    end
    infos.InterfaceLinkDiameter=interfaceLinkDiameter;
    
    disp('Post-Traitement : InterfaceLinkMeanDiameter')
    interfaceLinkMeanDiameter=zeros(1,nStep);
    for iStep=1:nStep
        interfaceLinkMeanDiameter(iStep)=mean(infos.InterfaceLinkDiameter{iStep}) ;
    end
    infos.InterfaceLinkMeanDiameter=interfaceLinkMeanDiameter;
    
    
    
    disp('Post-Traitement : CapillaryPressure')
    capillaryPressure=cell(1,nStep);
    for iStep=1:nStep
        [nelements,xcenters]=hist(outputInformation.clusters{iStep}.GetCriticalPressures,100);
        capillaryPressure{iStep}{1}=nelements;
        capillaryPressure{iStep}{2}=xcenters;
    end
    infos.CapillaryPressure=capillaryPressure;
    
    disp('Post-Traitement : DiffusionCoefficient')
    diffusionCoefficient=zeros(1,nStep);
    for iStep=1:nStep
         try
            transportPores = setdiff(1:network.GetNumberOfPores,outputInformation.clusters{iStep}.GetInvadedPores) ;

            diffusivity = 2e-5; % 02 in N2 at ambiant conditions
            conductancesDiffusion = LocalScaleComputeConductancesDiffusion(network,diffusivity);

            boundaryConditions.inletLink = outletLink;
            boundaryConditions.outletLink = inletLink;

            boundaryConditions.inletType = 'Dirichlet' ;
            boundaryConditions.outletType = 'Dirichlet' ;

            nInletLink = length(boundaryConditions.inletLink);
            nOutletLink = length(boundaryConditions.outletLink);

            boundaryConditions.inletValue = 1*ones(1,nInletLink);
            boundaryConditions.outletValue = 0.1*ones(1,nOutletLink);
            
            [ ~, ~, diffCoeff ] = ComputeLinearTransport( ...
                                                            network,transportPores, ...
                                                            conductancesDiffusion, ...
                                                            boundaryConditions  ...
                                                            );
        catch
            disp('Erreur lors du calcul de la diffusion');
            diffCoeff=0;
        end
        diffusionCoefficient(iStep)=diffCoeff;
    end
    infos.DiffusionCoefficient=diffusionCoefficient;
    
    
end