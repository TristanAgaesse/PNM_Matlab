
%% Scenario definition
scenario = 2;

if scenario == 1
    %Scenarios profil quadratique d'angle de contact
    minWettability=[100,100,100,100,100,100];
    maxWettability=[100,110,120,130,140,150];
    
    min=minWettability*pi/180;
    max=maxWettability*pi/180;
    
    plotAbsisse = maxWettability-minWettability;
    
    nParam=length(min);
    contactAngleOption=cell(1,nParam);
    for i=1:nParam
        contactAngleOption{i}.Type='quadratique';
        contactAngleOption{i}.MinContactAngle = min(i)  ;
        contactAngleOption{i}.MaxContactAngle = max(i) ;
    end

elseif scenario == 2 
    %Scenarios profil quadratique d'angle de contact
    minWettability=[130,120,110,100,95];
    maxWettability=[100,110,120,130,135];
    
    min=minWettability*pi/180;
    max=maxWettability*pi/180;
    
    plotAbsisse = maxWettability-minWettability;
    
    nParam=length(min);
    contactAngleOption=cell(1,nParam);
    for i=1:nParam
        contactAngleOption{i}.Type='quadratique';
        contactAngleOption{i}.MinContactAngle = min(i)  ;
        contactAngleOption{i}.MaxContactAngle = max(i) ;
    end
    
elseif scenario == 3
    %Scenarios PTFE_sur_petite_epaisseur
    minWettability=[130,120,120,120,110,110,110];
    maxWettability=[130,140,140,140,150,150,150];
    epaisseurRelative=[0.1,0.05,0.1,0.2,0.05,0.1,0.2];

    min=minWettability*pi/180;
    max=maxWettability*pi/180;

    nParam=length(min);
    plotAbsisse = 1:nParam;

    contactAngleOption=cell(1,nParam);
    for i=1:nParam
        contactAngleOption{i}.Type='PTFE_sur_petite_epaisseur';
        contactAngleOption{i}.MinContactAngle = min(i)  ;
        contactAngleOption{i}.MaxContactAngle = max(i) ;
        contactAngleOption{i}.EpaisseurRelative = epaisseurRelative(i);
    end

end

%% Simulations

diffusionCoefficient=zeros(1,nParam);
totalSaturation=zeros(1,nParam);
saturationProfile=cell(1,nParam);
cluster=cell(1,nParam);
concentrations=cell(1,nParam);

for i=1:nParam
        
    [diffCoeff,totalSat, satProfile,invCluster,concentr] = RunSimu_NonUniformWettability( network,contactAngleOption{i} );
    diffusionCoefficient(i)=diffCoeff;
    totalSaturation(i)=totalSat;
    saturationProfile{i}=satProfile;
    cluster{i}=invCluster;
    concentrations{i}=concentr;
end


%% Plot figures
figure
plot(plotAbsisse,diffusionCoefficient)
title('Relative Diffusion Coefficient')
xlabel('Ecart de mouillabilite entre la surface et l''interieur de la GDL')
ylabel('Relative Diffusion Coefficient')

figure
plot(plotAbsisse,totalSaturation)
title('Total saturation')
xlabel('Ecart de mouillabilite entre la surface et l''interieur de la GDL')
ylabel('Total saturation')

figure;
color=colormap(hsv(nParam));
for i=1:nParam
    satProf=saturationProfile{i};
    plot(satProf(:,1),satProf(:,2),'Color',color(i,:))
    hold on
end
title('Saturation profiles')
xlabel('Position dans l''epaisseur')
ylabel('Slice saturation')

