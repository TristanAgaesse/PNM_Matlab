# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 16:41:14 2015

@author: ta240184
"""
import matplotlib.pylab as _plt

import hdf5storage

folder ='/home/270.12-Modeling_PEMFC_Li/theseTristan/PSI_24BA/ResultsPNM/extractionResults/'

data=hdf5storage.loadmat(folder+'GDLonly_rassembledSizeDistributions.mat')


poreDiams=[];
poreDiams.append(data['poreDiam_AllMax'])
poreDiams.append(data['poreDiam_LM10'])
poreDiams.append(data['poreDiam_LM20'])
poreDiams.append(data['poreDiam_LM4'])
poreDiams.append(data['poreDiam_h4'])
poreDiams.append(data['poreDiam_h8'])
poreDiams.append(data['poreDiam_h12'])
poreDiams.append(data['poreDiam_h20'])

linkRadius=[];
linkRadius.append(data['linkDiam_AllMax'])
linkRadius.append(data['linkDiam_LM10'])
linkRadius.append(data['linkDiam_LM20'])
linkRadius.append(data['linkDiam_LM4'])
linkRadius.append(data['linkDiam_h4'])
linkRadius.append(data['linkDiam_h8'])
linkRadius.append(data['linkDiam_h12'])
linkRadius.append(data['linkDiam_h20'])

coordination=[];
coordination.append(data['coordination_AllMax'])
coordination.append(data['coordination_LM10'])
coordination.append(data['coordination_LM20'])
coordination.append(data['coordination_LM4'])
coordination.append(data['coordination_h4'])
coordination.append(data['coordination_h8'])
coordination.append(data['coordination_h12'])
coordination.append(data['coordination_h20'])




def ComputeCaracteristics(poreDiameters,linkRadii,poreVoisins):

    nPore = poreDiameters.size
    nLink = linkRadii.size
    
#    poreDiameters=poreDiameters(poreDiameters<2e-4);
    #[poreDiamBins,poreDiamBinCenters] = hist(poreDiameters,25);
    
    edges = transpose(linspace(0,2e-4,25));
    #edges = linspace(min(poreDiameters),max(poreDiameters),25);
    
    
    dp = network['pore.diameter']
    Vp = network['pore.volume']
    dt = network['throat.diameter']
    Vt = network['throat.volume']
    dmax = max(max(dp), max(dt))
    steps = _sp.linspace(0, dmax, 100, endpoint=True)
    vals = _sp.zeros_like(steps)
    for i in range(0, len(steps)-1):
        temp1 = dp > steps[i]
        temp2 = dt > steps[i]
        vals[i] = sum(Vp[temp1]) + sum(Vt[temp2])
    
    
    
    [poreDiamBinsOccurrences,bin] = histc(poreDiameters,edges);
    poreDiamBinCenters = (edges(1:end-1)+edges(2:end))/2;
    #volume=network.GetPoreData('Volume');
    volume=(4/3)*pi*(poreDiameters/2).^3;
    poreDiamBinsVolume=zeros(length(edges),1);
    for iPore in range(length(poreDiameters)):
        poreDiamBinsVolume(bin(iPore)) = poreDiamBinsVolume(bin(iPore))+volume(iPore);

    totalPoreVolume=sum(volume);
    poreDiamBinsVolume=poreDiamBinsVolume./totalPoreVolume;
    poreDiamBinsVolume=poreDiamBinsVolume(1:end-1);
    poreDiamBinsOccurrences=poreDiamBinsOccurrences(1:end-1);
    poreDiamData=[poreDiamBinsOccurrences,poreDiamBinsVolume,poreDiamBinCenters]
    
    
    
    linkDiameters=2*linkRadii
    linkDiameters=linkDiameters(linkDiameters<2e-4);
    [linkDiamBins,linkDiamBinCenters] = hist(linkDiameters,25);
    linkDiamData=[linkDiamBins,linkDiamBinCenters]
    
 
 
    poreVoisins=poreVoisins[poreVoisins<=25];
    [poreVoisinsBin,poreVoisinsBinCenters] =hist(poreVoisins,max(poreVoisins));
    poreVoisinsData=[poreVoisinsBin,poreVoisinsBinCenters]
    
    return nPore,nLink,poreDiamData,linkDiamData,poreVoisinsData
    
    
    
legende=['All Local Maxima','H-maxima with H-constrast=4 voxels',
         'H-maxima with H-constrast=8 voxels','H-maxima with H-constrast=12 voxels',
         'H-maxima with H-constrast=20 voxels','Local Maxima 4 voxels apart',
         'Local Maxima 10 voxels apart','Local Maxima 20 voxels apart']   


nPore_ = []
nLink_ = []
poreDiamBinsOccurrences_ = []
poreDiamBinCenters_ = []
linkDiamBins_ = []
linkDiamBinCenters_ = []
poreVoisinsBin_ = []
poreVoisinsBinCenters_ = []

for i in range(nNetwork):
    nPore,nLink,poreDiamData,linkDiamData,poreVoisinsData=ComputeCaracteristics(
                                        poreDiams[i],linkRadius[i],coordination[i])
    
    nPore_.append(                   nPore)
    nLink_.append(                   nLink)
    poreDiamBinsOccurrences_.append( poreDiamData[0])
    poreDiamBinsVolume_.append(      poreDiamData[1])
    poreDiamBinCenters_.append(      poreDiamData[2])
    linkDiamBins_.append(            linkDiamData[0])
    linkDiamBinCenters_.append(      linkDiamData[1])
    poreVoisinsBin_.append(          poreVoisinsData[0])
    poreVoisinsBinCenters_.append(   poreVoisinsData[1])


#%Plots
#
#%     colors={'green',...
#%         'blue',...
#%         'cyan',...
#%         'red',...
#%         'magenta',...
#%         'black',...
#%         'yellow' };


fig = _plt.figure()
dp = network['pore.diameter']
Vp = network['pore.volume']
dt = network['throat.diameter']
Vt = network['throat.volume']
dmax = max(max(dp), max(dt))
steps = _sp.linspace(0, dmax, 100, endpoint=True)
vals = _sp.zeros_like(steps)
for i in range(0, len(steps)-1):
    temp1 = dp > steps[i]
    temp2 = dt > steps[i]
    vals[i] = sum(Vp[temp1]) + sum(Vt[temp2])
yaxis = vals
xaxis = steps
_plt.semilogx(xaxis, yaxis, 'b.-')
_plt.xlabel('Pore & Throat Diameter (m)')
_plt.ylabel('Cumulative Volume (m^3)')
return fig


fig = _plt.figure()
for i in range(nNetwork):
#    %plot(poreDiamBinCenters_{i},poreDiamBins_{i}/nPore_{i},'Color',colors{i},'LineWidth',6);
#    %plot(poreDiamBinCenters_{i},poreDiamBinsVolume_{i},'Color',colors{i},'LineWidth',6);
    yaxis = poreDiamBinsVolume_[i]
    xaxis = poreDiamBinCenters_[i]
    _plt.semilogx(xaxis, yaxis, 'b.-')
    _plt.xlabel('Pore Diameter (m)')
    _plt.ylabel('Cumulative Volume (m^3)')
    
#    volumeContribution=poreDiamBinsOccurrences_{i}.*(4/3)*pi.*(poreDiamBinCenters_{i}/2).^3;
#    plot(poreDiamBinCenters_{i},volumeContribution,'Color',colors{i},'LineWidth',6);



title('Distribution of pore diameters')
axe = get(fig, 'Children');
legend(axe,legende)


fig=figure;hold on;
for i in range(nNetwork):
#    %plot(linkDiamBinCenters_{i},linkDiamBins_{i}/nLink_{i},'Color',colors{i},'LineWidth',6);
    plot(linkDiamBinCenters_{i},linkDiamBins_{i},'Color',colors{i},'LineWidth',6);

title('Distribution of link diameters')
axe = get(fig, 'Children');
legend(axe,legende)


fig=figure;hold on;
for i in range(nNetwork):
#    %plot(poreVoisinsBinCenters_{i},poreVoisinsBin_{i}/nPore_{i},'Color',colors{i},'LineWidth',6)
    plot(poreVoisinsBinCenters_{i},poreVoisinsBin_{i},'Color',colors{i},'LineWidth',6)

title('Distribution of coordination number')
axe = get(fig, 'Children');
legend(axe,legende)