

figName={'../ResultsFullMorphology/FullMorpho_theta115_VoxelMatchTotal_slices15-85.fig',...
         '../ResultsPNM/LocalMax1_theta115_voxelMatchTotal.fig', ...
         '../ResultsPNM/LocalMax10_theta115_voxelMatchTotal.fig', ...
         '../ResultsPNM/h8_theta115_voxelsMatchTotal.fig'  
         '../ResultsFullMorphology/FullMorpho_theta115_capPress_slices15-85_simuOnly.fig',...
         '../ResultsPNM/LocalMax1_theta115_capPressure_simuOnly.fig', ...
         '../ResultsPNM/LocalMax10_theta115_capPress_simuOnly.fig', ...
         '../ResultsPNM/h8_theta115_capPress_simuOnly.fig' 
         };
         

fig1 = open(figName{1});
fig2 = open(figName{2});
fig3 = open(figName{3});
fig4 = open(figName{4});
fig5 = open(figName{5});
fig6 = open(figName{6});
fig7 = open(figName{7});
fig8 = open(figName{8});

axe = get(fig1, 'Children');

axfig2 = get(fig2, 'Children');
for i = 1 : numel(axfig2)
    ax2Children = get(axfig2(i),'Children');copyobj(ax2Children, axe(i));
end

axfig3 = get(fig3, 'Children');
for i = 1 : numel(axfig3)
    ax2Children = get(axfig3(i),'Children');copyobj(ax2Children, axe(i));
end

axfig4 = get(fig4, 'Children');
for i = 1 : numel(axfig4)
    ax2Children = get(axfig4(i),'Children');copyobj(ax2Children, axe(i));
end

axfig5 = get(fig5, 'Children');
for i = 1 : numel(axfig5)
    ax2Children = get(axfig5(i),'Children');copyobj(ax2Children, axe(i));
end

axfig6 = get(fig6, 'Children');
for i = 1 : numel(axfig6)
    ax2Children = get(axfig6(i),'Children');copyobj(ax2Children, axe(i));
end

axfig7 = get(fig7, 'Children');
for i = 1 : numel(axfig7)
    ax2Children = get(axfig7(i),'Children');copyobj(ax2Children, axe(i));
end

axfig8 = get(fig8, 'Children');
for i = 1 : numel(axfig8)
    ax2Children = get(axfig8(i),'Children');copyobj(ax2Children, axe(i));
end

leg={'VoxelMatch Full Morphology', ...
    'VoxelMatch Pore Network - Markers = all local maxima', ...
    'VoxelMatch Pore Network - Markers = local maxima 10 voxels apart', ...
    'VoxelMatch Pore Network - Markers = H-maxima (8 voxels contrast)'
    'Saturation Full Morphology', ...
    'Saturation Pore Network - Markers = all local maxima', ...
    'Saturation Pore Network - Markers = local maxima 10 voxels apart', ...
    'Saturation Pore Network - Markers = H-maxima (8 voxels contrast)'
    };


legend(axe,leg)