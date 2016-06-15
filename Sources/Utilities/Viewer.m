classdef Viewer
%VIEWER Visualise un réseau avec Matlab
%   Main function : Viewer.View
    
    properties
        NetworkGeometry
    end
    
    methods
        function viewer=Viewer(output_struct)
            viewer.NetworkGeometry=output_struct;
        end
        
        function View(viewer,options,varargin)
        %Input : viewer,options, (+ varargin : facultative input)
        % 'Network': affichage du maillage
        % 'PoreList': affiche une liste de pores : liste_cellules=varargin{1}
        % 'PoreField' : affichage d'un champ scalaire sur les cellules, champ=varargin{1}
        % 'LinkField' : affichage d'un champ scalaire sur les liens, champ=varargin{1}
        % 'Boundaries': affichage des faces du maillage, avec code couleur pour les fronti�res
        % 'Edges' : affichage des aretes du maillage, avec code couleur pour les frontieres
        % 'Sommets' : affichage des fibres avec des disques
            
            dimension=viewer.NetworkGeometry.ATTRIBUTE.Dimension;
            vertices=viewer.NetworkGeometry.Vertices;
            cells_to_vertices=viewer.NetworkGeometry.CellsToVertices;
            faces=viewer.NetworkGeometry.Faces.Face;
            boundaries=viewer.NetworkGeometry.Boundaries.Boundary;
            edges=viewer.NetworkGeometry.Edges;
            nFace=length(faces);
            nCell=length(cells_to_vertices);
            nEdge=length(edges(:,1));
            
            %Affichage du maillage
            if strcmp('Network',options)
                if dimension==2
                    colors=cell(1,nCell);
                    for iCell=1:nCell
                        if isempty(varargin)
                            colors{iCell}=iCell;
                        else
                            colors{iCell}=varargin{iCell};
                        end
                    end
                    for iCell = 1:nCell
                        patch(vertices(cells_to_vertices{iCell},1),vertices(cells_to_vertices{iCell},2),colors{iCell});
                    end
                elseif dimension==3
                    colors=cell(1,nFace);
                    for iLink=1:nFace
                        if isempty(varargin)
                            colors{iLink}=iLink;
                        else
                            colors{iLink}=varargin{iLink};
                        end
                    end
                    for iLink=1:nFace
                        patch(vertices(faces{iLink},1),vertices(faces{iLink},2),vertices(faces{iLink},3),colors{iLink});            
                    end
                end
            end
            
            %Affichage liste de cellules
            if strcmp('PoreList',options)
                liste_cellules=varargin{1};
                for i = 1:length(liste_cellules)
                    foo=i/length(liste_cellules);
                    
                    if length(varargin)>1
                        if strcmp('blue',varargin{2})
                            color='blue';
                        elseif strcmp('gris',varargin{2})
                            color=[foo,foo,foo];
                        elseif strcmp('white',varargin{2})
                            color=[1 1 1];
                        end
                    else
                        color='blue';
                    end
                    
                    if dimension==2
                        patch(vertices(cells_to_vertices{liste_cellules(i)},1),vertices(cells_to_vertices{liste_cellules(i)},2),color);
                    elseif dimension==3
                        patch(vertices(cells_to_vertices{liste_cellules(i)},1),vertices(cells_to_vertices{liste_cellules(i)},2),vertices(cells_to_vertices{liste_cellules(i)},3),color);            
                    end
                end
            end
            
            %Affichage d'un champ scalaire sur les pores
            if strcmp('PoreField',options)
                champ_cellules=varargin{1};
                assert(length(champ_cellules)==length(cells_to_vertices),'L''input doit donner un scalaire par cell');
                
                colormap(jet);
                colors=viewer.GetCustomColorField(champ_cellules);
                
                for i = 1:length(champ_cellules)
                    if dimension==2
                        patch(vertices(cells_to_vertices{i},1),vertices(cells_to_vertices{i},2),colors{i});
                    elseif dimension==3
                        patch(vertices(cells_to_vertices{i},1),vertices(cells_to_vertices{i},2),vertices(cells_to_vertices{i},3),colors{i});            
                    end
                    hold on
                end
            end
            
            %Affichage d'un champ scalaire sur les liens
            if strcmp('LinkField',options)
                champ_cellules=varargin{1};
                assert(length(champ_cellules)==nEdge,'varargin{1} doit être une liste de taille nLink')
                
                inf=min(champ_cellules);
                sup=max(champ_cellules);
                if sup==inf
                   delta=1;
                else
                    delta=sup-inf;
                end
                
                colormap(jet(128));
                map=colormap;
                color_foo=zeros(1,length(champ_cellules));
                for i = 1:length(champ_cellules)
                    foo=(champ_cellules(i)-inf)/(delta);
                    color_foo(i)=foo;
                end
                colors=cell(1,length(champ_cellules));
                for i = 1:length(champ_cellules)
%                     if ceil(128*color_foo(i))>0
%                         colors{i}=map(ceil(128*color_foo(i)),:);
%                     else
%                         colors{i}=[1,1,1];
%                     end
                    colors{i}=map(floor(127*color_foo(i))+1,:);
                end
                
                for iEdge = 1:nEdge
                    if dimension==2
                        p=patch(vertices(edges(iEdge,:),1),vertices(edges(iEdge,:),2),floor(128*color_foo(i))+1);
                        set(p,'EdgeColor',colors{iEdge},'LineWidth',floor(4*color_foo(i))+1);
                    elseif dimension==3
                        p=patch(vertices(edges(iEdge,:),1),vertices(edges(iEdge,:),2),vertices(edges(iEdge,:),3),colors{iEdge});
                    end
                end
            end
            
            
            
            
            %Affichage des faces du maillage, avec code
            %couleur pour les fronti�res
            if strcmp('Boundaries',options)
                nBoundary=length(boundaries);
                colormap(hsv(nBoundary));
                colors=colormap;
                
                for iBoundary=1:nBoundary
                    startFace=boundaries(iBoundary).ATTRIBUTE.StartFace;
                    endFace=boundaries(iBoundary).ATTRIBUTE.StartFace+boundaries(iBoundary).ATTRIBUTE.NombreFaces-1;
                    for iLink=startFace:endFace
                        if dimension==2
                            p=patch(vertices(faces{iLink},1),vertices(faces{iLink},2),iBoundary);
                            set(p,'EdgeColor',colors(iBoundary,:),'LineWidth',3);
                        elseif dimension==3
                            p=patch(vertices(faces{iLink},1),vertices(faces{iLink},2),vertices(faces{iLink},3),iBoundary);
                            set(p,'FaceColor',colors(iBoundary,:))
                        end
                    end
                end
            end

            %Affichage des ar�tes du maillage
            if strcmp('Edges',options)
                epaisseurs=varargin{1};
                assert(length(epaisseurs)==nEdge,'varargin{1} doit être une liste d''epaisseurs de taille nEdge')
                black=[0 0 0];
                for iEdge = 1:nEdge
                    if dimension==2
                        p=patch(vertices(edges(iEdge,:),1),vertices(edges(iEdge,:),2),black);
                        set(p,'LineWidth',epaisseurs(iEdge));
                    elseif dimension==3
                        p=patch(vertices(edges(iEdge,:),1),vertices(edges(iEdge,:),2),vertices(edges(iEdge,:),3),black);
                    end
                end
            end
            
            %Affichage des fibres
            if strcmp('Sommets',options)
                taille_fibres=varargin{1};
                if dimension==2
                    assert(length(taille_fibres)==length(vertices(:,1)),'Input : tableau taille de chaque fibre');
                    for num_vertex=1:length(taille_fibres)
                        center=vertices(num_vertex,:);
                        radius=taille_fibres(num_vertex);
                        
                        r=radius;
                        N=8;
                        color=[0,0,0];
                        THETA=linspace(0,2*pi,N);
                        RHO=ones(1,N)*r;
                        [X,Y] = pol2cart(THETA,RHO);
                        X=X+center(1);
                        Y=Y+center(2);
                        patch(X,Y,color);
                    end
                end
            end



            %Affichage pas � pas d'une liste ordonn�e de cellules
%             if strcmp('cellules_pas_a_pas',options)
%                 liste_cellules=varargin{1};
%                 num_pas=varargin{2};
%                 assert(num_pas<=length(liste_cellules),'numero maximum d�pass�');
% 
%                 bleu=[0,0,1];
%                 if dimension==2
%                     patch(vertices(cells_to_vertices{liste_cellules(num_pas)},1),vertices(cells_to_vertices{liste_cellules(num_pas)},2),bleu);
%                 elseif dimension==3
%                     patch(vertices(cells_to_vertices{liste_cellules(num_pas)},1),vertices(cells_to_vertices{liste_cellules(num_pas)},2),vertices(cells_to_vertices{liste_cellules(num_pas)},3),bleu);            
%                 end
% 
%             end

        end
        
        
        
        
        
        
        
        
        
        
        function ViewNetwork(viewer,options,varargin)
            dimension=viewer.NetworkGeometry.ATTRIBUTE.Dimension;
            vertices=viewer.NetworkGeometry.Vertices;
            cells_to_vertices=viewer.NetworkGeometry.CellsToVertices;
            faces=viewer.NetworkGeometry.Faces.Face;
            boundaries=viewer.NetworkGeometry.Boundaries.Boundary;
            edges=viewer.NetworkGeometry.Edges;
            
            nFace=length(faces);
            nCell=length(cells_to_vertices);
            nEdge=length(edges(:,1));
            
            if dimension==2
                colors=cell(1,nCell);
                for iCell=1:nCell
                    if isempty(varargin)
                        colors{iCell}=iCell;
                    else
                        colors{iCell}=varargin{iCell};
                    end
                end
                for iCell = 1:nCell
                    patch(vertices(cells_to_vertices{iCell},1),vertices(cells_to_vertices{iCell},2),colors{iCell});
                end
            elseif dimension==3
                colors=cell(1,nFace);
                for iLink=1:nFace
                    if isempty(varargin)
                        colors{iLink}=iLink;
                    else
                        colors{iLink}=varargin{iLink};
                    end
                end
                for iLink=1:nFace
                    patch(vertices(faces{iLink},1),vertices(faces{iLink},2),vertices(faces{iLink},3),colors{iLink});            
                end
            end
        end
        
        function ViewPoreList(viewer,options,varargin)
            
            dimension=viewer.NetworkGeometry.ATTRIBUTE.Dimension;
            vertices=viewer.NetworkGeometry.Vertices;
            cells_to_vertices=viewer.NetworkGeometry.CellsToVertices;
            faces=viewer.NetworkGeometry.Faces.Face;
            boundaries=viewer.NetworkGeometry.Boundaries.Boundary;
            edges=viewer.NetworkGeometry.Edges;
            nFace=length(faces);
            nCell=length(cells_to_vertices);
            nEdge=length(edges(:,1));
            
            liste_cellules=varargin{1};
            
            for i = 1:length(liste_cellules)
                foo=i/length(liste_cellules);

                if length(varargin)>1
                    if strcmp('blue',varargin{2})
                        color='blue';
                    elseif strcmp('gris',varargin{2})
                        color=[foo,foo,foo];
                    elseif strcmp('white',varargin{2})
                        color=[1 1 1];
                    end
                else
                    color='blue';
                end

                if dimension==2
                    patch(vertices(cells_to_vertices{liste_cellules(i)},1),vertices(cells_to_vertices{liste_cellules(i)},2),color);
                elseif dimension==3
                    patch(vertices(cells_to_vertices{liste_cellules(i)},1),vertices(cells_to_vertices{liste_cellules(i)},2),vertices(cells_to_vertices{liste_cellules(i)},3),color);            
                end
            end
        
        end
        
        function ViewPoreData(viewer,poreData,colorOption)
            
            dimension=viewer.NetworkGeometry.ATTRIBUTE.Dimension;
            vertices=viewer.NetworkGeometry.Vertices;
            cells_to_vertices=viewer.NetworkGeometry.CellsToVertices;
            nPore = length(cells_to_vertices);
            assert(length(poreData)==nPore,'L''input doit donner un scalaire par cell');

            colormap(jet);
            if length(colorOption)==2
                caxis(colorOption)
                colors=cell(1,nPore);
                for i=1:nPore
                    colors{i}=poreData(i);
                end
            elseif strcmp(colorOption, 'MapTo01')
                colors=GetCustomColorField(poreData);
                caxis auto
            end    

            for i = 1:length(poreData)
                if dimension==2
                    patch(vertices(cells_to_vertices{i},1),vertices(cells_to_vertices{i},2),colors{i});
                elseif dimension==3
                    patch(vertices(cells_to_vertices{i},1),vertices(cells_to_vertices{i},2),vertices(cells_to_vertices{i},3),colors{i});            
                end
                hold on
            end
        end
        
        function ViewLinkField(viewer,options,varargin)
            
            dimension=viewer.NetworkGeometry.ATTRIBUTE.Dimension;
            vertices=viewer.NetworkGeometry.Vertices;
            cells_to_vertices=viewer.NetworkGeometry.CellsToVertices;
            faces=viewer.NetworkGeometry.Faces.Face;
            boundaries=viewer.NetworkGeometry.Boundaries.Boundary;
            edges=viewer.NetworkGeometry.Edges;
            nFace=length(faces);
            nCell=length(cells_to_vertices);
            nEdge=length(edges(:,1));
            
            champ_cellules=varargin{1};
            assert(length(champ_cellules)==nEdge,'varargin{1} doit être une liste de taille nLink')

            inf=min(champ_cellules);
            sup=max(champ_cellules);
            if sup==inf
               delta=1;
            else
                delta=sup-inf;
            end

            colormap(jet(128));
            map=colormap;
            color_foo=zeros(1,length(champ_cellules));
            for i = 1:length(champ_cellules)
                foo=(champ_cellules(i)-inf)/(delta);
                color_foo(i)=foo;
            end
            colors=cell(1,length(champ_cellules));
            for i = 1:length(champ_cellules)
%                     if ceil(128*color_foo(i))>0
%                         colors{i}=map(ceil(128*color_foo(i)),:);
%                     else
%                         colors{i}=[1,1,1];
%                     end
                colors{i}=map(floor(127*color_foo(i))+1,:);
            end

            for iEdge = 1:nEdge
                if dimension==2
                    p=patch(vertices(edges(iEdge,:),1),vertices(edges(iEdge,:),2),floor(128*color_foo(i))+1);
                    set(p,'EdgeColor',colors{iEdge},'LineWidth',floor(4*color_foo(i))+1);
                elseif dimension==3
                    p=patch(vertices(edges(iEdge,:),1),vertices(edges(iEdge,:),2),vertices(edges(iEdge,:),3),colors{iEdge});
                end
            end
        end
        
        
       
        
        function colors=GetCustomColorField(viewer,champ_cellules)
            
            inf=min(champ_cellules(champ_cellules~=0));
            sup=max(champ_cellules(champ_cellules~=0));
            delta=sup-inf;
            if abs(delta/mean(champ_cellules(champ_cellules~=0)))<1e-5;
                delta=1;
            end
            color_foo=zeros(1,length(champ_cellules));
            for i = 1:length(champ_cellules)
                if champ_cellules(i)~=0
                    %foo=1-(champ_cellules(i)-inf)/(sup-inf);
                    foo=(champ_cellules(i)-inf)/delta;
                    color_foo(i)=foo;
                else
                    color_foo(i)=1;
                end
            end
%               color_foo=imadjust ([color_foo',color_foo',color_foo']); %augmentation du contraste
%               color_foo=color_foo(:,1);
            colors=cell(1,length(champ_cellules));
            for i = 1:length(champ_cellules)
                if champ_cellules(i)~=0
                    colors{i}=color_foo(i);
                    %colors{i}=64*color_foo(i);
                else
                    colors{i}=[1,1,1];
                end
            end
        end
            
    end
    
end

