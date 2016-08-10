classdef PostProcessorVolumeAverage
    %PostprocessorVolumeAverage Tools to make averages of network fields on
    % custom spatial discretizations
    
    properties
        Network
    end
    
    methods
        
        function postProcessor = PostProcessorVolumeAverage(network)
            % Constructor function
            % Example : postProcessor = PostProcessorVolumeAverage(network)
            
            postProcessor.Network = network;
        end
        
        %% Upscaled regions look up table
        function poreLookUpTable = BuildPoreLookUpTable(postProcessor,geometryRegions)
            %returns a look up table to know which pores are in each block, and which blocks contains a pore
            %Input : - postProcessor
            %        - geometryRegions : cell containing infos on geometry
            %                           discretisation
            %Output : poreLookUpTable : array(1,nPore), poreLookUpTable(iPore)=iBlock
            % Example : 
            %   infos_Discretization=struct
            %   infos_Discretization.X = [xmin,xmax];
            %   infos_Discretization.Y = [ymin,ymean,ymax];
            %   infos_Discretization.Z = linspace(zmin,zmax,10);
            %   [geometryRegions,geometryRegions3DIndex] = postProcessor.BuildGeometryRegions3DDiscretized(infos_Discretization)
            %   region_x1y2z5 = geometryRegions{geometryRegions3DIndex(1,2,5)}
            %   poreLookUpTable = postProcessor.BuildPoreLookUpTable(geometryRegions)
            %   poresInRegion_x1y2z5 = find(ismember(poreLookUpTable,geometryRegions3DIndex(1,2,5)))
            
            nPore = postProcessor.Network.GetNumberOfPores;
            poreLookUpTable = zeros(1,nPore);
            
            nRegion=length(geometryRegions);
            poreCenters = postProcessor.Network.GetPoreCenter(1:nPore);
            
            for iRegion=1:nRegion
                thisRegion = geometryRegions{iRegion};
                
                xfilter = and(poreCenters(:,1)>=thisRegion.Xmin,poreCenters(:,1)<=thisRegion.Xmax) ;
                yfilter = and(poreCenters(:,2)>=thisRegion.Ymin,poreCenters(:,2)<=thisRegion.Ymax) ;
                zfilter = and(poreCenters(:,3)>=thisRegion.Zmin,poreCenters(:,3)<=thisRegion.Zmax) ;
                
                poreLookUpTable(and(and(xfilter,yfilter),zfilter))=iRegion;
            end
            
            assert(all(poreLookUpTable>0))
        end
        
        
        function blockAverage = GetFieldBlockAveraged(postProcessor,field,poreLookUpTable)
            % Input : postProcessor,field,poreLookUpTable
            % Output : blockAverage (mean value of field for each block)
            % Example : field = rand(nPore,1)
            %           postProcessor.GetFieldBlockAveraged(field,poreLookUpTable) 
            
            nBlock = max(poreLookUpTable);
            blockAverage = zeros(1,nBlock);
            poreVolume=postProcessor.Network.GetPoreData('Volume');
            field=field(:);
            for iBlock=1:nBlock
                poreFilter = (poreLookUpTable == iBlock);
                blockVolume = sum(poreVolume(poreFilter));
                fieldIntegral = sum(field(poreFilter).*poreVolume(poreFilter));
                blockAverage(iBlock) = fieldIntegral/blockVolume;
            end
        end
        
        
        
        function [geometryRegions,geometryRegions3DIndex] = BuildGeometryRegions3DDiscretized(postProcessor,infos_Discretization)
            % Put
            % Input : postProcessor,infos_Discretization
            % Output : geometryRegions,geometryRegions3DIndex
            % Example : 
            %   infos_Discretization=struct
            %   infos_Discretization.X = [xmin,xmax];
            %   infos_Discretization.Y = [ymin,ymean,ymax];
            %   infos_Discretization.Z = linspace(zmin,zmax,10);
            %   [geometryRegions,geometryRegions3DIndex] = postProcessor.BuildGeometryRegions3DDiscretized(infos_Discretization)
            %   region_x1y2z5 = geometryRegions{geometryRegions3DIndex(1,2,5)}
            
            
            nx = length(infos_Discretization.X)-1;
            ny = length(infos_Discretization.Y)-1;
            nz = length(infos_Discretization.Z)-1;
            
            geometryRegions3DIndex = zeros([nx,ny,nz]);
            nRegion=nx*ny*nz;
            geometryRegions=cell(1,nRegion);
            
            for iRegion=1:nRegion
                [i,j,k]=ind2sub([nx,ny,nz],iRegion);
                geometryRegions3DIndex(i,j,k)=iRegion;
                
                thisRegion=struct;
                thisRegion.Xmin = infos_Discretization.X(i);
                thisRegion.Xmax = infos_Discretization.X(i+1);
                thisRegion.Ymin = infos_Discretization.Y(j);
                thisRegion.Ymax = infos_Discretization.Y(j+1);
                thisRegion.Zmin = infos_Discretization.Z(k);
                thisRegion.Zmax = infos_Discretization.Z(k+1);
                
                geometryRegions{iRegion}=thisRegion;
            end
        end
        
        
        function [geometryRegions,geometryRegions2DIndex] = BuildGeometryRegions2DDiscretized(postProcessor,infos_Discretization)
            % Put
            % Input : postProcessor,infos_Discretization
            % Output : geometryRegions,geometryRegions3DIndex
            % Example : 
            %   infos_Discretization=struct
            %   infos_Discretization.X = [xmin,xmax];
            %   infos_Discretization.Y = [ymin,ymean,ymax];
            %   [geometryRegions,geometryRegions2DIndex] = postProcessor.BuildGeometryRegions2DDiscretized(infos_Discretization)
            %   region_x1y2 = geometryRegions{geometryRegions2DIndex(1,2)}
            
            
            infos_Discretization.Z=[-Inf,+Inf];
            [geometryRegions,geometryRegions3DIndex] = GetGeometryRegions3DDiscretized(postProcessor,infos_Discretization);
            
            geometryRegions2DIndex = geometryRegions3DIndex(:,:,1);
        end
        

        %% Plots
        
        function PlotProfile(postProcessor,field,direction, nStep)
            % input : - postProcessor
            %         - field : pore data which profile you want to plot 
            %         - direction : 'x','y' ou 'z' (direction of discretization)
            %         - nStep : discretization
            %
            %Example : postProcessor.PlotProfile(saturation,'z', 10) %plot saturation profile
            
            
            %Check inputs
            isX = strcmp(direction,'x');
            isY = strcmp(direction,'y');
            isZ = strcmp(direction,'z');
            assert( isX || isY || isZ,'direction must equal "x","y" or "z"')
            
            nx=1;
            ny=1;
            nz=1;
            assert(nStep>0 && ceil(nStep)==nStep,'nStep must be a positive interger')
            if isX
                nx=nStep;
            elseif isY
                ny=nStep;
            else
                nz=nStep;
            end
            
            %Build discretization look up table
            nPore = postProcessor.Network.GetNumberOfPores;
            dimension = postProcessor.Network.GetDimension;
            poreCenters = postProcessor.Network.GetPoreCenter(1:nPore);
            
            xmin = min(poreCenters(:,1));
            xmax = max(poreCenters(:,1));
            ymin = min(poreCenters(:,2));
            ymax = max(poreCenters(:,2));
            if dimension==2
                zmin = -Inf;
                zmax = Inf;
            elseif dimension==3
                zmin = min(poreCenters(:,3));
                zmax = max(poreCenters(:,3));
            end
            
            infos_Discretization=struct;
            infos_Discretization.X = linspace(xmin,xmax,nx+1);
            infos_Discretization.Y = linspace(ymin,ymax,ny+1);
            infos_Discretization.Z = linspace(zmin,zmax,nz+1);
            [geometryRegions,geometryRegions3DIndex] = postProcessor.BuildGeometryRegions3DDiscretized(infos_Discretization);
            poreLookUpTable = postProcessor.BuildPoreLookUpTable(geometryRegions);
            
            
            %Get data along the profile
            blockAverage = postProcessor.GetFieldBlockAveraged(field,poreLookUpTable);
            profile=zeros(1,nStep);
            for iStep=1:nStep
                [i,j,k]=ind2sub([nx,ny,nz],iStep);
                profile(iStep)=blockAverage(geometryRegions3DIndex(i,j,k));
            end
            
            
            % Plot data
            if isX
                abscisse = (infos_Discretization.X(1:end-1)+infos_Discretization.X(2:end))/2;
            elseif isY
                abscisse = (infos_Discretization.Y(1:end-1)+infos_Discretization.Y(2:end))/2;
            else
                abscisse = (infos_Discretization.Z(1:end-1)+infos_Discretization.Z(2:end))/2;
            end
            vname=@(x) inputname(1);
            fieldname = vname(field);
            
            plot(abscisse,profile)
            xlabel(direction)
            ylabel(fieldname)
            mytitle=strcat('Profile of ',fieldname);
            mytitle=strcat(mytitle,' along direction ');
            mytitle=strcat(mytitle,direction);
            title(mytitle) 
        end
        
        
    end
end

