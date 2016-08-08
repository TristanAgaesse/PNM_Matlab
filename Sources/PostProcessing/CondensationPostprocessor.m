classdef CondensationPostprocessor
    %CondensationPostprocessor Tools to postprocess condensation results 
    % for further custom analysis
    
    properties
        CondensationInfos
        CondensationClusters
        CondensationOptions
        Network
    end
    
    methods
        
        function postProcessor = CondensationPostprocessor(network,condensationOptions,condensationInfos,condensationClusters)
            % Constructor function
            % Example : [condensationInfos,condensationClusters]=ComputeCondensation(network, options)
            %           postProcessor = PostprocessCondensation(condensationInfos,condensationClusters,condensationOptions,network)
            
            postProcessor.Network = network;
            postProcessor.CondensationOptions = condensationOptions;
            postProcessor.CondensationInfos = condensationInfos;
            postProcessor.CondensationClusters = condensationClusters;
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
            %   [geometryRegions,geometryRegions3DIndex] = postProcessor.BuildGeometryRegions_XYZDiscretized(infos_Discretization)
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
                zfilter = and(poreCenters(:,3)>=thisRegion.Ymin,poreCenters(:,3)<=thisRegion.Ymax) ;
                
            end
            
        end
         
        function blockAverage = GetFieldBlockAveraged(postProcessor,field,poreLookUpTable)
            % Input : postProcessor,field,poreLookUpTable
            % Output : blockAverage (mean value of field for each block)
            % Example : field = rand(nPore,1)
            %           postProcessor.GetFieldBlockAveraged(field,poreLookUpTable) 
            
            nBlock = max(poreLookUpTable);
            blockAverage = zeros(1,nBlock);
            poreVolume=postProcessor.Network.GetPoreData('Volume');
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
            %   [geometryRegions,geometryRegions3DIndex] = postProcessor.BuildGeometryRegions_XYZDiscretized(infos_Discretization)
            %   region_x1y2z5 = geometryRegions{geometryRegions3DIndex(1,2,5)}
            
            
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
        

        
        
        %% Cluster evolution tree
        function [tree,treeInfos] = BuildClusterEvolutionTree(postProcessor)
            %returns a tree of the clusters evolution during CondendationClusterGrowth
            %Input : - postProcessor
            %Output : tree,treeInfos 
            %Example : [tree,treeInfos] = postProcessor.BuildClusterEvolutionTree
            
            error('Not Implemented')
            
            % detect cluster coalescence            
            
            %Correler tree a taille d amas, positions, taux de croissance ?
            
        end
        
        
        %% Getters on field evolution
        
        function field = GetFieldVaporConcentration(postProcessor,iStep)
            
            
            nPore = postProcessor.Network.GetNumberOfPores;
            field = zeros(nPore,1);
            
            %CondensationInfos.     ; separer nucleation et growth ?
            
            %put c=c_eq in liquid invaded pores 
            
        end
        
        function field = GetFieldLiquidSaturation(postProcessor,iStep)
            % saturation : 0 ou 1 pour chaque pore
            
            nPore = postProcessor.Network.GetNumberOfPores;
            field = zeros(nPore,1);
            
            %CondensationInfos.     ; separer nucleation et growth ?
            
        end
        
        %% Plots
        
        % Global liquid saturation (as a function of anisotropy degree or GDL penetration) 
        % In-plane and throuugh plane saturation profiles (comme dans la thèse de benjamin)
        % 3D patterns (cf thèse Benjamin)
        
    end
    
    
end

