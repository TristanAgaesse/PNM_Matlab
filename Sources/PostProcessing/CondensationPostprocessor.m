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
        
        
        
        %% Getters on field evolution
        
        function field = GetFieldVaporConcentration(postProcessor,iCondensationStep)
            
            
            nPore = postProcessor.Network.GetNumberOfPores;
            field = zeros(nPore,1);
            
            %CondensationInfos.     ; separer nucleation et growth ?
            
            %put c=c_eq in liquid invaded pores 
            
        end
        
        function field = GetFieldLiquidSaturation(postProcessor,iCondensationStep)
            % saturation : 0 ou 1 pour chaque pore
            
            nPore = postProcessor.Network.GetNumberOfPores;
            field = zeros(nPore,1);
            
            %CondensationInfos.     ; separer nucleation et growth ?
            
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
        
        
        %% Plots
        
        % Global liquid saturation (as a function of anisotropy degree or GDL penetration) 
        % In-plane and through plane saturation profiles (comme dans la thèse de benjamin)
        % 3D patterns (cf thèse Benjamin)
        
        function PlotSaturationProfile(postProcessor,iCondensationStep,direction,nDiscretization)
            
            saturation=postProcessor.GetFieldLiquidSaturation(iCondensationStep);
            
            volumePostProcessor = PostprocessorVolumeAverage(postProcessor.Network);
            volumePostProcessor.PlotProfile(saturation,direction, nDiscretization) 
        end
        
        
    end
    
    
end

