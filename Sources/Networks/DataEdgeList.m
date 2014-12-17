classdef DataEdgeList <  Data
    %DATAEDGELIST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        EdgeDatas
        NombreEdges %nombre d'ar�tes
    end
    
    methods
        function data_edge_list=DataEdgeList(nEdges)
            data_edge_list.NombreEdges=nEdges;
            data_edge_list.EdgeDatas=struct;
        end
        
        function AddData(data_edge_list,data,name)
            assert(length(data)==data_edge_list.NombreEdges,'Un edge data doit �tre un tableau de taille NombreEdges');
            data_edge_list.EdgeDatas.(name)=data;
        end
        
        function RemoveData(data_list,name)
            data_list.EdgeDatas=rmfield(data_list.EdgeDatas,name);
        end        
        
    end
    
end

