classdef DataVerticeList <  Data
    %DATAVERTICELIST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        VerticeDatas
        NombreVertice %nombre d'ar�tes
    end
    
    methods
        function data_Vertice_list=DataVerticeList(nVertice)
            data_Vertice_list.NombreVertice=nVertice;
            data_Vertice_list.VerticeDatas=struct;
        end
        
        function AddData(data_Vertice_list,data,name)
            assert(length(data)==data_Vertice_list.NombreVertice,'Un vertice data doit �tre un tableau de taille NombreVertice');
            data_Vertice_list.VerticeDatas.(name)=data;
        end
        
        function RemoveData(data_list,name)
            data_list.VerticeDatas=rmfield(data_list.VerticeDatas,name);
        end           
        
    end
    
end

