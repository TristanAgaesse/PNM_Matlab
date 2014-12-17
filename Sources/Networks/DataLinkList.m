classdef DataLinkList <  Data
    %DATALINKLIST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LinkDatas
        NombreLink %nombre de liens
    end
    
    methods
        function data_link_list=DataLinkList(nLink)
            data_link_list.NombreLink=nLink;
            data_link_list.LinkDatas=struct;
        end
        
        function AddData(data_link_list,data,name)
            assert(length(data)==data_link_list.NombreLink,'Un link data doit �tre un tableau de taille NombreLinks');
            data_link_list.LinkDatas.(name)=data;
        end
       
        function RemoveData(data_link_list,name)
            data_link_list.LinkDatas=rmfield(data_link_list.LinkDatas,name);
        end
        
    end
    
end

