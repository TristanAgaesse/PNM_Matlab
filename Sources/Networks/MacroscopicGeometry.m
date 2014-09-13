classdef MacroscopicGeometry < handle
    %MACROSCOPICGEOMETRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Dimension
        Vertices
        Blocks
        Boundaries
        AnisotropyBand
        NVertice
        NBlock
        NBoundary
        NAnisotropyBand
        NetworkType
        PavageType
        ConvertToMeters
        RandomSeed
    end
    
    methods
        
        function myGeometry=MacroscopicGeometry()
            myGeometry.Dimension=[];
            myGeometry.NVertice=0;
            myGeometry.NBlock=0;
            myGeometry.NBoundary=0;
            myGeometry.NAnisotropyBand=0;
            myGeometry.Vertices=[];
            myGeometry.Blocks=cell(1,0);
            myGeometry.Boundaries=cell(1,0);
            myGeometry.AnisotropyBand=cell(1,0);
            myGeometry.NetworkType='';
            myGeometry.PavageType='';
            myGeometry.ConvertToMeters=1;
            myGeometry.RandomSeed=0;
        end
        
        
        %Functions to access geometry
        
        function dimension=GetDimension(myGeometry)
            dimension=myGeometry.Dimension;
        end
        
        function networkType=GetNetworkType(myGeometry)
            networkType=myGeometry.NetworkType;
        end
        
        function pavageType=GetPavageType(myGeometry)
            pavageType=myGeometry.PavageType;
        end
        
        function nBlock=GetNumberOfBlocks(myGeometry)
            nBlock=myGeometry.NBlock;
        end
        
        function nBoundary=GetNumberOfBoundaries(myGeometry)
            nBoundary=myGeometry.NBoundary;
        end
        
        function nVertice=GetNumberOfVertices(myGeometry)
            nVertice=myGeometry.NVertice;
        end
        
        function ConvertToMeters = GetConvertToMeters(myGeometry)
            ConvertToMeters=myGeometry.ConvertToMeters;
        end
        
        function seed=GetRandomSeed(myGeometry)
            seed=myGeometry.RandomSeed;
        end
        
        function boundaryType=GetBoundaryType(myGeometry,numBoundary)
            boundaryType = myGeometry.Boundaries{numBoundary}.Type;
        end
        
        function vertices=GetBoundaryVertices(myGeometry,numBoundary)
            vertices=myGeometry.Vertices(myGeometry.Boundaries{numBoundary}.Face,:);
        end
        
        function numVertices=GetBoundaryVerticesNumber(myGeometry,numBoundary)
            numVertices = myGeometry.Boundaries{numBoundary}.Face;
        end
        
        
        function name=GetBoundaryName(myGeometry,numBoundary)
            name = myGeometry.Boundaries{numBoundary}.Name;
        end
        
        function rugosity=GetBoundaryRugosity(myGeometry,numBoundary)
            rugosity = myGeometry.Boundaries{numBoundary}.Rugosity;
        end
        
        function vertices=GetAllVertices(myGeometry)
            vertices=myGeometry.Vertices;
        end
        
        function vertice=GetVertice(myGeometry,iVertice)
            vertice=myGeometry.Vertices(iVertice,:);
        end
        
        
        
        function indice_block = FindBlockOfThisBoundary(myGeometry,iBoundary)            
            %Renvoie l'indice du bloc auquel appartient la frontiere
            %numero indice_boundary
            boundaryVerticesNumber = myGeometry.GetBoundaryVerticesNumber(iBoundary);
            block_found = false;
            for iBlock = 1:myGeometry.GetNumberOfBlocks
                if block_found == false
                        blockVerticesNumber = myGeometry.GetBlockVerticesNumber(iBlock);
                    if  isequal(ismember(boundaryVerticesNumber,blockVerticesNumber),ones(1,length(boundaryVerticesNumber)))
                        block_found = true;
                        indice_block = iBlock;
                    end
                end
            end
        end 
        
        
        
        
        function vertices=GetBlockVertices(myGeometry,numBlock)
            vertices=myGeometry.Vertices(myGeometry.Blocks{numBlock}.VertexNumbers,:);
        end
        
        function numVertices=GetBlockVerticesNumber(myGeometry,numBlock)
            numVertices=myGeometry.Blocks{numBlock}.VertexNumbers;
        end
        
        function parameters=GetBlockFillingParameters(myGeometry,numBlock)
            parameters=myGeometry.Blocks{numBlock}.Filling;
        end
        
        function nBand=GetNumberOfAnisotropyBand(myGeometry)
            nBand=myGeometry.NAnisotropyBand;
        end
        
        function direction = GetAnisotropyBandDirection(myGeometry,numBand)
            direction=myGeometry.AnisotropyBand{numBand}.Direction;
        end
                
        function intensity = GetAnisotropyBandIntensity(myGeometry,numBand) 
            intensity=myGeometry.AnisotropyBand{numBand}.Intensity;
        end
                
        function location=GetAnisotropyBandLocation(myGeometry,numBand)
            location=myGeometry.AnisotropyBand{numBand}.Location;
        end
        
        
        %Functions to construct a geometry
        
        function AddVertices(myGeometry,vertices)
            myGeometry.Vertices=vertcat(myGeometry.Vertices,vertices);
            myGeometry.NVertice=myGeometry.NVertice+length(vertices(:,1));
        end
        
        function RemoveVertices(myGeometry,numVertices)
            myGeometry.Vertices=myGeometry.Vertices(setdiff(1:myGeometry.NVertice,numVertices),:);
            myGeometry.NVertice=myGeometry.NVertice-length(numVertices);
        end
        
        function AddBlock(myGeometry,block)
            assert(isa(block,'struct'),'block needs to be a struct')
            myGeometry.Blocks{myGeometry.NBlock+1}=block;
            myGeometry.NBlock=myGeometry.NBlock+1;
        end
        
        function RemoveBlock(myGeometry,numBlock)
            myGeometry.Blocks=myGeometry.Blocks( setdiff(1:myGeometry.NBlock,numBlock) );
            myGeometry.NBlock=myGeometry.NBlock-length(numBoundary);
        end
        
        function AddBoundary(myGeometry,boundary)
            assert(isa(boundary,'struct'),'boundary needs to be a struct')
            myGeometry.Boundaries{end+1}=boundary;
            myGeometry.NBoundary=myGeometry.NBoundary+1;
        end
        
        function RemoveBoundary(myGeometry,numBoundary)
            myGeometry.Boundaries=myGeometry.Boundaries(setdiff( 1:myGeometry.NBoundary,numBoundary ));
            myGeometry.NBoundary=myGeometry.NBoundary-length(numBoundary);
        end
        
        function AddAnisotropyBand(myGeometry,anisotropyBand)
            assert(isa(anisotropyBand,'struct'),'anisotropyBand needs to be a struct');
            myGeometry.AnisotropyBand{end+1}=anisotropyBand;
            myGeometry.NAnisotropyBand=myGeometry.NAnisotropyBand+1;
        end
        
        function RemoveAnisotropyBand(myGeometry,numBoundary)
            myGeometry.AnisotropyBand=myGeometry.AnisotropyBand(setdiff( 1:myGeometry.NBoundary,numBoundary ));
            myGeometry.NAnisotropyBand=myGeometry.NAnisotropyBand-1;
        end 
        
        function SetDimension(myGeometry,dimension)
            myGeometry.Dimension=dimension;
        end
        
        function SetNetworkType(myGeometry,networkType)
            myGeometry.NetworkType=networkType;
        end
        
        function SetPavageType(myGeometry,pavageType)
            myGeometry.PavageType=pavageType;
        end
        
        function SetConvertToMeters(myGeometry,convertToMeters)
            myGeometry.ConvertToMeters=convertToMeters;
        end
        
        function SetRandomSeed(myGeometry,seed)
            myGeometry.RandomSeed=seed;
        end        
        
        %Functions to process a geometry
        
        function ConvertScale(myGeometry)
            
            cvtmrs=myGeometry.GetConvertToMeters;
            
            myGeometry.Vertices = cvtmrs*myGeometry.Vertices;
            
            for iBlock=1:myGeometry.GetNumberOfBlocks
                if isfield(myGeometry.GetBlockFillingParameters(iBlock),'FiberThickness')
                    myGeometry.Blocks{iBlock}.Filling.FiberThickness=cvtmrs*myGeometry.GetBlockFillingParameters(iBlock).FiberThickness;
                end
            end
            
            for iBand=1:myGeometry.GetNumberOfAnisotropyBand
                myGeometry.AnisotropyBand{iBand}.Location=cvtmrs*myGeometry.GetAnisotropyBandLocation(iBand);
            end
            
            myGeometry.SetConvertToMeters(1);
        end 
        
        function ChangeVertice(myGeometry,iVertice,newCoordinates)
            myGeometry.Vertices(iVertice,:)=newCoordinates;
        end
        
        
        %Utilities to handle a geometry

        function myNewGeometry=CopyGeometry(myGeometry)
            
            myNewGeometry=MacroscopicGeometry();
            
            myNewGeometry.Dimension=myGeometry.Dimension;
            myNewGeometry.Vertices=myGeometry.Vertices;
            myNewGeometry.Blocks=myGeometry.Blocks;
            myNewGeometry.Boundaries=myGeometry.Boundaries;
            myNewGeometry.AnisotropyBand=myGeometry.AnisotropyBand;
            myNewGeometry.NVertice=myGeometry.NVertice;
            myNewGeometry.NBlock=myGeometry.NBlock;
            myNewGeometry.NBoundary=myGeometry.NBoundary;
            myNewGeometry.NAnisotropyBand=myGeometry.NAnisotropyBand;
            myNewGeometry.NetworkType=myGeometry.NetworkType;
            myNewGeometry.PavageType=myGeometry.PavageType;
            myNewGeometry.ConvertToMeters=myGeometry.ConvertToMeters;
            myNewGeometry.RandomSeed=myGeometry.RandomSeed;
        end
        
        function View(myGeometry,option,varargin)
            %'Vertices', 'Blocks', 'Boundaries'
            
        end
        
        function LoadGeometry(myGeometry,filename)
            %input : myGeometry,filename
            
            reader=FileReaderXML(filename);
            infos=reader.Read;
            
            geometry = infos.MacroscopicGeometry;
            
            myGeometry.SetDimension(geometry.ATTRIBUTE.Dimension);
            myGeometry.SetNetworkType(geometry.ATTRIBUTE.Type);
            
            if isfield(geometry.ATTRIBUTE,'Pavage')
                myGeometry.SetPavageType(geometry.ATTRIBUTE.Pavage);
            end
            if isfield(geometry.ATTRIBUTE,'ConvertToMeters')
                myGeometry.SetConvertToMeters(geometry.ATTRIBUTE.ConvertToMeters);
            end
            
            myGeometry.AddVertices(geometry.Vertices);
            
            if isfield(geometry,'Block')
                nBlock=length(geometry.Block); 
            else
                nBlock=0;
            end
            for iBlock=1:nBlock
                block=struct;
                block.VertexNumbers = geometry.Block(iBlock).VertexNumbers;
                block.Filling = geometry.Block(iBlock).Filling.ATTRIBUTE;
                myGeometry.AddBlock(block);
            end
            
            if isfield(geometry,'Boundary')
                nBoundary=length(geometry.Boundary);  
            else
                nBoundary=0;
            end
            nBoundary=length(geometry.Boundary);
            for iBoundary=1:nBoundary
                boundary=struct;
                
                boundary.Face = geometry.Boundary(iBoundary).Face ;
                boundary.Name = geometry.Boundary(iBoundary).ATTRIBUTE.Name ;
                boundary.Rugosity = geometry.Boundary(iBoundary).ATTRIBUTE.Rugosity ;
                boundary.Type = geometry.Boundary(iBoundary).ATTRIBUTE.Type ;
                
                myGeometry.AddBoundary(boundary);
            end
            
            if isfield(geometry,'AnisotropyBand')
                nAnisotropyBand=length(geometry.AnisotropyBand);  
            else
                nAnisotropyBand=0;
            end
            for iBand=1:nAnisotropyBand
                band=struct;
                
                band.Direction = geometry.AnisotropyBand(iBand).ATTRIBUTE.Direction ;
                band.Intensity = geometry.AnisotropyBand(iBand).ATTRIBUTE.Intensity ;           
                band.Location = str2num(geometry.AnisotropyBand(iBand).ATTRIBUTE.Location) ;
                
                myGeometry.AddAnisotropyBand(band);
            end
            
        end
        
        function WriteGeometryFile(myGeometry,filename)
            %input : myGeometry,filename
            
            outputStruct=myGeometry.PrivateOutputStruct;
            
            writer=FileWriterXML(filename);
            writer.Write(outputStruct);
            writer.delete;
        end
        
        function outputStruct=PrivateOutputStruct(myGeometry)
            
            outputStruct=struct;
            outputStruct.ATTRIBUTE.Type='MacroscopicGeometry';
            
            geometry=struct;
            geometry.ATTRIBUTE.Dimension=myGeometry.GetDimension;
            geometry.ATTRIBUTE.Type=myGeometry.GetNetworkType;
            geometry.ATTRIBUTE.Pavage=myGeometry.GetPavageType;
            geometry.ATTRIBUTE.ConvertToMeters=myGeometry.GetConvertToMeters;
            
            geometry.Vertices=myGeometry.GetAllVertices;
            
            nBlock=myGeometry.GetNumberOfBlocks;
            geometry.Block(nBlock,1)=struct;
            for iBlock=1:nBlock
                geometry.Block(iBlock).VertexNumbers=myGeometry.GetBlockVerticesNumber(iBlock);
                geometry.Block(iBlock).Filling.ATTRIBUTE=myGeometry.GetBlockFillingParameters(iBlock);
            end
            
            nBoundary=myGeometry.GetNumberOfBoundaries;
            geometry.Boundary(nBoundary,1)=struct;
            for iBoundary=1:nBoundary                
                geometry.Boundary(iBoundary).Face = myGeometry.GetBoundaryVerticesNumber(iBoundary);
                geometry.Boundary(iBoundary).ATTRIBUTE.Name = myGeometry.GetBoundaryName(iBoundary);
                geometry.Boundary(iBoundary).ATTRIBUTE.Rugosity = myGeometry.GetBoundaryRugosity(iBoundary);
                geometry.Boundary(iBoundary).ATTRIBUTE.Type = myGeometry.GetBoundaryType(iBoundary);
            end
            
            nAnisotropyBand=myGeometry.GetNumberOfAnisotropyBand;
            geometry.AnisotropyBand(nAnisotropyBand,1)=struct;
            for iBand=1:nAnisotropyBand
                geometry.AnisotropyBand(iBand).ATTRIBUTE.Direction=myGeometry.GetAnisotropyBandDirection(iBand);
                geometry.AnisotropyBand(iBand).ATTRIBUTE.Intensity=myGeometry.GetAnisotropyBandIntensity(iBand);           
                geometry.AnisotropyBand(iBand).ATTRIBUTE.Location=myGeometry.GetAnisotropyBandLocation(iBand);
            end

            outputStruct.MacroscopicGeometry=geometry;
            
        end
        
        
        
    end
    
end

