function sommets_ordonnes=PolygoneOrdonnerSommetsDansPlan(coordonnees_planes)
%ORDONNERSOMMETSPOLYGONEDANSPLAN Ordonne les sommets d'un polygone convexe 
%   2D dans l'ordre trigonométrique.
%   input : -coordonnees_planes: array(nombre_sommets_du_polygone,2), 
%            coordonnees X,Y des sommets d'un polygone convexe  
%   output : - sommets_ordonnes : array(1,nombre_sommets_du_polygone), num
%              des sommets ordonnes dans l'ordre trigonometrique
%Algo : on part du sommet S1 le plus à gauche. On prend S0 à
%   la verticale de S1 sur le bas.On calcule les produits scalaires 
%   <S0S1,SiS1> renormalisés. Ils varient de 1 à -1 lorsqu'on fait le tour
%   des sommets dans l'ordre trigonométrique (tous les sommets sont au
%   dessus de S0S1). En ordonnant les produits scalaires on ordonne les
%   sommets.
    
    [~,index]=sortrows(coordonnees_planes);
    S1=coordonnees_planes(index(1),:);
    S0=[S1(1),(S1(2)-1)];
    
    nSommets=length(coordonnees_planes(:,1));
    produits_scalaire=zeros(1,nSommets);
    produits_scalaire(index(1))=-1;
    for iSommet=1:nSommets
        if iSommet~=index(1)
            vect=coordonnees_planes(iSommet,:)-S1;
            vect=vect/norm(vect);
            produits_scalaire(iSommet)=dot((S1-S0),vect);
        end
    end
    [~,ind]=sort(produits_scalaire);
    sommets_ordonnes=ind;
end

