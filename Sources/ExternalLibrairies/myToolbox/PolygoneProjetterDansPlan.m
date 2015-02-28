function coordonnees_planes=PolygoneProjetterDansPlan(coordonnees_polygone)
%ProjetterPolygoneDansPlan Renvoie les coordonn�es 2D des sommets d'un
%   polygone plan plonge dans l'espace 3D.
%   input: - coordonnees_polygone=array(nombre_sommets,dimension)
    
    dimension=length(coordonnees_polygone(1,:));
    assert(dimension==3);
    centre=mean(coordonnees_polygone);
    %coordonnees_polygone=transpose(coordonnees_polygone);
    %centre=transpose(centre);
    vect1=coordonnees_polygone(1,:)-centre;
    vect2=coordonnees_polygone(2,:)-centre;
    a = vect1/norm(vect1);
    b = vect2/norm(vect2);
    normal = [a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]; %cross product
    
    axe_rotation= [normal(2),-normal(1),0] ;  % cross(normal,[0 0 1])
    if norm(axe_rotation)>1e-6
        axe_rotation=axe_rotation\norm(axe_rotation);
        a=axe_rotation(1);
        b=axe_rotation(2);
        c=axe_rotation(3);
        th=acos(c);
        M=[a^2*(1-cos(th))+cos(th)  a*b*(1-cos(th))-c*sin(th) a*c*(1-cos(th))+b*sin(th); ...
           a*b*(1-cos(th))+c*sin(th) b^2*(1-cos(th))+cos(th)   b*c*(1-cos(th))-a*sin(th); ...
            a*c*(1-cos(th))-b*sin(th) b*c*(1-cos(th))+a*sin(th)  c^2*(1-cos(th))+cos(th)];

        foo=transpose(M*transpose(coordonnees_polygone));
        coordonnees_planes=horzcat(foo(:,1),foo(:,2));
    else
        coordonnees_planes=horzcat(coordonnees_polygone(:,1),coordonnees_polygone(:,2));
    end
    
end

