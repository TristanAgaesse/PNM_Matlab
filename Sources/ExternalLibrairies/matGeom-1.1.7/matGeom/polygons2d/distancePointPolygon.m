function varargout = distancePointPolygon(point, poly)
%DISTANCEPOINTPOLYGON  Compute shortest distance between a point and a polygon
%
%   DIST = distancePointPolygon(POINT, POLYGON)
%   Computes the shortest distance between the point POINT and the polygon
%   boundary POLYGON. POINT is a 1-by- row vector, and POLYGON is a N-by-2
%   array containing vertex coordinates.
%
%   Example
%   distancePointPolygon
%
%   See also
%   polygons2d, points2d, distancePointPolyline
%   distancePointEdge, projPointOnPolyline
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2009-04-30,    using Matlab 7.7.0.471 (R2008b)
% Copyright 2009 INRA - Cepia Software Platform.

% eventually copy first point at the end to ensure closed polygon
if sum(poly(end, :) == poly(1,:)) ~= 2
    poly = [poly; poly(1,:)];
end

% call to distancePointPolyline 
minDist = distancePointPolyline(point, poly);

% process output arguments
if nargout <= 1
    varargout{1} = minDist;
end
