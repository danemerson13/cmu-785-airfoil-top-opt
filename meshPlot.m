function h4 = meshPlot(ElementList, NodeList, ElementList_Active_Bool, NodeList_FixedDOF_Bool, fillValue)
    %Setup data to send to quadplot
    %ElementList assumed to be [m x 5]. m is number of elements
    %[element#, node#1, node#2, node#3, node#4]
    % NodePos is assumed to be [n x 3] and sorted in numerical order by node number.
    % n is number of nodes
    % [node#, x#, y#]
    
    % meshPlot(ElementList,NodeList)
    % meshPlot(ElementList,NodeList, ElementList_Active_Bool, NodeList_FixedDOF_Bool)
    
    quad = ElementList(:, 2:end);
    x = NodeList(:,2);
    y = NodeList(:,3);
    
    if (nargin == 2)
        quadplot(quad,x,y)
    elseif (nargin < 5)
        quadplot(quad,x,y, 'k', ElementList_Active_Bool, NodeList_FixedDOF_Bool);
    else
%         quadplot(quad,x,y, 'k', ElementList_Active_Bool, NodeList_FixedDOF_Bool, fillValue, clim, colormap)
        h4 = quadplot(quad,x,y, 'k', ElementList_Active_Bool, NodeList_FixedDOF_Bool, fillValue);

    end
end

function h4 = quadplot(quad, varargin)
%TRIPLOT Plots a 2D triangulation
%   QUADPLOT(QUAD,X,Y) displays the quadrilaterals defined in the
%   M-by-4 matrix QUAD.  A row of QUAD contains indices into X,Y that
%   define a single quadrilateal. The default line color is blue.
%
%   QUADPLOT(...,COLOR) uses the string COLOR as the line color.
%
%   H = QUADPLOT(...) returns a vector of handles to the displayed 
%   quadrilaterals
%
%   QUADPLOT(...,'param','value','param','value'...) allows additional
%   line param/value pairs to be used when creating the plot.
%
%   See also TRISURF, TRIMESH, DELAUNAY, TriRep, DelaunayTri.
%
%   Script code based on copyrighted code from mathworks for TRIPLOT.
%   Allan P. Engsig-Karup, apek@imm.dtu.dk.
error(nargchk(1,inf,nargin,'struct'));
% start = 1;
x = varargin{1};
y = varargin{2};
quads = quad;
if (nargin == 3)
    c = 'blue';
    ElementList_Active_Bool = [];
    NodeList_FixedDOF_Bool = [];
%     start = 3;
else
    c = varargin{3};
    ElementList_Active_Bool = varargin{4};
    NodeList_FixedDOF_Bool = varargin{5};
%     start = 4;
end
  
if (nargin > 6)
    fillValue = varargin{6};
%     clim = varargin{7};
%     colormap = varargin{8};
end

d = quads(:,[1 2 3 4 1]);

if (nargin > 3)
    e = d(logical(ElementList_Active_Bool(:,2)),:);
end

d = d';

if (nargin > 3)
    e = e';
    f = logical(NodeList_FixedDOF_Bool(:,2));
end

% h = plot(x(d), y(d),c,varargin{start:end});
h = plot(x(d), y(d),c);

if (nargin > 3)
    hold on
    h2 = plot(x(e), y(e),'red');
    h3 = plot(x(f), y(f),'*');
end

h4 = 0;
if (nargin > 6)
    h4 = patch(x(d), y(d), fillValue);
end

axis equal
xlim([-0.025, 0.45])
ylim([-0.1, 0.1])
hold off

if nargout == 1, hh = h; end
end