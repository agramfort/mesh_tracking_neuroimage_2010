function [labels,maxflow,borders,de,se_space,se_time,face_labels] = mesh_graph_cut(points, faces, data_cost, lambda_space, lambda_time, weight_type,verbose)
%
%   Performing Graph Cut energy minimization operations on a 3D mesh with time
%
%   Usage:
%       [LABELS,MAXFLOW,BORDERS,DE,SE_SPACE,SE_TIME] = 
%               MESH_GRAPH_CUT(POINTS, FACES, DATA_COST, LAMBDA_SPACE, LAMBDA_TIME, WEIGHT_TYPE);
% 
%   weight_type == 0 : no weighting
%   weight_type == 1 : point based weighting
%   weight_type == 2 : face based weighting
%
%   DE is the data term in the energy
%   SE is obtained with SE = MAXFLOW - DE
% 
%   SE_SPACE : space term without lambda_space
%   SE_TIME : time term without lambda_time
% 
%   We should have :
%       MAXFLOW = DE + SE
%               = DE + LAMBDA_SPACE*SE_SPACE + LAMBDA_TIME*SE_TIME
%
%   This wrapper for Matlab was written by Alexandre Gramfort (firsname.lastname@sophia.inria.fr).
%   Odyssee ENPC/INRIA/Ens Ulm
%   http://www-sop.inria.fr/odyssee/team/
%
%	The core cpp application was written by Vladimir Kolmogorov
%   http://www.adastral.ucl.ac.uk/~vladkolm/software.html
%
%   Matlab Wrapper for Graph Cut.
%        Alexandre Gramfort.
%        in http://www-sop.inria.fr/members/Alexandre.Gramfort/, April 2007.
%
%   This software can be used only for research purposes, you should  cite ALL of
%   the aforementioned papers in any resulting publication.
%
%   The Software is provided "as is", without warranty of any kind.
%

% $Id: mesh_graph_cut.m 171 2009-10-22 13:23:06Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-10-22 15:23:06 +0200 (Thu, 22 Oct 2009) $
% $Revision: 171 $

if nargin < 7
    verbose = true;
end
verbose = single(verbose);

npoints = size(points,1);
if any( size(points,2) ~= [3] )
    size(points,2)
    error('mesh_graph_cut: points has incorrect size');
end

% faces
if any( size(faces,2) ~= [3] )
    error('mesh_graph_cut: faces has incorrect size');
end
if ( (min(faces(:)) <= 0) || (max(faces(:) > npoints)) )
    min(faces(:))
    max(faces(:))
    npoints
    error('mesh_graph_cut: faces have incorrect indexes')
end

% Data cost
if ndims(data_cost) ~= 3 & ndims(data_cost) ~= 2
    error('mesh_graph_cut: data cost has incorect dimensionality');
end
nlabels = size(data_cost,2);
ntimes = size(data_cost,3);
if size(data_cost,1) ~= npoints
    size(data_cost,1)
    npoints
    error('mesh_graph_cut: data cost has incorrect size')
end

if weight_type>=2
    nfaces = size(faces,1);
    data_cost = reshape(data_cost(faces,:,:),[nfaces 3 2 ntimes]);
    data_cost = squeeze(sum(data_cost,2) ./ 3);
end

% Weight data cost
nnodes = npoints; % number of nodes per time frame
if weight_type==3
    areas = ones(nfaces,1);
    nnodes = nfaces;
elseif weight_type==2
    areas = mesh_areas(points,faces);
    nnodes = nfaces;
elseif weight_type==1
    areas = mesh_cell_areas(points,faces);
else % weight_type==0
    areas = ones(nnodes,1);
end

% ===============================
% = Hack : for scale invariance =
% ===============================
scale_invariance = true;
if scale_invariance
    areas = sqrt(areas);
end

data_cost = data_cost .* repmat(areas(:),[1,nlabels,ntimes]);
data_cost = single(data_cost);
data_cost = permute(data_cost,[2 1 3]);
data_cost = data_cost(:);

% Space Smoothness param
if any( size(lambda_space) ~= [1 1] )
    error('mesh_graph_cut: lambda_space has incorrect size');
end
lambda_space = single(lambda_space);

% Time Smoothness param
if any( size(lambda_time) ~= [1 1] )
    error('mesh_graph_cut: lambda_time has incorrect size');
end

if weight_type > 0
    time_edges = areas;
else
    time_edges = ones(nnodes,1);
end

weighted_time_edges = lambda_time .* time_edges;
weighted_time_edges = single(weighted_time_edges(:));

% Edges
if weight_type>=2 
    space_edges = mesh_edge_faces(points,faces);
    if weight_type==3
        space_edges(:,3) = 1;
    end
elseif weight_type==1
    space_edges = mesh_edge_weighted(points,faces);
else % weight_type==0
    [ii,jj] = find(mesh_edges(faces));
    space_edges = [ii,jj,ones(length(ii),1)];
end
gidx = find(space_edges(:,1)>space_edges(:,2));
space_edges = space_edges(gidx,:);
nedges = size(space_edges,1);

weighted_space_edges = space_edges;
weighted_space_edges(:,1:2) = weighted_space_edges(:,1:2) - 1;
weighted_space_edges(:,3) = lambda_space.*weighted_space_edges(:,3);
weighted_space_edges = single(weighted_space_edges');
weighted_space_edges = weighted_space_edges(:);

% no negative dataterm (capacities are all positive)
data_cost = data_cost-min(data_cost(:));

[labels,maxflow] = mesh_graph_cut_mex(nnodes, nedges, ntimes, data_cost, weighted_space_edges, weighted_time_edges,verbose);

% compute dataterm in the energy
de = sum(data_cost(1:2:end).*(labels(:)==0))+sum(data_cost(2:2:end).*(labels(:)==1));

% Removing constant from datacost to make sure that the smallest datacost is 0
maxflow = maxflow - sum(min(data_cost(2:2:end),data_cost(1:2:end)));
de = de - sum(min(data_cost(2:2:end),data_cost(1:2:end)));

% Compute time regularity of the solution
se_time = 0;
diff_labels = abs(diff(labels,1,2));
for k=1:size(diff_labels,2)
    se_time = se_time + sum(time_edges(find(diff_labels(:,k))));
end

% Compute space regularity of the solution
se_space = 0;
for k=1:size(labels,2)
    se_space = se_space + sum(space_edges(:,3).*labels(space_edges(:,1),k).*(1-labels(space_edges(:,2),k)));
end

if weight_type < 2
    borders = mesh_isolevels(faces,labels);
else
    borders = mesh_face_isolevels(faces,labels);
end

face_labels = [];
if weight_type>=2
    face_labels = labels;
    [vertex_numbering,I] = sort(double(faces(:))); % sorted point numbers
    faces_numbering = rem(I-1,nfaces)+1; % triangle number for each Vertex
    face_to_point = sparse(vertex_numbering,faces_numbering,1,npoints,nfaces);
    labels = face_to_point*labels;
    labels = labels > 0;
end

