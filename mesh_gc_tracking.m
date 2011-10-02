function [labels,borders,maxflow,de,se_space,se_time,face_labels] = mesh_gc_tracking(points,faces,data,thresh,lambda_space,lambda_time,options)
%   MESH_GC_TRACKING   Compute tracking on a 3D triangulated mesh
%
%       [LABELS,BORDERS,MAXFLOW,DE,SE_SPACE,SE_TIME] = MESH_GC_TRACKING(POINTS,FACES,DATA,THRESH,LAMBDA_SPACE,LAMBDA_TIME,OPTIONS)
%
%   Compute tracking on a 3D triangulated mesh using a binary cut of a graph
%
%   Created by Alexandre Gramfort on 2008-05-28.
%   Copyright (c) 2007 Alexandre Gramfort. All rights reserved.
%

% $Id: mesh_gc_tracking.m 171 2009-10-22 13:23:06Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-10-22 15:23:06 +0200 (Thu, 22 Oct 2009) $
% $Revision: 171 $

me = 'MESH_GC_TRACKING';

% ===========
% = Options =
% ===========

if nargin<7
    options.null = 0;
end

if ~isfield(options, 'weight_type')
    options.weight_type = 2; % use face weighted areas and edge length for graph weights
end
weight_type = options.weight_type;

if ~isfield(options, 'verbose')
    options.verbose = true;
end
verbose = options.verbose;

if nargin == 0
    eval(['help ',lower(me)])
    options = rmfield(options,'null')
    return
end

% ========
% = Core =
% ========

npoints = size(data,1);
nwin = size(data,2);

data_cost = thresh*ones([npoints,2,nwin]);
data_cost(:,2,:) = data;

if verbose
    disp(['---- Using tracking graph-cut']);
end
[labels,maxflow,borders,de,se_space,se_time,face_labels] = mesh_graph_cut(points, faces, data_cost, lambda_space, lambda_time, weight_type, verbose);

se = maxflow-de;;
if verbose
    disp(sprintf('  maxflow  =  data_term + lambda_space *  se_space  + lambda_time * se_time'));
    disp(sprintf('%01.2e =  %01.2e +      %01.2e  * %01.2e +    %01.2e   *  %01.2e', ...
        [maxflow de lambda_space se_space lambda_time se_time]))
    disp(sprintf('%01.2e =  %01.2e +                %01.2e          +                %01.2e', ...
        [maxflow, de, lambda_space*se_space, lambda_time*se_time]))
    % disp(['---- Energy  : ',num2str(maxflow)]);
    % disp(['---- DE + SE : ',num2str(de,'%10.3e'),' + ',num2str(se,'%10.3e')]);
end

end % function