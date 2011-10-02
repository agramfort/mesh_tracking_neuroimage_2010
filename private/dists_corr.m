function [dists] = dists_corr(template,data)

% DISTS_CORR   Compute distance between template and data points using correlation
%
%   
%   
%   SYNTAX
%       [DISTS] = DISTS_CORR(TEMPLATE,DATA)
%   

%
%   Created by Alexandre Gramfort on 2008-03-27.
%   Copyright (c) 2007 Alexandre Gramfort. All rights reserved.
%

% template = template;
% dists = repmat(template,1,size(data,1)) .* data';
% % dists = 2 .* (1 - abs(sum(dists)) ./ norm(template)^2); % take relative distance
% dists = (1 - abs(sum(dists)) ./ norm(template)^2); % take relative distance

[dists] = dists_relative(template,data);
dists(dists > 1) = 2 - dists(dists > 1);
