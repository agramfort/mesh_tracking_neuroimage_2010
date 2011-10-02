function [dists] = dists_relative(template,data)

% DISTS_RELATIVE   Compute relative distances between a template and data points
%
%   
%   
%   SYNTAX
%       [DISTS] = DISTS_RELATIVE(TEMPLATE,DATA)
%   

%
%   Created by Alexandre Gramfort on 2008-03-27.
%   Copyright (c) 2007 Alexandre Gramfort. All rights reserved.
%

energy_max = norm(template);
dists = repmat(template,1,size(data,1))-data';
dists = sqrt(sum(dists .* dists)) / energy_max; % take relative distance
