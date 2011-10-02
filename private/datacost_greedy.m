function [datacost,templates_idx] = datacost_greedy(M,options)
%   DATACOST_GREEDY   Short description
%       [DATACOST] = DATACOST_GREEDY(M,)
% 
%   Long description
%   
%   Created by Alexandre Gramfort on 2008-06-02.
%   Copyright (c) 2007 Alexandre Gramfort. All rights reserved.

if nargin<2
    options.null = 0;
end

if ~isfield(options, 'use_corr')
    options.use_corr = true;
end
use_corr = options.use_corr;

if ~isfield(options, 'pct_max_decrease')
    options.pct_max_decrease = 0.2;
end
pct_max_decrease = options.pct_max_decrease;

if ~isfield(options, 'debug_mode')
    options.debug_mode = false;
end
debug_mode = options.debug_mode;

if ~isfield(options, 'nb_max_atoms')
    options.nb_max_atoms = 1;
end
nb_max_atoms = options.nb_max_atoms;

%%%%%%%%%

energies = sqrt(sum(M .* M,2));
[energy_max,template_idx] = max(energies);
templates_idx = template_idx;
template_energy = energy_max;
orig_template = M(template_idx,:)';

if use_corr
    dists = dists_corr(orig_template,M);
else
    dists = dists_relative(orig_template,M);
end

last_template_normalized = orig_template ./ norm(orig_template);

last_energy_max = energy_max;
for k=1:(nb_max_atoms-1)
    M = M - (M * last_template_normalized) * last_template_normalized';
    energies = sqrt(sum(M .* M,2));
    [new_energy_max,new_template_idx] = max(energies);
    new_template = M(new_template_idx,:)';
    if new_energy_max > (pct_max_decrease * energy_max)
        disp(['    Adding one template, pct_decrease : ',num2str(fix(100*new_energy_max/energy_max))]);

        if use_corr
            new_dists = dists_corr(new_template,M);
        else
            new_dists = dists_relative(new_template,M);
        end

        templates_idx = [templates_idx,new_template_idx];

        last_energy_max = new_energy_max;
        dists = min(dists,new_dists); % take relative distance
        last_template_normalized = new_template ./ norm(new_template);
    else
        break
    end
end

datacost = dists;

end %  function