[scales, params, funcs] = default_values();

params.lambda = 1e-3;

mrates  = 0:0.25:12;
Ks = 0:0.25:6;

% mrates  = 6:0.25:8.75;
% Ks = 0:0.25:10;

mrates  = 0:0.25:20;
for l = 2:length(mrates)
    funcs.meltf = @(x) mrates(l)*ones(size(x));
    find_steady_state(sprintf('steady_states/SS_m%g_K0.mat', mrates(l)), sprintf('steady_states/SS_m%g_K0.mat', mrates(l-1)),XA,params,funcs);
end

for l = 1:length(mrates)
funcs.meltf = @(x) mrates(l)*ones(size(x));
load("variable_grid_300.mat","XA");

for ll = 2:length(Ks)
    params.K = Ks(ll);  

    initfilename= sprintf('steady_states/SS_m%g_K%g.mat', mrates(l), Ks(ll-1));
    if isfile(initfilename)
    find_steady_state(sprintf('steady_states/SS_m%g_K%g.mat', mrates(l), Ks(ll)), initfilename, XA,params,funcs);
    end    

end
end

for l = 1:length(mrates)
funcs.meltf = @(x) mrates(l)*ones(size(x));
load("variable_grid_300.mat","XA");

for ll = 2:length(Ks)
    params.K = Ks(ll);  

    initfilename= sprintf('steady_states/SS_m%g_K%g.mat', mrates(l), Ks(ll-1));
    if isfile(initfilename)
    find_steady_state(sprintf('steady_states/SS_m%g_K%g.mat', mrates(l), Ks(ll)), initfilename, XA,params,funcs);
    end    

end
end