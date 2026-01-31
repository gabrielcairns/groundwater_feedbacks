function [scales, params, funcs] = default_values()

s_per_yr = 365.25 * 24 * 3600; % seconds per year
scales.g = 9.8; 
scales.rhoi = 900; 
rhow = 1000; % gravity, density
%eta = 1e-3; % Pa s

xd = 5e5; 
scales.xdkm = xd/1e3; % lengthscale

ad = 0.3/s_per_yr; % accumulation scale


Cw = 7.624e6; 
%A = 4.227e-25;
A = 1e-25; 

params.mW = 1/3; % weertman law
params.n = 3; % glen's law

scales.md_mmyr = 10;
md = scales.md_mmyr /(1e3*s_per_yr) ; % 10 mm/yr 

scales.Hd = (Cw*xd^(params.mW+1)*ad^params.mW/(scales.rhoi*scales.g)).^(1/(params.mW+2));

scales.td = (Cw*xd^(params.mW+1)/(scales.rhoi*scales.g*ad^2)).^(1/(params.mW+2));  
scales.td_kyr = scales.td/(1e3*s_per_yr);

ud = xd/scales.td;  
scales.ud_myr = ud*s_per_yr;

params.delta = 1-scales.rhoi/rhow;
params.eps = (2*scales.rhoi*scales.g*scales.Hd*(A*scales.td)^(1/params.n))^(-1);


params.E = 0.05;

Nd = params.E*(scales.rhoi*scales.g*scales.Hd);
scales.Nd_MPa = Nd*1e-6;
% Nd = Ed/hd;
Ed = 1e5; % Made up - could be 1e6
hd = Ed/Nd;

%hd = (eta*xd^2*md)/(k*scales.rhoi*g*scales.Hd);

params.lambda = hd/(md*scales.td);

% To be decided whether we include some conversion factors for K and Sigma
% into kd and Ss, and whether we include implied values

%k = 0.05*eta*xd^2*md/Ed;
%k = (eta*xd^2*md)/(hd*scales.rhoi*scales.g*scales.Hd);

params.C = 5;
%Cc = (scales./scales.E)*scales.Hd/xd;



%kd = 1.6e-12;
%kd = k*scales.Hd/scales.Hd;

%K = kd*scales.Hd/(k*scales.Hd);

%xi = 0.2;
%Ss = 1e-4;
%Sigma = (scales.rhoi*Ss*(1-xi)*scales.Hd^2)/(rhow*md*td);

arate = 1;

funcs.bedf = @(x) -100/scales.Hd - .1*x;  
funcs.bedxf = @(x) - .1*ones(size(x));  
funcs.Hsbf = @(x) 2000/scales.Hd*ones(size(x));

funcs.accf = @(x) arate*ones(size(x));
funcs.intaccf = @(x) arate*x;

end