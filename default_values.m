function [scales, params, funcs] = default_values()
% Creates structs containing the default scalings, parameter values and
% functions (b, Hsb etc.) to use in our solutions
% Values correspond to SI units unless denoted otherwise with a suffix (e.g. _km, _mmyr etc.)

% Number of seconds in a year (for unit conversions)
s_per_yr = 365.25 * 24 * 3600;

% Gravity and densities
scales.g = 9.8; 
scales.rhoi = 900; 
rhow = 1000; 

% Lengthscale
xd = 5e5; 
scales.xd_km = xd/1e3;

% Accumulation scale (0.3 m /yr)
ad = 0.3/s_per_yr; 

% Weertman law
Cw = 7.624e6; 
params.mW = 1/3; 

% Glen's law
A = 1e-25; 
params.n = 3; 

% Melt rate (10 mm/yr)
scales.md_mmyr = 10;
md = scales.md_mmyr /(1e3*s_per_yr) ;

% Finding scales derived from ice sheet model
scales.Hd = (Cw*xd^(params.mW+1)*ad^params.mW/(scales.rhoi*scales.g)).^(1/(params.mW+2));

scales.td = (Cw*xd^(params.mW+1)/(scales.rhoi*scales.g*ad^2)).^(1/(params.mW+2));  
scales.td_kyr = scales.td/(1e3*s_per_yr);

ud = xd/scales.td;  
scales.ud_myr = ud*s_per_yr;

% Parameters from ice sheet model
params.delta = 1-scales.rhoi/rhow;
params.eps = (2*scales.rhoi*scales.g*scales.Hd*(A*scales.td)^(1/params.n))^(-1);

% Parameters and scales for the subglacial hydrology
params.E = 0.05;

Nd = params.E*(scales.rhoi*scales.g*scales.Hd);
scales.Nd_MPa = Nd*1e-6;
Ed = 1e5;
hd = Ed/Nd;

params.lambda = hd/(md*scales.td);
params.C = 5;

% Functions for geometry of bed and sedimentary basin
funcs.bedf = @(x) -100/scales.Hd - .1*x;  
funcs.bedxf = @(x) - .1*ones(size(x));  
funcs.Hsbf = @(x) 2000/scales.Hd*ones(size(x));

% Accumulation functions (default to dimensionless a=1)
funcs.accf = @(x) ones(size(x));
funcs.intaccf = @(x) x;

end