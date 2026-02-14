function make_rough_steady_state(params,funcs,XH,XU,ICfilename)
% Code to create an approximate steady state using asymptotic solutions

% Load the parameters and functions from provided structs
delta = params.delta;   eps = params.eps;   E = params.E;
n = params.n;   mW = params.mW;
C = params.C;   K = params.K;  

bedf = funcs.bedf;  
bedxf = funcs.bedxf; 
Hsbf = funcs.Hsbf;
intaccf = funcs.intaccf;
intmeltf = funcs.intmeltf;    

% Create an initial guess for xg using Tsai approximation (see Tsai et al. 2015)
fxg = @(x) ( 0.61*delta^(n-1)*E/( 8^(n-1)*C*eps^n*(1-delta)^(n+2) ) ) * (-bedf(x)).^(n+2) - intaccf(x);
xg_R = fsolve(@(x) fxg(x), 2, optimoptions('fsolve','Display','off')); 

% Preallocate arrays for H and N, which we will find using direct integration
H_R = ones(size(XH));
N_R = zeros(size(XH));

% Floatation boundary condition for H
H_R(end) =  -bedf(xg_R)/(1-delta);
% (Note N=0 at grounding line already enforced by our preallocation)

% Integration from xg inwards
for k = length(XH):-1:2

    % Integrating ODE for H, ignoring extensional stress terms
    HxT = -(intaccf(xg_R*XH(k))).^mW./H_R(k).^(mW+1) - bedxf(xg_R*XH(k)) ;
    H_R(k-1) = H_R(k)- xg_R*(XH(k)-XH(k-1))*HxT;

    % Integrating ODE for N given H
    NxT = ( HxT + bedxf(xg_R*XH(k))/(1-delta) + ( N_R(k)./(1 + K*Hsbf(xg_R*XH(k))*N_R(k)) ).*intmeltf(xg_R*XH(k)) )/E ;
    N_R(k-1) = N_R(k)- xg_R*(XH(k)-XH(k-1))*NxT;
end

% U found from ice flux
U_R = [0; intaccf(xg_R*XU(2:end))./H_R];

% Concatenate our rough solution into a vector
HUxN = [H_R; U_R; xg_R; N_R];

% Save under a filename that reminds us this is a rough approximation
save([ICfilename '_rough.mat'],'HUxN');
end