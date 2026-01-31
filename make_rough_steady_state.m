function make_rough_steady_state(params,funcs,XH,XU,ICfilename)

delta = params.delta;   eps = params.eps;   E = params.E;
n = params.n;   mW = params.mW;
C = params.C;   K = params.K;  


bedf = funcs.bedf;  
bedxf = funcs.bedxf; 
Hsbf = funcs.Hsbf;
intaccf = funcs.intaccf;
intmeltf = funcs.intmeltf;    

% Create an initial guess using Tsai approx
fxg = @(x) ( 0.61*delta^(n-1)*E/( 8^(n-1)*C*eps^n*(1-delta)^(n+2) ) ) * (-bedf(x)).^(n+2) - intaccf(x);
xg_R = fsolve(@(x) fxg(x), 2, optimoptions('fsolve','Display','off')); 

% Asymptotic H and N by direct integration
H_R = ones(size(XH));
N_R = zeros(size(XH));

H_R(end) =  -bedf(xg_R)/(1-delta);

for k = length(XH):-1:2
    HxT = -(intaccf(xg_R*XH(k))).^mW./H_R(k).^(mW+1) - bedxf(xg_R*XH(k)) ;
    H_R(k-1) = H_R(k)- xg_R*(XH(k)-XH(k-1))*HxT;

    NxT = ( HxT + bedxf(xg_R*XH(k))/(1-delta) + ( N_R(k)./(1 + K*Hsbf(xg_R*XH(k))*N_R(k)) ).*intmeltf(xg_R*XH(k)) )/E ;
    N_R(k-1) = N_R(k)- xg_R*(XH(k)-XH(k-1))*NxT;
end
U_R = [0; intaccf(xg_R*XU(2:end))./H_R];

HUxN = [H_R; U_R; xg_R; N_R]; % Have to declare the variable

save([ICfilename '_rough.mat'],'HUxN');
end