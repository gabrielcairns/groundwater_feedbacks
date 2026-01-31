function fun = OBJ_FUNC_ice_hyd(HUxN, HxNI, XA,  dt, funcs, params)
%fun = OBJ_FUNC_ice_hyd(HUxN, HxNI, XA, bedf, accf, mf, Hdf, dt, params)
% Parameters
%[dt, J, delta, eps, E, lambda, n, m, C, K, Sigma] = deal(params{:});
J = (length(XA)-3)/2;

delta = params.delta;   eps = params.eps;   E = params.E;
lambda = params.lambda;     n = params.n;   mW = params.mW;
C = params.C;   K = params.K;   Sigma = params.Sigma;

bedf = funcs.bedf;  accf = funcs.accf;
meltf = funcs.meltf;    Hsbf = funcs.Hsbf;

% Mesh grids
XH = XA(2:2:2*J+2);
XU = XA(1:2:2*J+3);

% Current values to solve for
H = HUxN(1:J+1);
U = HUxN(J+2:2*J+3);
xg =  HUxN(2*J+4);
N = HUxN(2*J+5:3*J+5);

% Previous value
HI = HxNI(1:J+1);
xgI = HxNI(J+2);
NI = HxNI(J+3:2*J+3);

% Time derivatives
xgt = (xg - xgI)/dt;
dNdt = (N - NI)/dt;

a = accf(xg*XH);
bed = bedf(xg*XH);

% Momentum balance:
Ux = (xg^(-1))*(U(2:J+2)-U(1:J+1))./(XU(2:J+2)-XU(1:J+1));
exstr = H.* (Ux.^2 + 1e-4).^((1-n)/(2*n)).*Ux ;

tauex = 4*eps* (xg^(-1))* (exstr(2:J+1)-exstr(1:J))./(XH(2:J+1)-XH(1:J)) ;
taudr = .5*(H(2:J+1)+H(1:J) ) .*  (xg^(-1)).* (H(2:J+1)-H(1:J) + bed(2:J+1)-bed(1:J))./(XH(2:J+1)-XH(1:J));

Nep = .5*(N(1:J) + N(2:J+1)); % Midpoint method
if isinf(C)
    tausl = U(2:J+1).*abs(U(2:J+1)).^(mW-1);
else 
    tausl = C*Nep.* sign(U(2:J+1)).*(abs(U(2:J+1))./(abs(U(2:J+1)) + (C*Nep).^(1/mW) )).^mW ;
end

% Flux condition
fun(1) = U(1)+U(2); 

% Momentum balance
fun(2:J+1) = tauex - tausl - taudr;

% Stress balance at GL
fun(J+2) = U(J+2)-U(J+1) - xg*(XU(J+2)-XU(J+1))*(delta*H(J+1)/(8*eps)).^n;

% Mass equation:
Hp = [H(1)+bed(2)-bed(1); H]; % minus-oneth grid cell
XHp = [-XH(2); XH];

% Advection (downwind method works best, by experience)
HU = Hp.*U(1:J+2); 
divflux = (xg^(-1))*(HU(2:J+2)-HU(1:J+1))./(XU(2:J+2)-XU(1:J+1)); 

% Accounting for moving grid
GLdrift = (xgt/xg).*XHp(1:J+1).*( (Hp(2:J+2)-Hp(1:J+1))./(XHp(2:J+2)-XHp(1:J+1)));

% Mass balance
dHdt = GLdrift + a - divflux;

% Mass equation
fun(J+3:2*J+3) = H-HI-dt*dHdt;

% Flotation condition
fun(2*J+4) = H(J+1)+bed(J+1)/(1-delta);

% Hydrology equation
% Hydraulic potential
Psi = H - E*N + bed/(1-delta) ;

% Flux (encapsulating no flux, but including melt into first cell)
flux1 = (xg^(-1))*(1./N(1:J) + K*Hsbf(XU(2:J+1)) ).* (Psi(2:J+1)-Psi(1:J))./((XH(2:J+1)-XH(1:J)));
flux = [xg*(XU(2)-XH(1))*( meltf(xg*XH(1)) - Sigma*Hsbf(xg*XH(1))*(dHdt(1) - GLdrift(1)) ) ; flux1];

% Accounting for moving grid
GLdriftN = (xgt/xg).*XH(1:J).*( (N(2:J+1)-N(1:J))./(XH(2:J+1)-XH(1:J)));

fun(2*J+5:3*J+4) =  (-lambda*N(1:J).^(-2) -  E*Sigma*Hsbf(xg*XH(1:J)) ).*(dNdt(1:J)-GLdriftN) - ...
    (meltf(xg*XH(1:J)) - Sigma*Hsbf(xg*XH(1:J)).*(dHdt(1:J) - GLdrift(1:J)) ...
    +(xg^(-1))*(flux(2:J+1)-flux(1:J))./(XU(2:J+1)-XU(1:J)) );

% Zero effective pressure at GL
fun(3*J+5) = N(J+1);
end