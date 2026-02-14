% Objective function to solve the equations for ice and subglacial
% hydrology
function fun = OBJ_FUNC_ice_hyd(HUxN, HxNI, XA,  dt, funcs, params)
% HUxN = [H; U; xg; N] contains the variables to be found at current step
% HxNI = [HI; xgI; NI] contains initial values i.e. previous timestep
% XA is the meshgrid (interleaving whole points XH with half points XU)
% dt is the timestep
% funcs and params are structs containing functions (b, Hsb, a, m etc.) and
% parameters used in finding the solution

% Number of cells (minus one)
J = (length(XA)-3)/2;

% Parameters
delta = params.delta;   eps = params.eps;   E = params.E;
lambda = params.lambda;     n = params.n;   mW = params.mW;
C = params.C;   K = params.K;   Sigma = params.Sigma;

% Functions
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

% Initial values
HI = HxNI(1:J+1);
xgI = HxNI(J+2);
NI = HxNI(J+3:2*J+3);

% Time derivatives
xgt = (xg - xgI)/dt;
dNdt = (N - NI)/dt;

a = accf(xg*XH);
bed = bedf(xg*XH);

%% Ice momentum balance
% Stresses are found at half mesh points (XU)

% Extensional stress
Ux = (xg^(-1))*(U(2:J+2)-U(1:J+1))./(XU(2:J+2)-XU(1:J+1));
exstr = H.* (Ux.^2 + 1e-4).^((1-n)/(2*n)).*Ux ; 
% (Includes regularisation to avoid NaN error when Ux = 0)
tauex = 4*eps* (xg^(-1))* (exstr(2:J+1)-exstr(1:J))./(XH(2:J+1)-XH(1:J)) ;

% Driving stress (midpoint method)
taudr = .5*(H(2:J+1)+H(1:J) ) .*  (xg^(-1)).* (H(2:J+1)-H(1:J) + bed(2:J+1)-bed(1:J))./(XH(2:J+1)-XH(1:J));

% Sliding stress
Nsl = .5*(N(1:J) + N(2:J+1)); % Midpoint method
if isinf(C)
    % (We only use finite C, but included for completeness)
    tausl = U(2:J+1).*abs(U(2:J+1)).^(mW-1);
else 
    tausl = C*Nsl.* sign(U(2:J+1)).*(abs(U(2:J+1))./(abs(U(2:J+1)) + (C*Nsl).^(1/mW) )).^mW ;
end

% No flux condition
fun(1) = U(1)+U(2); 

% Momentum balance
fun(2:J+1) = tauex - tausl - taudr;

% Stress balance at grounding line
fun(J+2) = U(J+2)-U(J+1) - xg*(XU(J+2)-XU(J+1))*(delta*H(J+1)/(8*eps)).^n;

%% Ice mass balance

% Variables extended into minus-oneth grid cell
Hp = [H(1)+bed(2)-bed(1); H]; 
XHp = [-XH(2); XH];

% Advection terms (downwind method)
HU = Hp.*U(1:J+2); 
divflux = (xg^(-1))*(HU(2:J+2)-HU(1:J+1))./(XU(2:J+2)-XU(1:J+1)); 

% Term due to time-dependent grid stretching
GLdrift = (xgt/xg).*XHp(1:J+1).*( (Hp(2:J+2)-Hp(1:J+1))./(XHp(2:J+2)-XHp(1:J+1)));

% Mass balance
dHdt = GLdrift + a - divflux;
fun(J+3:2*J+3) = H-HI-dt*dHdt;

% Flotation condition
fun(2*J+4) = H(J+1)+bed(J+1)/(1-delta);

%% Subglacial hydrology equation

% Hydraulic potential
Psi = H - E*N + bed/(1-delta) ;

% Flux (encapsulating no flux, but including melt into first cell)
flux1 = (xg^(-1))*(1./N(1:J) + K*Hsbf(XU(2:J+1)) ).* (Psi(2:J+1)-Psi(1:J))./((XH(2:J+1)-XH(1:J)));
flux = [xg*(XU(2)-XH(1))*( meltf(xg*XH(1)) - Sigma*Hsbf(xg*XH(1))*(dHdt(1) - GLdrift(1)) ) ; flux1];

% Term due to time-dependent grid stretching
GLdriftN = (xgt/xg).*XH(1:J).*( (N(2:J+1)-N(1:J))./(XH(2:J+1)-XH(1:J)));

% Mass balance
fun(2*J+5:3*J+4) =  (-lambda*N(1:J).^(-2) -  E*Sigma*Hsbf(xg*XH(1:J)) ).*(dNdt(1:J)-GLdriftN) - ...
    (meltf(xg*XH(1:J)) - Sigma*Hsbf(xg*XH(1:J)).*(dHdt(1:J) - GLdrift(1:J)) ...
    +(xg^(-1))*(flux(2:J+1)-flux(1:J))./(XU(2:J+1)-XU(1:J)) );

% Zero effective pressure at grounding line
fun(3*J+5) = N(J+1);
end