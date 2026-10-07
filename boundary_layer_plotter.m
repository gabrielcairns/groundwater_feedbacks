% Auxiliary function to plot cases from the parameter sweep, zoomed in on
% the boundary layer
% Can take filename explicitly, or this field can be left empty and a
% default case is chosen
function boundary_layer_plotter(filename,mrate,K)

if isempty(filename)
    filename = sprintf('steady_states/SS_m%g_K%g',mrate,K);
end

[scales, params, funcs] = default_values();
Hd = scales.Hd;   xd_km = scales.xd_km; Nd_MPa = scales.Nd_MPa; ud_myr = scales.ud_myr;
delta = params.delta;  E = params.E;  
params.Sigma = 0;
delta = params.delta;   eps = params.eps;   E = params.E;
lambda = params.lambda;     n = params.n;   mW = params.mW;
C = params.C;   


bedf = funcs.bedf;
bedxf = funcs.bedxf;
Hsbf = funcs.Hsbf;

load("variable_grid_300.mat","XH","XU");
J = 300;

funcs.meltf = @(x) mrate*ones(size(x));

HUxN = load(filename,'HUxN').HUxN;

H = HUxN(1:J+1);
U = HUxN(J+2:2*J+3);
xg =  HUxN(2*J+4);
N = HUxN(2*J+5:3*J+5);

bed = bedf(xg*XH);

% Momentum balance for the ice
Ux = (xg^(-1))*(U(2:J+2)-U(1:J+1))./(XU(2:J+2)-XU(1:J+1));
exstr = H.* (Ux.^2 + 1e-4).^((1-n)/(2*n)).*Ux ; 
tauex = 4*eps* (xg^(-1))* (exstr(2:J+1)-exstr(1:J))./(XH(2:J+1)-XH(1:J)) ;
taudr = .5*(H(2:J+1)+H(1:J) ) .*  (xg^(-1)).* (H(2:J+1)-H(1:J) + bed(2:J+1)-bed(1:J))./(XH(2:J+1)-XH(1:J));
Nsl = .5*(N(1:J) + N(2:J+1)); 
tausl = C*Nsl.* sign(U(2:J+1)).*(abs(U(2:J+1))./(abs(U(2:J+1)) + (C*Nsl).^(1/mW) )).^mW ;

% Approximate solutions for N
P = HUxN(1:J+1)+ bedf(xg*XH)/(1-delta);
G = xg^(-1)*(P(2:J+1)-P(1:J))./(XH(2:J+1)-XH(1:J));
Nx = xg^(-1)*(N(2:J+1)-N(1:J))./(XH(2:J+1)-XH(1:J));
Nouter = -G./(xg*mrate*XU(2:J+1) - K*Hsbf(xg*XU(2:J+1)).*G);
Ninner = P/E;

figure(1); clf;
aH = make_nice_axes(subplot(5,1,1),'$x-x_g$','$z$'); 
atau = make_nice_axes(subplot(5,1,2),'$x-x_g$','Stress');
asl = make_nice_axes(subplot(5,1,3),'$x-x_g$','$\tau_b$');
aN = make_nice_axes(subplot(5,1,4),'$x-x_g$','$N$');
aNeq = make_nice_axes(subplot(5,1,5),'$x-x_g$','');

plot(aH, xg*(XH-1), H + bed ,'LineWidth', 2 );
plot(aH, xg*(XH-1), bed ,'k' ,'LineWidth', 2);

plot(atau, xg*(XU(2:J+1)-1),  -tauex,'LineWidth', 2 );
plot(atau, xg*(XU(2:J+1)-1),  -taudr,'LineWidth', 2 );
plot(atau, xg*(XU(2:J+1)-1),  tausl,'LineWidth', 2 );
legend(atau,{'$-\tau_{ex}$','$-\tau_{dr}$','$\tau_{sl}$'},'FontSize',14,'Interpreter','latex');

plot(asl, xg*(XU(2:J+1)-1),  tausl,'LineWidth', 2 );
plot(asl, xg*(XU(2:J+1)-1),  C*Nsl,'LineWidth', 2 );
plot(asl, xg*(XU(2:J+1)-1),  abs(U(2:J+1)).^mW,'LineWidth', 2 );
legend(asl,{'$\tau_{sl}$','$CN$','${u_i}^{m_W}$'},'FontSize',14,'Interpreter','latex');

plot(aN, xg*(XH-1), N,'LineWidth', 2 );
plot(aN, xg*(XU(2:J+1)-1), Nouter,'LineWidth', 2 );
plot(aN, xg*(XH(end-30:end)-1), Ninner(end-30:end),'LineWidth', 2 );
legend(aN,{'$N$','$E=0$','$\psi=0$'},'FontSize',14,'Interpreter','latex');

plot(aNeq,xg*(XU(2:J+1)-1), G,'LineWidth', 2 );
plot(aNeq,xg*(XU(2:J+1)-1), E*Nx,'LineWidth', 2 );
plot(aNeq,xg*(XH-1), mrate*xg*XH.*N./(1+K*Hsbf(xg*XH).*XH),'LineWidth', 2 );

% Internal axis plotting function
function ax = make_nice_axes(ax,xlab,ylab)
    ax.XScale = 'log';
    xlim([-0.2 0]);
    cla;
    hold(ax,'on');
    set(ax,'FontSize',14,'TickLabelInterpreter','latex');
    xlabel(ax,xlab,'FontSize',14,'Interpreter','latex');
    ylabel(ax,ylab,'FontSize',14,'Interpreter','latex');
end
end