function find_steady_state(filename,initfilename,XA,params,funcs)
if ~isfile(filename)
    load(initfilename,'HUxN');
    params.Sigma = 0;

    J = (length(XA)-3)/2;

    dt= 0.01;
    maxnsteps = 1000;
    tol = 1e-2;

    %params = num2cell([dt J delta eps E lambda n m C K 0]);
    % params.K = K;
    % funcs.meltf = @(x) mrate+0*x;
    % funcs.Hsbf = Hsbf;

    ex = 1;
    step = 1;
    change = inf;
    while ( (step<maxnsteps) && (ex>0) &&  (change> tol*dt) && (dt>1e-5))
        [HUxNn, fval, ex] = fsolve( ...
            @(huxn) OBJ_FUNC_ice_hyd(huxn, HUxN([1:J+1 2*J+4:3*J+5]), ...
            XA, dt, funcs, params),...
            HUxN, optimoptions('fsolve','Display','off'));
        change = norm(HUxNn - HUxN);
        if ex>0
            HUxN = HUxNn;
            % If it worked, we proceed
        else
            dt=dt/2;
            ex = 4;
            % If not, we reduce the timestep and try again
        end
        step = step + 1;
    end

    if ex<1
        disp(['Failed to converge on step' num2str(step-1)]);
        figure(8); plot(fval)
    elseif step>=maxnsteps
        disp('Max no. steps reached');
    elseif dt<=1e-5
        disp('Timestep reduced'); % Better way to do this?
    elseif change<tol
        disp('Steady state found'); 
        save(filename,'HUxN','params','funcs','-mat')
    end
end
end