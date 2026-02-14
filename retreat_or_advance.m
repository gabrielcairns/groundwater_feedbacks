function retreat_or_advance(mode,mrate,K,Sigma)

% Runs either the grounding line retreat or advance experiment for given parameter 
% values m, K, Sigma, and saves the result in an appropriately named file

% If a file whose name matches the desired case already exists, this code does nothing

arguments
    mode {mustBeMember(mode,["a","r"])} 
    mrate {mustBeNumeric}
    K {mustBeNumeric}
    Sigma {mustBeNumeric}
end  

% Folder in which solutions are saved
if ~isfolder('retreat_advance'); mkdir('retreat_advance'); end

% Filename under which solution will be saved, and name of file to be
% loaded for initial condition
if mode == 'r'
    filename = sprintf('retreat_advance/retreat_m%g_K%g_Sig%g',mrate,K,Sigma);
    ICfilename = sprintf('retreat_advance/max_xg_m%g_K%g',mrate,K);
elseif mode == 'a'
    filename = sprintf('retreat_advance/advance_m%g_K%g_Sig%g',mrate,K,Sigma);
    ICfilename = sprintf('retreat_advance/min_xg_m%g_K%g',mrate,K);
end


if ~isfile([filename '.mat'])

    % Use default parameter values and functions, other than those specified
    [scales, params, funcs] = default_values();
    Hd = scales.Hd;     rhoi = scales.rhoi;     g = scales.g;
    td = scales.td;     n = params.n;  

    % Update the parameter and function structs with our desired values
    params.K = K;
    params.Sigma = Sigma;
    funcs.meltf = @(x) mrate*ones(size(x));    
    funcs.intmeltf = @(x) mrate*x;

    % Meshgrid
    load("variable_grid_300.mat","XA","XH","XU");
    J = 300;
    
    % The values of A which we switch between
    Amin = 1e-25;
    Amax = 2e-25;
    
    % Find the steady state for the max xg / minimum A, if it doesn't
    % already exists
    params.eps = (2*rhoi*g*Hd*(Amin*td)^(1/n))^(-1);
    if ~isfile(sprintf('retreat_advance/max_xg_m%g_K%g.mat',mrate,K))
        % Finding an initial condition. We first try cases where m or K is
        % a little lower, or failing that we load or make a rough guess from
        % scratch
        if isfile(sprintf('retreat_advance/max_xg_m%g_K%g.mat',mrate-.5,K))
            initfilename = sprintf('retreat_advance/max_xg_m%g_K%g.mat',mrate-.5,K);
        elseif isfile(sprintf('retreat_advance/max_xg_m%g_K%g.mat',mrate,K-.5))
            initfilename = sprintf('retreat_advance/max_xg_m%g_K%g.mat',mrate,K-.5);
        else
            initfilename = sprintf('retreat_advance/max_xg_m%g_K%g_rough.mat',mrate,K);
            if ~isfile(initfilename)
                make_rough_steady_state(params,funcs,XH,XU,sprintf('retreat_advance/max_xg_m%g_K%g',mrate,K));
            end
        end
        
        % Having obtained an initial condition, we find the steady state
        find_steady_state(sprintf('retreat_advance/max_xg_m%g_K%g.mat',mrate,K),initfilename,XA,params,funcs);
    end
    
    % Same as above, but for min xg/ maximum A
    params.eps = (2*rhoi*g*Hd*(Amax*td)^(1/n))^(-1);
    if ~isfile(sprintf('retreat_advance/min_xg_m%g_K%g.mat',mrate,K))
        % Finding an initial condition
        if isfile(sprintf('retreat_advance/min_xg_m%g_K%g.mat',mrate-.5,K))
            initfilename = sprintf('retreat_advance/min_xg_m%g_K%g.mat',mrate-.5,K);
        elseif isfile(sprintf('retreat_advance/min_xg_m%g_K%g.mat',mrate,K-.5))
            initfilename = sprintf('retreat_advance/min_xg_m%g_K%g.mat',mrate,K-.5);
        else
            initfilename = sprintf('retreat_advance/min_xg_m%g_K%g_rough.mat',mrate,K);
            if ~isfile(initfilename)
                make_rough_steady_state(params,funcs,XH,XU,sprintf('retreat_advance/min_xg_m%g_K%g',mrate,K));
            end
        end
        % Having obtained an initial condition, we find the steady state
        find_steady_state(sprintf('retreat_advance/min_xg_m%g_K%g.mat',mrate,K),initfilename,XA,params,funcs);
    end

    % One of the states above will be our initial condition, depending on
    % whether we are advancing or retreating

    % Value of A over time
    if mode == 'r'
        Afunc = @(t) Amax + 0*t;
    elseif mode == 'a'
        Afunc = @(t) Amin + 0*t;
    end 
    
    % Variable time meshgrid (much smaller steps initially)
    t= linspace(0,sqrt(2),200).^2;

    % Load initial condition
    HUxNw = load([ICfilename '.mat'],'HUxN').HUxN;

    % Preallocating arrays in which to record the solution
    ks = 1:length(t);
    Hs = nan(J+1,length(ks));
    Us = nan(J+2,length(ks));
    Ns = nan(J+1,length(ks));
    xgs = nan(1,length(t));
    
    % For progress updates
    percents = round(linspace(length(t)/10,length(t),10));

    % Record the initial state
    Hs(:,1) = HUxNw(1:J+1);
    Us(:,1) = HUxNw(J+2:2*J+3);    
    Ns(:,1) = HUxNw(2*J+5:3*J+5);
    xgs(1) = HUxNw(2*J+4);

    ex = 1;
    for k = 1:length(t)-1

        % Update epsilon for this timestep
        A = Afunc(t(k));
        params.eps = (2*rhoi*g*Hd*(A*td)^(1/n))^(-1);

        % We continue while fsolve is converging. If at any point it fails,
        % we note the step and move on
        if ex>0

            % Timestep forward using fsolve
            [HUxNw, ~, ex] = fsolve( ...
                @(huxn) OBJ_FUNC_ice_hyd(huxn, HUxNw([1:J+1 2*J+4:3*J+5]), ...
                XA, t(k+1)-t(k), funcs, params),...
                HUxNw,  optimoptions('fsolve','Display','off'));

            % Note which step we're on, in case it fails
            kc = k;
            
            % Record all the values
            xgs(k+1) = HUxNw(2*J+4);
            Hs(:,k==ks) = HUxNw(1:J+1);
            Us(:,k==ks) = HUxNw(J+2:2*J+3);
            Ns(:,k==ks) = HUxNw(2*J+5:3*J+5);
            
            % Give updates on progress
            if any(k==percents)
                disp([num2str(10*find(k==percents)) '% complete']);
            end
            
            % We also create a special error flag if N becomes negative somehow
            if any(HUxNw(2*J+5:3*J+5)<-1e-6) % (Allowing some tolerance)
                ex = -5; 
            end
        end
    end

    % After exiting the loop
    if ex<1    
        % If it failed at a certain timestep, display an error message
        fprintf('Failure creating %s: failed to converge on step %g, xg=%g, t=%g with exit flag%g \n',filename,kc,xgs(kc),t(kc),ex);
    else   
        % Otherwise, print a success message and save the solution
        fprintf('Success creating %s: converged at all steps \n',filename);
        save([filename '.mat'],'Hs','Us','Ns','xgs','t','ks') 
    end
else
    
    % Output message if the file already exists
    fprintf('File %s already exists, no action taken \n',filename)    
end

end