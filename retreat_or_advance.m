function retreat_or_advance(mode,mrate,K,Sigma)

% Runs either the grounding line retreat or advance experiment, and saves
% the result under a specified filename

arguments
    mode {mustBeMember(mode,["a","r"])} 
    mrate {mustBeNumeric}
    K {mustBeNumeric}
    Sigma {mustBeNumeric}
end  

if ~isfolder('retreat_advance'); mkdir('retreat_advance'); end

% Filename under which solution will be saved
if mode == 'r'
    filename = sprintf('retreat_advance/retreat_m%g_K%g_Sig%g',mrate,K,Sigma);
    ICfilename = sprintf('retreat_advance/max_xg_m%g_K%g',mrate,K);
elseif mode == 'a'
    filename = sprintf('retreat_advance/advance_m%g_K%g_Sig%g',mrate,K,Sigma);
    ICfilename = sprintf('retreat_advance/min_xg_m%g_K%g',mrate,K);
end

if ~isfile([filename '.mat'])

    % Setup
    options = optimoptions('fsolve','Display','off');
    
    [scales, params, funcs] = default_values();
    Hd = scales.Hd;     rhoi = scales.rhoi;     g = scales.g;
    td = scales.td;     n = params.n;  
    % delta = params.delta;   C = params.C;   E = params.E;  
    
    % Topography
    % bedf = @(x) -100/Hd - .1*x;  
    % bedxf = @(x)  -.1 +0*x;
    % Hdf = @(x) 2000/Hd + 0*x;

    % bedf = funcs.bedf;
    % bedxf = funcs.bedxf;
    % Hsbf = funcs.Hsbf;

    params.K = K;
    params.Sigma = Sigma;
    funcs.meltf = @(x) mrate*ones(size(x));    
    funcs.intmeltf = @(x) mrate*x;

    load("variable_grid_300.mat","XA","XH","XU");
    J = 300;
    
    Amin = 1e-25;
    Amax = 2e-25;
    
    %% Creation of initial and final states
    params.eps = (2*rhoi*g*Hd*(Amin*td)^(1/n))^(-1);
    if ~isfile(sprintf('retreat_advance/max_xg_m%g_K%g.mat',mrate,K))
        if isfile(sprintf('retreat_advance/max_xg_m%g_K%g.mat',mrate-.5,K))
            initfilename = sprintf('retreat_advance/max_xg_m%g_K%g.mat',mrate-.5,K);
        elseif isfile(sprintf('retreat_advance/max_xg_m%g_K%g.mat',mrate,K-.5))
            initfilename = sprintf('retreat_advance/max_xg_m%g_K%g.mat',mrate,K-.5);
        else
            %HUxNw = load_or_make([ICfilename '.mat'],[ICfilename '_rough.mat']);
            initfilename = sprintf('retreat_advance/max_xg_m%g_K%g_rough.mat',mrate,K);
            if ~isfile(initfilename)
                make_rough_steady_state(params,funcs,XH,XU,sprintf('retreat_advance/max_xg_m%g_K%g',mrate,K));
            end
        end
        find_steady_state(sprintf('retreat_advance/max_xg_m%g_K%g.mat',mrate,K),initfilename,XA,params,funcs);
    end

    params.eps = (2*rhoi*g*Hd*(Amax*td)^(1/n))^(-1);
    if ~isfile(sprintf('retreat_advance/min_xg_m%g_K%g.mat',mrate,K))
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
        find_steady_state(sprintf('retreat_advance/min_xg_m%g_K%g.mat',mrate,K),initfilename,XA,params,funcs);
    end

    if mode == 'r'
      %  dt = 2.5e-4;
        Afunc = @(t) Amax + 0*t;
    elseif mode == 'a'
     %   dt = 1e-3;
        Afunc = @(t) Amin + 0*t;
    end 
    
    %dt = 1e-4;
    %t = 0:dt:2;
    t= linspace(0,sqrt(2),200).^2;

    HUxNw = load([ICfilename '.mat'],'HUxN').HUxN;

    % Preallocating arrays in which to record the answer
    %kgif = 1:10:length(t);
    kgif = 1:length(t);

    percents = round(linspace(length(t)/10,length(t),10));
    Hgif = nan(J+1,length(kgif));
    Ugif = nan(J+2,length(kgif));
    Ngif = nan(J+1,length(kgif));
    xggif = nan(1,length(t));
    
    % Record the initial state
    Hgif(:,1) = HUxNw(1:J+1);
    Ugif(:,1) = HUxNw(J+2:2*J+3);    
    Ngif(:,1) = HUxNw(2*J+5:3*J+5);
    xggif(1) = HUxNw(2*J+4);

    ex = 1;

    params.K = K;
    params.Sigma = Sigma;
    funcs.meltf = @(x) mrate + 0*x;
    for k = 1:length(t)-1

        % A / eps as a smooth function
        A = Afunc(t(k));
        params.eps = (2*rhoi*g*Hd*(A*td)^(1/n))^(-1);

        if ex>0
            HxNI =  HUxNw([1:J+1 2*J+4:3*J+5]);
            [HUxNw, ~, ex] = fsolve( ...
                @(huxn) OBJ_FUNC_ice_hyd(huxn, HxNI, ...
                XA, t(k+1)-t(k), funcs, params),...
                HUxNw, options);

            kc = k;
            
            % Rcord all the values
                xggif(k+1) = HUxNw(2*J+4);
            if any(k==kgif)
                Hgif(:,k==kgif) = HUxNw(1:J+1);
                Ugif(:,k==kgif) = HUxNw(J+2:2*J+3);
                Ngif(:,k==kgif) = HUxNw(2*J+5:3*J+5);
            end
            if any(k==percents)
                disp([num2str(10*find(k==percents)) '% complete']);
            end
            
            % Special exit message for negative N
            if any(HUxNw(2*J+5:3*J+5)<-1e-6) % Allow some tolerance
                ex = -5; % Special error for negative N
            end
        end
    end

    % After exiting the loop
    if ex<1    
        % If it errored, give error message
        fprintf('Failed to converge on step %g, xg=%g, t=%g with exit flag%g\n',kc,xggif(kc),t(kc),ex);
    else   
        % Otherwise save 
        disp('Converged at all steps');
        save([filename '.mat'],'Hgif','Ugif','Ngif','xggif','t','kgif') 
    end
end

end