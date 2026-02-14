function find_steady_state(filename,initfilename,XA,params,funcs)
% Solves the steady state for given parameters and functions (b, Hsb, a, m etc.), provided in
% the structs 'params' and 'funcs'
% The solution is saved under 'filename'
% 'initfilename' should be a .mat file containing HUxN=[H; U; xg; N], which
% is used as the initial value
% XA is the meshgrid (whole points XH interspersed with half points XU)

% If there's already a file with the name 'filename', this code does nothing

% We find the steady state by time-evolving, with a variable timestep, until the change in the norm of
% HUxN is less than a specified tolerance (0.01) times the timestep.
% If the solution fails to converge despite reducing the timestep beyond a
% minimum (1e-5), or a max number of steps (1000), we give up

if ~isfile(filename)

    % Load initial value
    load(initfilename,'HUxN');
    params.Sigma = 0;

    J = (length(XA)-3)/2;
    
    % Initial timestep, max no. of steps, tolerance, minimum timestep
    dt= 0.01;
    maxnsteps = 1000;
    tol = 1e-2;
    mindt = 1e-5;

    step = 1;
    change = inf;

    % Evolve forward in time until an exit criterion is reached
    while ( (step<maxnsteps) &&  (change> tol*dt) && (dt>mindt))

        % Attempt ind the new values, given current ones
        % Record the exit flag of fsolve (i.e. whether it worked)
        [HUxNn, ~, ex] = fsolve( ...
            @(huxn) OBJ_FUNC_ice_hyd(huxn, HUxN([1:J+1 2*J+4:3*J+5]), ...
            XA, dt, funcs, params),...
            HUxN, optimoptions('fsolve','Display','off'));

        % Calculate change in norm
        change = norm(HUxNn - HUxN);

        if ex>0
            % If fsolve worked, we update the values and proceed
            HUxN = HUxNn;
        else
            % If not, we reduce the timestep and try again
            dt=dt/2;
        end

        % Update number of steps
        step = step + 1;
    end
    
    % Once we exit the loop, print a message saying why
    if step>=maxnsteps
        fprintf('Failed creating %s: max number of steps reached \n',filename);
    elseif dt<=mindt
        fprintf('Failed creating %s: timestep reduced beyond minimum tolerance \n',filename); 
    elseif change<tol
        fprintf('Success creating %s \n', filename); 

        % If we were successful, we save the solution
        % We also save the parameters and functions used to find the
        % solution for reproducibility
        save(filename,'HUxN','params','funcs','-mat')
    end
else
    
    % Output message if the file already exists
    fprintf('File %s already exists, no action taken \n',filename)
end
end