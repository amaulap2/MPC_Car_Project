classdef MpcControl_lon < MpcControlBase
    
    methods
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   x0           - initial state (estimate)
            %   V_ref, u_ref - reference state/input
            %   d_est        - disturbance estimate
            %   x0other      - initial state of other car
            % OUTPUTS
            %   u0           - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_segs = ceil(mpc.H/mpc.Ts); % Horizon steps
            N = N_segs + 1;              % Last index in 1-based Matlab indexing

            [nx, nu] = size(mpc.B);
            
            % Targets
            V_ref = sdpvar(1);
            u_ref = sdpvar(1);

            % Disturbance estimate (Ignore this before Todo 4.1)
            d_est = sdpvar(1);

            % Initial states
            x0 = sdpvar(nx, 1);
            x0other = sdpvar(nx, 1); % (Ignore this before Todo 5.1)

            % Input to apply to the system
            u0 = sdpvar(nu, 1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

            % System dynamics and steady-state offsets
            A = mpc.A;
            B = mpc.B;
            xs = mpc.xs;
            us = mpc.us; 
            
            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            Q = [0,0;0,100]; % Weight for state tracking error
            R = 10*eye(nu);  % Weight for control effort
            
            T_settle = 10; % Settling time in seconds
            N_settle = ceil(T_settle / mpc.Ts); % Number of steps for settling

            % Prediction horizon constraints and objective
            M = [1;-1];
            m = [1;1];

            x = sdpvar(nx,N,'full');
            u = sdpvar(nu,N-1,'full');

            % Initialize constraints and objective
            con = (x(:,1) == x0); % Initial state constraint
            con = con + (u0 == u(:,1));
            obj = (x(:,1) - V_ref)' * Q * (x(:,1) - V_ref); % Initial state cost

            % Add constraints and objective for the prediction horizon
            for i = 1:N-1
                % Dynamics constraints
                u(:,i) = u(:,i);
                x(:,i) = x(:,i);
                con = con + (x(:,i+1) == xs + A * (x(:,i)-xs) + B * (u(:,i)-us) + B * d_est);
                % Input constraints
                con = con + (M * u(:,i) <= m);
                % Cost function
                obj = obj + (x(:,i) - V_ref)' * Q * (x(:,i) - V_ref) + (u(:,i) - u_ref)' * R * (u(:,i) - u_ref);
            end

            % Terminal cost
            obj = obj + (x(:,N) - V_ref)' * Q * (x(:,N) - V_ref);

            % [u, X, U] = mpc_lon.get_u(x0, ref);
            % with debugVars = {X_var, U_var};
            debugVars = {x,u};

            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {x0, V_ref, u_ref, d_est, x0other}, {u0, debugVars{:}});
        end
        
        % Computes the steady state target which is passed to the
        % controller
        function [Vs_ref, us_ref] = compute_steady_state_target(mpc, ref, d_est)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            %   d_est  - disturbance estimate (Ignore before Todo 4.1)
            % OUTPUTS
            %   Vs_ref, us_ref - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Steady-state subsystem
            A = mpc.A(2, 2);
            B = mpc.B(2, 1);

            % Subsystem linearization steady-state
            xs = mpc.xs(2);
            us = mpc.us;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % Reference velocity (steady-state)
            Vs_ref = ref; % Steady-state velocity matches the reference velocity

            % Compute steady-state input
            if abs(B) > 1e-6 % Avoid division by zero
                us_ref = (1-A)*(Vs_ref-xs)/B + us - d_est;
            else
                us_ref = us; % Default to zero if dynamics are degenerate
            end
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end
