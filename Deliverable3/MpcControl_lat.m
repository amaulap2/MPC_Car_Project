classdef MpcControl_lat < MpcControlBase
    
    methods
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   x0           - initial state (estimate)
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   u0           - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_segs = ceil(mpc.H/mpc.Ts); % Horizon steps
            N = N_segs + 1;              % Last index in 1-based Matlab indexing

            [nx, nu] = size(mpc.B);
            
            % Targets
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);

            % Initial states
            x0 = sdpvar(nx, 1);
            x0other = sdpvar(nx, 1); % (Ignore this, not used)

            % Input to apply to the system
            u0 = sdpvar(nu, 1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D
            %       are the DISCRETE-TIME MODEL of your system
            %       You can find the linearization steady-state
            %       in mpc.xs and mpc.us.
            
            % System dynamics and steady-state offsets
            A = mpc.A;
            B = mpc.B;
            xs = mpc.xs;
            us = mpc.us; 
            
            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            Q = 10*eye(nx); % Weight for state tracking error
            R = 2*eye(nu);  % Weight for control effort


            T_settle = 3; % Settling time in seconds
            N_settle = ceil(T_settle / mpc.Ts); % Number of steps for settling

            % Prediction horizon constraints and objective
            M = [1;-1]; 
            m = [1; 1];
            m_y = [3.5;0.5];
            m_theta = [deg2rad(5);deg2rad(5)];
            m_delta = [deg2rad(30);deg2rad(30)];

            slack = sdpvar(N-1, 1);  % Slack variable for input constraint
            
            F = [1 0; 0 1; -1 0; 0 -1]; f = [m_y(1); m_theta(1); m_y(2); m_theta(1)];

            [K,Qf,~] = dlqr(A,B,Q,R);
            K = -K; 

            % Compute maximal invariant set
            Xf = polytope([F;M*K],[f;m]);
            Acl = [A+B*K];
            while 1
                prevXf = Xf;
                [T,t] = double(Xf);
                preXf = polytope(T*Acl,t);
                Xf = intersect(Xf, preXf);
                if isequal(prevXf, Xf)
                    break
                end
            end
            [Ff,ff] = double(Xf);

            x = sdpvar(nx,N,'full');
            u = sdpvar(nu,N-1,'full');

            % Initialize constraints and objective
            con = (x(:,1) == x0); % Initial state constraint
            con = con + (u0 == u(:,1)); % Initial input constraint
            con = con + (slack(1) >= 0);
            con = con + (M * x(1,1) <= m_y) + (M * x(2,1) <= m_theta + slack(1));
            obj = (x(:,1) - x_ref)' * Q * (x(:,1) - x_ref); % Initial state cost
            obj = obj + 1000 * slack(1)^2;  % Penalize constraint violation

            % Add constraints and objective for the prediction horizon
            for i = 1:N-1
                % Dynamics constraints
                con = con + (x(:,i+1) == A * x(:,i) + B * u(:,i));
                % State constraints
                con = con + (slack(i) >= 0);
                con = con + (M * x(1,i) <= m_y) + (M * x(2,i) <= m_theta + slack(i) );
                obj = obj + 10000 * sum(slack(i)^2);
                % Input constraints
                con = con + (M * u(:,i) <= m_delta);
                % Cost function
                % Attention us
                obj = obj + (x(:,i) - x_ref)' * Q * (x(:,i) - x_ref) + (u(:,i) - u_ref)' * R * (u(:,i) - u_ref);
            end

            con = con + (Ff*x(:,N) <= ff);
            % Terminal cost
            obj = obj + (x(:,N) - x_ref)' * Qf * (x(:,N) - x_ref);


            % Replace this line and set u0 to be the input that you
            % want applied to the system. Note that u0 is applied directly
            % to the nonlinear system. You need to take care of any 
            % offsets resulting from the linearization.
            % If you want to use the delta formulation make sure to
            % substract mpc.xs/mpc.us accordingly.

            % Pass here YALMIP sdpvars which you want to debug. You can
            % then access them when calling your mpc controller like
            % [u, X, U] = mpc_lon.get_u(x0, ref);
            % with debugVars = {X_var, U_var};
            debugVars = {x, u};
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {x0, x_ref, u_ref, x0other}, {u0, debugVars{:}});
        end
        
        % Computes the steady state target which is passed to the
        % controller
        function [xs_ref, us_ref] = compute_steady_state_target(mpc, ref)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            % OUTPUTS
            %   xs_ref, us_ref - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Steady-state subsystem
            A = mpc.A(2, 2);
            B = mpc.B(2, 1);

            % Subsystem linearization steady-state
            xs = mpc.xs;
            us = mpc.us;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % Reference velocity (steady-state)
            xs_ref = zeros(2,1);
            xs_ref(1) = ref; % always y_ref
            xs_ref(2) = 0; % theta needs to be 0 to drive straight
            % Compute steady-state input
            if abs(B) > 1e-6 % Avoid division by zero
                us_ref = (1 - A) * xs_ref(2) / B;
            else
                us_ref = 0; % Default to zero if dynamics are degenerate
            end
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end
