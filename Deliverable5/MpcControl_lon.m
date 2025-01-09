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
            load('Q.mat','Q')
            load('R.mat','R')
            load('K.mat','K')
            load('Qf','Qf')
            
            load('X_tight.mat','X_tight')
            load('U_tight.mat', 'U_tight')
            load('x_safe.mat','x_safe')
            load('Xf.mat','Xf')
            load('Eps.mat', 'Eps')

            F = X_tight.A;
            f = X_tight.b;
            M = U_tight.A;
            m = U_tight.b;
            FF = Xf.A;
            Ff = Xf.b;

            x = sdpvar(nx,N,'full');
            x_lead = sdpvar(nx,N,'full');

            Delta = sdpvar(nx,N,'full');

            u = sdpvar(nu,N-1,'full');
            mu_tube = sdpvar(nu,N-1,'full');


            % Initialize constraints and objective
            con = (x(:,1) == x0);
            con = con + (x_lead(:,1) == x0other);
            % Slide 33 chap 8 
            FEps = Eps.A;
            fEps = Eps.b;
            con = con + (FEps * ((x_lead(:,1)-x(:,1)-[x_safe;0]) - Delta(:,1)) <= fEps);  % <-->

            con = con + (mu_tube(:,1) == u0);
            con = con + (mu_tube(:,1) == K*((x_lead(:,1)-x(:,1)-[x_safe;0]) - Delta(:,1)) + u(:,1));

            obj = (Delta(:,1))' * Q * (Delta(:,1)) + u(:,1)*R*u(:,1); % Initial state cost

            % Add constraints and objective for the prediction horizon
            for i = 1:N-1
                con = con + (F * Delta(:,i) <= f);
                con = con + (Delta(:,i+1) == A*(Delta(:,i)) - B * u(:,i));
                % Input constraints
                con = con + (M * u(:,i) <= m);
                % Cost function
                obj = obj + (Delta(:,i))' * Q * (Delta(:,i)) + u(:,i)*R*u(:,i);
            end

            % Terminal cost
            con = con + (FF * Delta(:,N) <= Ff);
            obj = obj + Delta(:,N)' * Qf * Delta(:,N);

            % [u, X, U] = mpc_lon.get_u(x0, ref);
            % with debugVars = {X_var, U_var};
            debugVars = {};

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
                us_ref = (1-A)*(Vs_ref-xs)/B + us;
            else
                us_ref = us; % Default to zero if dynamics are degenerate
            end
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end
