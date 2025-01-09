classdef NmpcControl < handle

    properties
        % The NMPC problem
        opti

        % Problem parameters
        x0, ref, x0other

        % Most recent problem solution
        sol

        % The input that you want to apply to the system
        u0

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add any variables you would like to read to debug here
        % and then store them in the NmpcControl function below.
        % e.g., you could place X here and then add obj.X = X
        % in the NmpcControl function below.
        % 
        % After solving the problem, you can then read these variables 
        % to debug via
        %   nmpc.sol.value(nmpc.X)
        % 
        X, U
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end

    methods
        function obj = NmpcControl(car, H)

            import casadi.*

            N_segs = ceil(H/car.Ts); % Horizon steps
            N = N_segs + 1;          % Last index in 1-based Matlab indexing

            nx = 4;
            nu = 2;

            % Define the NMPC optimization problem
            opti = casadi.Opti();
            
            % Parameters (symbolic)
            obj.x0 = opti.parameter(nx, 1);       % initial state
            obj.ref = opti.parameter(2, 1);       % target y, velocity
            obj.x0other = opti.parameter(nx, 1);  % initial state of other car

            % SET THIS VALUE TO BE YOUR CONTROL INPUT
            obj.u0 = opti.variable(nu, 1);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

            % Define your problem using the opti object created above
            % Variables
            X = opti.variable(nx, N); 
            U = opti.variable(nu, N - 1); 
            slack = opti.variable(1,N);

            % Vehicle parameters
            lr = car.lr;
            lf = car.lf; 
            Pmax = car.Pmax;
            Cd = car.Cd;
            Af = car.Af; 
            Cr = car.Cr; 
            m = car.mass;
            rho = car.rho; 
            g = car.g; 

            opti.subject_to(X(:, 1) == obj.x0);
            opti.subject_to( obj.u0 == U(:,1) );
            
            cost = 0;

            Q = 5*eye(2);
            R = 10*eye(nu);

            % Dynamics and constraints
            for k = 1:N-1
                beta = atan((lr * tan(U(1, k))) / (lr + lf));

                Fmotor = (U(2, k) * Pmax) / max(X(4, k), 0.1); % Prevent division by zero
                Fdrag = 0.5 * rho * Cd * Af * X(4, k)^2;
                Froll = Cr * m * g;

                x_dot = [  
                    X(4, k) * cos(X(3, k) + beta); 
                    X(4, k) * sin(X(3, k) + beta);       
                    X(4, k) / lr * sin(beta);          
                    (Fmotor - Fdrag - Froll) / m      
                    ];

                opti.subject_to(X(:, k+1) == X(:, k) + car.Ts * x_dot);

                opti.subject_to(-0.5 <= X(2, k)); % y constraints
                opti.subject_to(X(2, k) <= 3.5);
                opti.subject_to(-(deg2rad(5)+slack(k)) <= X(3, k)); % theta constraints
                opti.subject_to(X(3, k) <= deg2rad(5)+slack(k));
                opti.subject_to(-1 <= U(2, k)); % u_T constraints
                opti.subject_to(U(2, k) <= 1);
                opti.subject_to(-deg2rad(30) <= U(1, k)); % delta constraints
                opti.subject_to(U(1, k) <= deg2rad(30));

                cost = cost + ((X([2,4], k)-obj.ref)' * Q * (X([2,4], k)-obj.ref));
                cost = cost + ((U(:,k)-U(:,k-1))' * R * (U(:,k)-U(:,k-1)));
                cost = cost + (10000*slack(k)^2);
                
            end
            

            opti.minimize(cost);

            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Store the defined problem to solve in get_u
            obj.opti = opti;

            % Setup solver
            options = struct;
            options.ipopt.print_level = 0;
            options.print_time = 0;
            options.expand = true;
            obj.opti.solver('ipopt', options);
        end

        function u = get_u(obj, x0, ref, x0other)

            if nargin < 4
                x0other = zeros(4, 1);
            end

            % Compute solution from x0
            obj.solve(x0(1:4), ref, x0other(1:4));

            u = obj.sol.value(obj.u0);
        end

        function solve(obj, x0, ref, x0other)

            % Pass parameter values
            obj.opti.set_value(obj.x0, x0);
            obj.opti.set_value(obj.ref, ref);
            obj.opti.set_value(obj.x0other, x0other);

            obj.sol = obj.opti.solve();   % actual solve
            
            % Set warm start for next solve
            obj.opti.set_initial(obj.sol.value_variables());
            obj.opti.set_initial(obj.opti.lam_g, obj.sol.value(obj.opti.lam_g));
        end
    end
end
