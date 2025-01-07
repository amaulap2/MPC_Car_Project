classdef LonEstimator
    properties
        % continous sub-system
        sys
        % Extended linearization points
        xs_hat, us_hat
        % Extended system matrices
        A_hat, B_hat, C_hat
        % Observer gain matrix
        L
    end
    
    methods
        % Setup the longitudinal disturbance estimator
        function est = LonEstimator(sys, Ts)

            xs = sys.UserData.xs;
            us = sys.UserData.us;
            
            % Discretize the system and extract the A,B,C,D matrices
            [~, Ad, Bd, Cd, ~] = Car.c2d_with_offset(sys, Ts);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

            % Extend the system to include disturbance dynamics
            nx = size(Ad, 1);  % Number of states
            ny = size(Cd, 1);  % Number of outputs

            Bd_hat = Bd; % As Bd (2x1)

            Ob = obsv(Ad, Cd);



            % Extended state dynamics
            est.xs_hat = [xs; 0];  % Extended state including disturbance
            est.us_hat = us;
            est.A_hat = [Ad, Bd_hat; zeros(1, nx), 1];
            est.B_hat = [Bd; 0];
            est.C_hat = [Cd, zeros(ny, 1)];
            
            % Place poles slightly faster than the system dynamics
            poles = eig(est.A_hat)*0.8;  % Pole placement
            est.L = place(est.A_hat', est.C_hat', poles)';
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        % This function takes in the the estimate, input, and measurement
        % at timestep i and predicts the next (extended) state estimate at
        % the next timestep i + 1.
        function z_hat_next = estimate(est, z_hat, u, y)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   z_hat      - current estimate (V, dist)
            %   u          - longitudinal input (u_T)
            %   y          - longitudinal measurement (V)
            % OUTPUTS
            %   z_hat_next - next time step estimate
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

            A_hat_reduced = est.A_hat(2:end,2:end);
            B_hat_reduced = est.B_hat(2:end);
            C_hat_reduced = est.C_hat(2,2:end);
            L_hat_reduced = est.L(2:end,2);
            xs_hat_reduced = est.xs_hat(2:end);
            us_hat_reduced = est.us_hat;

            % Estimation equation
            z_hat_next = xs_hat_reduced + A_hat_reduced * (z_hat-xs_hat_reduced) + B_hat_reduced * (u - us_hat_reduced) + L_hat_reduced * (y - C_hat_reduced * z_hat);

            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end
