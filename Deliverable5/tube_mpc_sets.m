function Eps = tube_mpc_sets(Ad, Bd, K, max_iter, tol)


    
    
    % Define the 1D disturbance set W for u_T ± 0.5
    W = polytope([1 0;0 1; -1 0; 0 -1],[Bd*0.5; Bd*0.5]);  % 1D polytope for disturbance

    
    % Initialize the invariant set with the origin in 2D space
    Eps = polytope(zeros(2, 1), zeros(2, 1));  % {0} in 2D
    
    % Iterative computation to grow the invariant set
    for i = 1:max_iter
        % Propagate the disturbance through the system (A_d + B_d K)^i * W
        preEps = (Ad + Bd * K)^i * W;
        
        % Compute Minkowski sum to accumulate disturbances
        Eps = Eps + preEps;
        
        % Reduce polytope complexity by minimizing its H-representation
        Eps.minHRep();
        
        % Terminate if the norm of (A + BK)^i is sufficiently small
        if norm((Ad + Bd * K)^i, 'fro') < tol
            fprintf('Converged after %d iterations.\n', i);
            break;
        end
    end
    
    % Plot the resulting minimal robust invariant set
    figure;
    plot(Eps);
    title('Minimal Robust Invariant Set (mRPI)');
    xlabel('State Dimension 1');
    ylabel('State Dimension 2');
end
