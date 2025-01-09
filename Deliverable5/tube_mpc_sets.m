function Eps = tube_mpc_sets(Ad, Bd, Q, R, uT_s, max_iter, tol)

% Design LQR Controller 
% Method 1
[K,Qf,~] = dlqr(Ad, -Bd, Q, R);
K = -K;
save('K.mat', 'K');
save('Qf.mat','Qf');

% Define the 1D disturbance set W and map through Bd
W = Polyhedron([-1; 1], [0.5; 0.5]);  % u_T ∈ [-0.5, 0.5]
W = Bd * W;  % Map to 2D via Bd

% Initialize the invariant set at the origin in 2D space
Eps = W;


% Iterative computation to grow the invariant set
for i = 1:max_iter
    % Initialize propagated disturbance
    preEps = Eps;
    Eps = preEps + (Ad - Bd * K)^i * W;
    
    % Reduce polytope complexity by minimizing its H-representation
    Eps.minHRep();
    disp(i);
    % Terminate if the norm of (A + BK)^i is sufficiently small
    if norm((Ad - Bd * K)^i) < tol
        fprintf('Converged after %d iterations.\n', i);
        break;
    end
end

% Plot the minimal invariant set
figure;
plot(Eps);
save('Eps.mat', 'Eps');
title('Minimal Robust Invariant Set (mRPI)');
xlabel('State (x)');
ylabel('State (V)');

x_safe = 7;
save('x_safe',"x_safe")

X = Polyhedron([-1 0],-(6-x_safe));  % the space between the two cars should be at least 6 which is 1.7m before the two cars touche.
U = Polyhedron([1; -1], [1; 1]);          % Input constraints

Eps.minHRep();

KE = K * Eps;
KE.minHRep();

X_tight = X - Eps;   % X ⊖ E
U_tight = U - KE;    % U ⊖ (K * E)

save('X_tight.mat','X_tight');
save('U_tight.mat','U_tight');

% Plot the maximal terminal set
F = X_tight.A;
f = X_tight.b;
M = U_tight.A;
m = U_tight.b;

Xf = Polyhedron([F;M*K],[f;m]);
Acl = [Ad-Bd*K];
while 1
    prevXf = Xf;
    T = Xf.A;
    t = Xf.b;
    preXf = Xf*Acl;
    Xf = intersect(Xf, preXf);
    plot(Xf)
    if prevXf.volume - Xf.volume < 1e-4
        break
    end
end

plot(Xf);
save('Xf.mat','Xf');






end


