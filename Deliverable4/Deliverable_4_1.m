%% Check observability
clc;clear;
Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us] = car.steady_state(80 / 3.6);
sys = car.linearize(xs, us);
[~, Ad, Bd, Cd, ~] = Car.c2d_with_offset(sys, Ts);

Ob = obsv(Ad,Cd);

% Check rank
n = size(Ad,1);
rank_Ob = rank(Ob);

if rank_Ob == n
    disp('System is observable.');
else
    disp(['System is NOT observable. Rank: ', num2str(rank_Ob)]);
end


%%
clc;clear;
Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us] = car.steady_state(80 / 3.6);
sys = car.linearize(xs, us);
[sys_lon, sys_lat] = car.decompose(sys);

H_lon = 15; % Horizon length in seconds

mpc_lon = MpcControl_lon(sys_lon, Ts, H_lon);
mpc_lat = MpcControl_lat(sys_lat, Ts, H_lon);
mpc = car.merge_lin_controllers(mpc_lon, mpc_lat);

estimator = LonEstimator(sys_lon, Ts);

x0 = [0 0 0 80/3.6]'; % (x, y, theta, V) 
ref1 = [0 80/3.6]'; % (y ref, V ref) 
ref2 = [3 50/3.6]'; % (y ref, V ref)
params = {}; 
params.Tf = 15; 
params.myCar.model = car; 
params.myCar.x0 = x0;
params.myCar.est_fcn = @estimator.estimate;
params.myCar.est_dist0 = 0;
params.myCar.u = @mpc.get_u; 
params.myCar.ref = car.ref_step(ref1, ref2, 2); % delay reference step by 2s;
result = simulate(params);
visualization(car, result);

%% Plotting all closed loop dynamics
X = result.myCar.X;
U = result.myCar.U;

figure;

plot(0:(length(X)-1), X(4,:), '-o', 'LineWidth', 2);
xlabel('Time t (deciseconds)');
ylabel('Velocity V (m/s)');
yline(ref2(2),'--','LineWidth', 2);
title('Velocity');
legend('V','Reference','Location','southeast');
grid on;
