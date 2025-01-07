%% Check observability
clc;clear;
Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us] = car.steady_state(80 / 3.6);
sys = car.linearize(xs, us);
sysd = c2d(sys, Ts);
[Ad, Bd, Cd, ~] = ssdata(sysd);

Q = 10*eye(2);
R = 2;

K = 

Ad_lon = Ad([1,4],[1,4]);
Bd_lon = Bd([1,4],2);
Cd_lon = Cd([1,4],[1,4]);

max_iter = 100;
tol = 1e-3;

Eps = tube_mpc_sets(Ad_lon, Bd_lon, us(2), max_iter, tol);


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

X(4,end-20:end)