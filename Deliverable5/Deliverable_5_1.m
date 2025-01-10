%% Check observability
clc;clear;
Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us] = car.steady_state(80 / 3.6);
sys = car.linearize(xs, us);
sysd = c2d(sys, Ts);
[Ad, Bd, Cd, ~] = ssdata(sysd);

Q = 10*eye(2);
R = 30;

save('Q.mat','Q')
save('R.mat','R')

Ad_lon = Ad([1,4],[1,4]);
Bd_lon = Bd([1,4],2);
Cd_lon = Cd([1,4],[1,4]);

max_iter = 100;
tol = 1e-2;

Eps = tube_mpc_sets(Ad_lon, Bd_lon, Q, R, us(2), max_iter, tol);
%%
load("K.mat")
load("Eps.mat")
load('U_tight.mat')
load('X_tight.mat')
load("Xf.mat")
plot(U_tight,'Color','r')
hold on;
% plot(Xf,'Color','b')

%%
clc;clear;
Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us] = car.steady_state(80 / 3.6);
sys = car.linearize(xs, us);
[sys_lon, sys_lat] = car.decompose(sys);

H_lon = 15; % Horizon length in seconds

ref = [0 120/3.6]';
mpc_lon = MpcControl_lon(sys_lon, Ts, H_lon);
mpc_lat = MpcControl_lat(sys_lat, Ts, H_lon);
mpc = car.merge_lin_controllers(mpc_lon, mpc_lat);


otherRef = 100 / 3.6;

params = {};
params.Tf = 25;
params.myCar.model = car;
params.myCar.x0 = [0 0 0 80/3.6]';
params.myCar.u = @mpc.get_u;
params.myCar.ref = ref;
params.otherCar.model = car;
params.otherCar.x0 = [15 0 0 otherRef]';
params.otherCar.u = car.u_const(otherRef);
result = simulate(params);
visualization(car, result);

%% second test 
clc;clear;
Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us] = car.steady_state(80 / 3.6);
sys = car.linearize(xs, us);
[sys_lon, sys_lat] = car.decompose(sys);

H_lon = 25; % Horizon length in seconds

ref = [0 120/3.6]';
mpc_lon = MpcControl_lon(sys_lon, Ts, H_lon);
mpc_lat = MpcControl_lat(sys_lat, Ts, H_lon);
mpc = car.merge_lin_controllers(mpc_lon, mpc_lat);

params = {}; 
params.Tf = 25; 
params.myCar.model = car; 
params.myCar.x0 = [0 0 0 115/3.6]'; 
params.myCar.u = @mpc.get_u;
params.myCar.ref = ref;
params.otherCar.model = car; 
params.otherCar.x0 = [8 0 0 120/3.6]'; 
params.otherCar.u = car.u_fwd_ref(); 
params.otherCar.ref = car.ref_robust();
result = simulate(params);
visualization(car, result);