%%
clc;clear;
Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us] = car.steady_state(120 / 3.6);
sys = car.linearize(xs, us);
[sys_lon, sys_lat] = car.decompose(sys);

H_lon = 15; % Horizon length in seconds

mpc_lon = MpcControl_lon(sys_lon, Ts, H_lon);
mpc_lat = MpcControl_lat(sys_lat, Ts, H_lon);
mpc = car.merge_lin_controllers(mpc_lon, mpc_lat);

x0 = [0 0 0 80/3.6]'; % (x, y, theta, V)
ref1 = [0 80/3.6]'; % (y ref, V ref)
ref2 = [3 120/3.6]'; % (y ref, V ref)
params = {};
params.Tf = 15;
params.myCar.model = car;
params.myCar.x0 = x0;
params.myCar.u = @mpc.get_u;
params.myCar.ref = car.ref_step(ref1, ref2, 5); % delay reference step by 5s, makes following ref 1 for the 5 first seconds then wsitch to ref 2
result = simulate(params);
visualization(car, result);

%% Plotting all closed loop dynamics
X = result.myCar.X;
U = result.myCar.U;

figure;
subplot(2,2,1);
plot(0:(length(X)-1), X(1,:), '-o', 'LineWidth', 2);
xlabel('Time t (deciseconds)');
ylabel('Position x (m)');
title('Longitudinal Position x');
legend('x');
grid on;

subplot(2,2,2);
plot(0:(length(X)-1), X(2,:), '-o', 'LineWidth', 2);
xlabel('Time t (deciseconds)');
ylabel('Position y (m)');
yline(ref2(1),'--','LineWidth', 2);
title('Lateral Position');
legend('y', 'Reference', 'Location','southeast');
grid on;

subplot(2,2,3);
plot(0:(length(X)-1), X(3,:), '-o', 'LineWidth', 2);
xlabel('Time t (deciseconds)');
ylabel('Angle \theta (rad)');
title('Heading Angle');
legend('\theta');
grid on;

subplot(2,2,4);
plot(0:(length(X)-1), X(4,:), '-o', 'LineWidth', 2);
xlabel('Time t (deciseconds)');
ylabel('Velocity V (m/s)');
yline(ref2(2),'--','LineWidth', 2);
title('Velocity');
legend('V','Reference','Location','southeast');
grid on;

%% plot invariant set lateral
clc;clear;
Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us] = car.steady_state(120 / 3.6);
sys = car.linearize(xs, us);
[sys_lon, sys_lat] = car.decompose(sys);

A = sys_lat.A;
B = sys_lat.B;
nx = size(B,1);
nu = size(B,2);
Q = 10*eye(nx); % Weight for state tracking error
R = 2*eye(nu);

% Constraints
M = [1;-1]; m = [1; 1];
F = [1 0; 0 1; -1 0; 0 -1]; f = [3.5; deg2rad(5); 0.5; deg2rad(5)];

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

%%
figure;
Xf.plot();
xlabel('x_1');
ylabel('x_2');
yline(deg2rad(-5), '-.', 'LineWidth',2);
yline(deg2rad(5), '-.', 'LineWidth',2);
xline(-0.5, '-.', 'LineWidth',2, 'Color', 'Blue');
xline(3.5, '-.', 'LineWidth',2, 'Color', 'Blue');
title('Terminal Set X_f');
axis([-1,4,-0.1,0.1])
legend('X_f','theta\_bounds','','y\_bounds')
grid on;



