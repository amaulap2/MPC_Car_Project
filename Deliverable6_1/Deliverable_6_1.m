clc;clear;
Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us] = car.steady_state(80 / 3.6);

H = 15; mpc = NmpcControl(car, H);
x0 = [0 0 0 80/3.6]'; 
ref = [3 100/3.6]';
u = mpc.get_u(x0, ref); % check if the open−loop prediction is reasonable

params = {}; 
params.Tf = 10; 
params.myCar.model = car; 
params.myCar.x0 = [0 0 0 80/3.6]'; 
params.myCar.u = @mpc.get_u;
ref1 = [0 80/3.6]'; 
ref2 = [3 100/3.6]';
params.myCar.ref = car.ref_step(ref1, ref2, 2);
result = simulate(params); 
visualization(car, result);