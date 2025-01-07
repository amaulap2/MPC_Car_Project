Ts = 1/10;
car = Car(Ts);
u = [0.1, 0.5]';       % Slight steering left(+) right(-) , moderate throttle accel(+) decel(-)
x = [0, 0, 0.1, 20]';  % Start at origin, slight angle, 20 m/s speed
x_dot = car.f(x, u);
disp(x_dot);

%% Multiple steps for 2s simulation
car = Car(Ts);
Tf = 2.0; % Simulation end time
x0 = [0, 0, deg2rad(-2), 20/3.6]'; % (x, y, theta, V) Initial state
u = [deg2rad(-1), 0.7]'; % (delta, u T) Constant input
params = {}; % Setup simulation parameter struct
params.Tf = Tf;
params.myCar.model = car;
params.myCar.x0 = x0;
params.myCar.u = u;
result = simulate(params); % Simulate nonlinear model
visualization(car, result);

%% linearization around xs,us
car = Car(Ts);
Vs = 120/3.6; % 120 km/h
[xs, us] = car.steady_state(Vs); % Compute steady−state for which f s(xs,us) = 0
sys = car.linearize(xs, us); % Linearize the nonlinear model around xs, us

%% decomposition
[sys_lon, sys_lat] = car.decompose(sys);
Ts = 1/10;
[fd_xs_us, Ad, Bd, Cd, Dd] = Car.c2d_with_offset(sys, Ts);
delta_X = x - xs;
delta_U = u - us;
