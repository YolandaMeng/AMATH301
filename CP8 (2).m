%%% Problem 1
%%% Initialize t, and x_true

t = 0:0.1:10;
x_true = @(t) 0.5*(cos(t)+sin(t)+exp(-t));
x_sol = x_true(t);
x0 = 1; 
t0 = 0; 
t_end = 10;
h = 0.1; 
n = (t_end - t0)/ h;
x = zeros(1,n);
t = zeros(1,n);
t(1) = t0;
x(1) = x0;
f = @(t,x) cos(t)-x;

for k = 1:n
  x(k+1) = x(k) + h * f(t(k),x(k));
  t(k+1) = t(k) + h;
  A1 = x;  
end

A2 = abs(A1 - x_sol); 

for k = 1:n
    t(k+1)= t(k) + h;
    x(k+1)=fzero(@(x) x(k) + h * (cos(t(k+1)) - x) - x,x(k));
    A3 = x;
end

A4 = abs(A3-x_sol);


funDiff = @(t,x) cos(t)-x;
X0 = 1;
t_span = 0:0.1:10;
[t,x_sol] = ode45(@(t,x) funDiff(t,x),t_span,X0);
A5 = x_sol';
A6 = abs(A5-x_sol);


%%% Problem 2
%%% Initialize the parameters

a = (8.1)^2;
t_span = 0:0.01:2; 
x_true = @(t) 2 * atan(exp(a*t) / (1 + sqrt(2)));
dt = 0.01;
t = 0:dt:2;
xf = zeros(1, length(t));
xf(1) = pi / 4;

% Forward Euler Scheme

for k = 2:length(t)
    xf(k) = xf(k-1) + a * dt * sin(xf(k-1));
end

A8 = max(abs(xf - x_true(t)));
A7 = xf;

dt = 0.001;
t = 0:dt:2;
xf1 = zeros(1, length(t));
xf1(1) = pi / 4;

% Forward Euler Scheme

for k = 2:length(t)
    xf1(k) = xf1(k-1) + a * dt * sin(xf1(k-1));
end

A9 = A8 / max(abs(xf1 - x_true(t)));

%% Backward Euler Method with Predictor-Corrector
dt = 0.01;
t = 0:dt:2;
xb = zeros(1, length(t));
xb(1) = pi / 4;

for k = 2:length(t)
    p = xb(k-1) + a * dt * sin(xb(k-1));
    xb(k) = xb(k-1) + a * dt * sin(p);
end

A11 = max(abs(xb - x_true(t)));
A10 = xb;

dt = 0.001;
t = 0:dt:2;
xb1 = zeros(1, length(t));
xb1(1) = pi / 4;

for k = 2:length(t)
    p1 = xb1(k-1) + a * dt * sin(xb1(k-1));
    xb1(k) = xb1(k-1) + a * dt * sin(p1);
end

A12 = A11 / max(abs(xb1 - x_true(t)))

%% Builtin Solver ode45
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_ode45, x_ode45] = ode45(@(t, x) a * sin(x), t_span, pi / 4);


A14 = max(abs(x_ode45 - x_true(t_ode45)));
A13 = x_ode45;

ts = 0:0.001:2;
[t_ode45_new, x_ode45_new] = ode45(@(t, x) a * sin(x), ts, pi / 4);

A15 = A14 / max(abs(x_ode45_new - x_true(t_ode45_new)));

%%%  If you want to write local functions, put them here
