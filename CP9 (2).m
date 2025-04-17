%%% Problem 1
%%% Use ode45 to solve the Fitzhugh-Nagumo IVP
%%% For the maxima use the plot to narrow down the times you use to search
%%% for the maximum.

I = @(t)(5+sin(pi*t/10))/10;
a = 0.7;
b = 1;
tau = 12;
ode = @(t,y)[y(1)-y(1)^3/3-y(2)+I(t);(a+y(1)-b*y(2))/tau];
y0 = [1;0];
t = 0:0.5:100;
[t,y] = ode45(ode,t,y0);
v = y(:,1);
A1 = v;
%plot solution
plot(t,v)
xlabel('t')
ylabel('v')
title('FitzHugh-Nagumo')
peaks = v(1:end-2)<v(2:end-1)& v(3:end)<v(2:end-1);
peaks = logical([0;peaks;0]);
t_peaks = t(peaks);
A2 = t_peaks(1);
A3 = t_peaks(3);
A4 = 1/(A3-A2);


%%% Problem 2
%%% Use ode45 to solve the Chua equation
%%% You can tell something is chaotic if it is seemingly random
%%% If it looks like all solutions tend toward a point or a circle it is
%%% not chaotic.

alpha = 16;
beta = 30;

% Set the initial conditions
initial_conditions = [0.1; 0.2; 0.3];

% Set the time span
t_span = [0 100];

% Set the time step
dt = 0.05;

% Solve the system of ODEs using ode45
[t, y] = ode45(@(t, y) myODEs(t, y, alpha, beta), t_span, initial_conditions);

% Plot the phase space
figure;
plot3(y(:, 1), y(:, 2), y(:, 3));
xlabel('x');
ylabel('y');
zlabel('z');
title('Phase Space Plot');

% Check for chaos and save the result
A5 = 1; % Assume chaotic, update if non-chaotic conditions are met
% Add logic to check for non-chaotic conditions and update A5 if necessary

% Save A5 and A6
save('A5.mat', 'A5');

% Set the new initial conditions
new_initial_conditions = [0.1; 0.2 + 1e-5; 0.3];

% Solve the system of ODEs with the new initial conditions using ode45
[t_new, y_new] = ode45(@(t, y) myODEs(t, y, alpha, beta), t_span, new_initial_conditions);

% Interpolate the solutions to a common set of time points
A6 = interp1(t, y, t_new, 'linear');

% Calculate the maximum difference
A7 = max(max(abs(A6 - y_new)));

% Display the maximum difference
disp(['Maximum difference in absolute value between A6 and the new solution (A7): ', num2str(A7)]);

% Save A7
save('A7.mat', 'A7');

% Set the new value for beta
beta_new = 100;

% Set the initial conditions for the new scenario
initial_conditions_new = [0.1; 0.2; 0.3];

% Solve the system of ODEs with the new parameters using ode45
[t_new, y_new] = ode45(@(t, y) myODEs(t, y, alpha, beta_new), t_span, initial_conditions_new);

% Plot the phase space for the new scenario
figure;
plot3(y_new(:, 1), y_new(:, 2), y_new(:, 3));
xlabel('x');
ylabel('y');
zlabel('z');
title('Phase Space Plot (New Scenario)');

% Check for chaos and save the result
A8 = 0; % Assume chaotic, update if non-chaotic conditions are met
% Add logic to check for non-chaotic conditions and update A8 if necessary

% Save A8 and the solution for the new scenario
save('A8.mat', 'A8');
save('A9.mat', 'y_new');



%%% Problem 3

%%% Part 1: Finite Differences
%%% Use finite differences to solve the BVP
%%% Be careful about the shape of the vectors, you may have to transpose to
%%% get the correct shape.  It's a good idea to print the solutions out to
%%% make sure the shape is correct.

%% 
t0 = 0;  
tE = 6;   
N = 61;   
dt = (tE - t0) / (N - 1); 
t = linspace(t0, tE, N);  
A = zeros(N-2, N-2); 
b = zeros(N-2, 1);

for i = 1:(N-2)
    A(i, i) = -2/dt^2 + 1;
    if i > 1
        A(i, i-1) = 1/dt^2;
    end
    if i < N-2
        A(i, i+1) = 1/dt^2;
    end
    b(i) = 5 * cos(4 * t(i+1));
end

b(1) = b(1) - 1/dt^2; 
b(N-2) = b(N-2) - 0.5/dt^2; 
A9 = A; 
A10 = b; 
x_interior = A\b;
x_full = [1; x_interior; 0.5];
A11 = x_full;
true_solution =zeros(61,1);
C1 = ((1/2) + (1/3) * cos(24) - (4/3) * cos(6)) / sin(6);
C2 = 4/3;
for i = 0:60
    true_solution(i+1) = C1 * sin(i/10) + C2 * cos(i/10) - (1/3) * cos(4 * i/10);
end

A12 = max(abs(true_solution - A11));





%%
%%% Part 2: Shooting Method via Bisection
%%% Use the shooting method to solve the BVP
%%% It's a good idea to test out a few in the command window first to make
%%% sure that your initial conditions gets you to different sides of the right
%%% boundary condition.
%%% Use the plot to help you figure out what your choices of initial
%%% conditions should be
x0 = 1;
v1 = 1;
v2 = 3;
xT = 0.5
v_mid = (v1 + v2)/2;
t = [0:0.1:6];
[T, Y] = ode45(@(t, x) myODE(t,x), [0:0.1:6], [x0 v1]);
t_a = T; % for plotting
x_a = Y(:, 1);
[T, Y] = ode45(@(t, x) myODE(t,x), [0:0.1:6], [x0 v2]);
t_b = T; % for plotting
x_b = Y(:, 1);
[T, Y] = ode45(@(t, x) myODE(t,x), [0:0.1:6], [x0 v_mid]);
t_mid = T; % for plotting
x_mid = Y(:, 1);
plot(t_a, x_a, t_mid, x_mid, '--', t_b, x_b, 'linewidth', 4)

for i = 1:100
    if abs(x_mid(end)-0.5) < 1e-8
        break
    elseif sign(x_mid(end) - 0.5) ~= sign(x_a(end)-0.5)
        v2 = v_mid;
        [T, Y] = ode45(@(t, x) myODE(t,x), [0:0.1:6], [x0 v2]);
        x_b = Y(:, 1);
    else
        v1 = v_mid;
        [T, Y] = ode45(@(t, x) myODE(t,x), [0:0.1:6], [x0 v1]);
        x_a = Y(:, 1);
end

    v_mid = (v1 + v2)/2;
[T, Y] = ode45(@(t, x) myODE(t,x), [0:0.1:6], [x0 v_mid]);
x_mid = Y(:, 1);
end

A13 = x_mid;
A14 = max(max(abs(A13-xT)));
A15 = max(abs(A13-A11));

%%% You can set up your ODEs as functions here if you like
function dydt = myODEs(t, y, alpha, beta)
    dydt = [
        alpha * (y(2) + (1/6) * y(1) - (1/16) * y(1)^3);
        y(1) - y(2) + y(3);
        -beta * y(2);
    ];
end

function dydt = odesys(t, y)
    % y(1) corresponds to x
    % y(2) corresponds to dx/dt
    dydt = [y(2); 5*cos(4*t) - y(1)];
end

function dx = myODE(t,x)
    dx1 = x(2);
    dx2 = 5*cos(4*t)-x(1);
    dx = [dx1; dx2];
end
