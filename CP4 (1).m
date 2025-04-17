
%%%  Problem 1
%%%  First go to the end of your m-file to create your Jacobi and
%%%  Gauss-Seidel functions.

%%% Once you have created your functions return here.
%%% Initialize your matrix A and RHS b
A = [1.1 0.2 -0.2 0.5; 0.2 0.9 0.5 0.3; 0.1 0 1 0.4; 0.1 0.1 0.1 1.2];
b = [1;0;1;0];


%%% Use your Jacobi and Gauss-Seidel functions to find A1 through A4.
epsilons = [10^-2, 10^-4, 10^-6, 10^-8];
Tj = zeros(4, 1);
Ej = zeros(4, 1);
Tgs = zeros(4, 1);
Egs = zeros(4, 1);
  for i=1:4
[Tj(i), Ej(i)] = jacobi(A, b, epsilons(i));
[Tgs(i), Egs(i)] = gauss_seidel(A, b, epsilons(i));
  end
A1 = Tj';
A2 = Ej';
A3 = Tgs';
A4 = Egs';




%%%  Problem 2
%%%  Initialize your Day 0 vector x
x0 = [0.9; 0.09; 0.01];


      
%%%  Part 1: without a vaccine
%%%  Make sure to have p = 0
%%%  Initialize the SIR matrix M, and save it as A5
p = 0;
recover_rate = 1/1000;
infected_rate = 1/200;
mutation_rate = 1/10000;
M = [1-infected_rate 0 mutation_rate; infected_rate 1-recover_rate 0; 0 recover_rate 1-mutation_rate];
A5 = M;



%%%  Create a loop to find the day that the number of infected
%%%  individuals hits 50% and another loop for the steady state of the
%%%  infected population
%%%  There is a way to put everything under one loop if you make clever use
%%%  of conditionals
day = 0;
while x0(2)<=0.5
    x0 = M*x0;
    day = day+1;
end
D0 = day;

for i=1:100000
    x1=x0;
    x0 = M*x0;
    if abs(x0(2)- x1(2))< 1e-8
        break;
    end
end
F0 = x0(2);


%%% Save the days and steady state in a row vector A6

A6 = [D0, F0];



%%%  Reinitialize your Day 0 vector x

x0 = [0.9; 0.09; 0.01];


%%%  Part 2: with a vaccine
%%%  Make sure to have p = 2/1000
%%%  Initialize the SIR matrix M, and save it as A7
p = 2/1000;
M = [1-infected_rate-p 0 mutation_rate; infected_rate 1-recover_rate 0; p recover_rate 1-mutation_rate];
A7 = M;



%%%  Create a loop to find the day that the number of infected
%%%  individuals hits 50% and another loop for the steady state of the
%%%  infected population
%%%  There is a way to put everything under one loop if you make clever use
%%%  of conditionals
day = 0;
while x0(2)<=0.5
    x0 = M*x0;
    day = day+1;
end
D1 = day;

for i=1:100000
    x1=x0;
    x0 = M*x0;
    if abs(x0(2)- x1(2)) < 1e-8
        break;
    end
end
F1 = x0(2);



%%% Save the days and steady state in a row vector A8

A8 = [D1, F1];
 
 
%%%  Problem 3
  
%%%  Initialize your 114x114 tridiagonal matrix A
n = 114;
A9 = diag(2*ones(n,1)) - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1);



%%%  Initialize your 114x1 RHS column vector rho
A10 = zeros(114,1);
for j_values = 1:n
    rho = 2 * (1 - cos((53*pi)/115)) * sin((53*pi*j_values)/115);
    A10(j_values) = rho;
end


%%%  Implement Jacobi's method for this system.
%%%  Don't use the function you created before because that was designed for
%%%  a specific task, and will not work here.
phi = ones(114,1);
tolerance = 1e-5;
D = diag(A9);
L = tril(A9, -1);
U = triu(A9, 1);
M = -(L+U)./D;
c = A10./D;

for i = 2:100000
    phi_previous = phi;
    phi = M * phi + c;
    if max(abs(phi - phi_previous)) < tolerance
        A11 = phi;
        A12 = i;
        break;
    end
end



%%%  Create a column vector phi that contains the exact solution given in
%%%  the assignment file
for j = 1:114
    phi(j) = sin(53*pi*j/115);
end

true_solution = phi;


%%%  Save the difference of the Jacobi solution and the exact solution as
%%%  A13.  Use the maximal entry in absolute value to calculate this error.
A13 = max(abs(A11 - true_solution));


%%%  Implement Gauss-Seidel for this system.
%%%  Don't use the function you created before because that was designed for
%%%  a specific task, and will not work here.
phi = ones(114, 1);
tolerance = 1e-5;
LpD = tril(A9);
U = triu(A9,+1);
M = -LpD\U;
c = LpD\A10;

for i= 2:100000
    phi_previous = phi;
    phi = M * phi + c;
    if max(abs(phi - phi_previous)) < tolerance
        A14 = phi;
        A15 = i;
        break;
    end
end



%%%  Save the difference of the Gauss-Seidel solution and the exact solution as
%%%  A13.  Use the maximal entry in absolute value to calculate this error.
A16 = max(abs(A14 - true_solution));



  

%%% Jacobi and Gauss Seidel Iteration functions
%%% Create your functions here
%%% Both functions will need two outputs and three inputs
%%% The code within the function will be very similar to
%%% Week 4 coding lecture 2

%% function of jacobi method.
function [T, E] = jacobi(A, b, eps)
n = length(A);
y = zeros(n, 1);
T = 0;
while true
T = T + 1;
y_new = y;
for i = 1:n
    y_new(i) = (b(i) - A(i, [1:i-1,i+1:n])*y([1:i-1,i+1:n])) / A(i, i);
end
    E = max(abs(A*y_new - b));
if E < eps
break;
end
y = y_new;
end
end
%% function of gauss_seidel method.
function [T, E] = gauss_seidel(A, b, eps)
n = length(A);
y = zeros(n, 1);
T = 0;
while true
T = T + 1;
 for i = 1:n
y(i) = (b(i) - A(i, [1:i-1,i+1:n])*y([1:i-1,i+1:n])) / A(i, i);
 end
E = max(abs(A*y - b));
if E < eps
        break;
       end
    end
end