%%% Problem 1
%%% First model the problem, and then solve it using ode45.
input = 2; 
time = 0:0.01:6; 
dt = 0.01; 
dvdt1 = @(t, v) (input - v);
[t, V1] = ode45(dvdt1, time, 1); 
A1 = V1; 
dvdt2 = @(t, v) (input - exp(-t)); 
[t, V2] = ode45(dvdt2, time, 1); 
A2 = V2;
index = find(A2 > 10, 1); 
if isempty(index)
    A3 = max(time); 
else
    A3 = (index - 1) * dt; 
end
A4 = 2;
S1 = @(t, y) (input - y); 
[t, y] = ode45(S1, time, input); 
A5 = y; 
increase = (input * (time - time(1))) + 2; 
A6 = increase'; 
A7 = A6(index) / A2(index); 

%%% Problem 2
%%% Use finite differences for boundary value problems and loop to iterate
%%% each timstep
C = 1;
dx = 0.01; 
dt = 0.01;
x = linspace(-1, 1, 2/dx + 1);
t = 0:0.01:1;
PIC = @(x) exp(1 - 1./(1 - x.^2)); 
PBC = [0, 0];

P = PIC(x);
P(1) = PBC(1);
P(end) = PBC(2);
N = length(x) - 2;
l= C * dt / dx^2;
A = spdiags([ones(N,1)*(-l/2), ones(N,1)*(1+l), ones(N,1)*(-l/2)], -1:1, N, N);
b = zeros(N, 1);
A8 = full(A);
A9 = zeros(N, 1);
A10 = zeros(N, 1);
for k = 2:length(t)
    b(1) = l/2 * (P(1) + P(3)) + (1 - l) * P(2);
    for i = 2:N-1
        b(i) = l/2 * (P(i+2) + P(i)) + (1 - l) * P(i+1);
    end
    b(N) = l/2 * (P(N) + P(N+2)) + (1 - l) * P(N+1);
    b(1) = b(1) + l/2 * PBC(1);
    b(N) = b(N) + l/2 * PBC(2);
    
    if k == 2
        A9 = b;
    end
    P(2:end-1) = A\b;
    if k == length(t)
        A10 = b;
    end
end
A11 = P(x == 0.5);
A12 = 0;

%%% Problem 3

load CP10_M1.mat
load CP10_M2.mat
load CP10_M3.mat
%%%%%%%%%%%%%%%%%
Matrix1 = double(M1);
Matrix2 = double(M2);
Matrix3 = double(M3);
Data = Matrix2 * Matrix3 * Matrix1';
A13 = numel(Data) * 8 / 1e6;
[m, n] = size(Data);

k=0;
while [(m*k+k+n*k)/(m*n)]<=0.99
    k=k+1;
end
A14 = k;
data = Data / max(Data(:));
imshow(data, 'Colormap', gray, 'DisplayRange', [0, 1]);
A15 = (m * k + k + n * k) * 8 / 1e6;
A16 = 17;

