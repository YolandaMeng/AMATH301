%%% Problem 1
data = readmatrix('population.csv');
t = data(1, :);
N = data(2, :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Determine your stepsize dt from the vector t
h = t(2) - t(1);
t_2020 = find(t == 120);  % Assuming 120 represents the year 2020
A1 = (N(t_2020 + 1) - N(t_2020 - 1)) / (2 * h);
t_1880 = find(t == -20);  % Assuming -20 represents the year 1880
A2 = (N(t_1880 + 1) - N(t_1880 - 1)) / (2 * h);
t_1790 = find(t == -110);  % Assuming -110 represents the year 1790
A3 = (N(t_1790 + 1) - N(t_1790)) / h;
A4 = zeros(1, 24);
A4(1) = A3;
A4(24) = A1;

for i = 2:23
    A4(i) = (N(i + 1) - N(i - 1)) / (2 * h);
end

A5 = A4./ N;
A6 = mean(A5);




%%% For dN/dt you will need to use a combination of the above differences,
%%% but the choice will be obvious based on which direction you can/cannot
%%% go in the horizontal axis.  Whenever possible use central difference;
%%% only use forward or backward when central is not possible.




%%% Problem 2
data = readmatrix('brake_pad.csv');
r = data(1, :);
T = data(2, :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Determine your stepsize dr from the vector r
d_r = 0.017;


%%% Use the LHR formula from the coding lecture

b = 0.7051;
a = sum(T(1:end -1).* r(1:end -1).*b)* d_r;
A7 = a;
A8 = A7/ sum(r(1:end -1).*b)* d_r;



%%% Use the RHR formula from the coding lecture

A9 = sum(T(2:end).* r(2:end).*b)* d_r;
A10 = A9/ sum(r(2:end).*b)* d_r;




%%% Use the Trapezoid rule formula or the trapz function from the coding lecture

A11 = (A7 + A9) /2;
A12 = A11 / (d_r)/2 *(r(1) + 2* sum(r(2:end -1)+ r(end) *b ));




%%% Problem 3
%%% You'll have to use anonymous functions here.  You can see the syntax in
%%% the Numerical Integration coding lecture where the builtin function
%%% "integrand" is used.

f = @(x)(x^2/2 - x^3/ 3);
m = @(x,z) ((x^2* z^2)/2 - (x^3 * z^3/3));
n = @(x,z) (x/ sqrt(f(x) - m(x,z)));
l = 0;
p = 1;
T = @(x)integral(@(z)n(x,z), l , p)
A13 = T(0.95);
A14 = T(0.5);
A15 = T(0.01);

