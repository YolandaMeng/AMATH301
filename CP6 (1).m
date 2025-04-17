
%%% Problem 1
data = readmatrix('lynx.csv');
t = data(1, :);
pop = data(2, :);
%%% Don't delete the lines above when submitting to gradescope

%%% Replace the value of the population for the years given in the assignment file and save it as A1
pop(1956-1946+1)=34;
pop(1974-1946+1)=27;
A1=pop;

%%% Calculate the value of the cubic spline interpolation of the data at t = 24.5 using the interp1 function.  Save this as A2.
interp_value = spline(t, A1, 24.5);
A2 = interp_value;

%%% Use polyfit to calculate the coefficients for A3, A5, and A7
%%% Use norm to calculate the error for A4, A6, and A8

A3 = polyfit(t, A1, 1);
A4 = norm(polyval(A3, t) - A1);
A5 = polyfit(t, A1, 2);
A6 = norm(polyval(A5, t) - A1);
A7 = polyfit(t, A1, 10);
A8 = norm(polyval(A7, t) - A1);



%%% Problem 2
data = readmatrix('CO2_data.csv');
t = data(1, :);
co2 = data(2, :);
%%% Don't delete the lines above when submitting to gradescope

%%% Use polyfit to calculate the coefficients for A9
%%% Use norm to calculate the error for A10
A9 = polyfit(t,co2,1);
A10 = norm(polyval(A9,t) - co2);

%%% Fit the exponential

a = 260;
n = co2 - a;
lnn = log(n);
x = [t', ones(length(t),1)];
b = x \ lnn';
c = exp(b(2));
r = b(1);
n = c* exp(r* t)+a;

A11 = [c,r,a];
A12 = norm(n-co2);


%%% Fit the sinusoidal
%%% There are a few different ways to do this, and we will refrain from giving away the answer to this part.  The class has been doing loops for a while now, so this part should be doable, albeit a little tricky.  We can however check to see if there are any bugs that we can spot.

% Seasonal Adjustments
% Assume the period of seasonal oscillations is one year, so B = 2pi
B = 2 * pi;

% Initialize sum of amplitudes
s = 0;

% Loop through the data by years to calculate the average amplitude
for year = 1:62
    start = (year-1)*12 + 1;
    endin = year*12;
    yearly = co2(start:endin) - (a * exp(r * t(start:endin)) + b);
    maxVal = max(yearly);
    minVal = min(yearly);
    s = s + (maxVal - minVal);
end

% Calculate average amplitude
A = s / (2 * 62); % Dividing by 2 to get average amplitude from peak-to-peak

A13 = [A, B];

% Calculate the fit including the seasonal component
f = a * exp(r * t) + b + A * sin(B * t);
error_seasonal = norm(f - co2, 2); % l2 norm of the error
A14 = error_seasonal;