%%% Problem 1
%%% Implement the Bisection method as we did in the Week 5 Coding Lecture
%%
t = linspace(0, 30, 100); % Adjust the time range accordingly
c = 1.3 * (exp(-t/11) - exp(-4*t/3));
plot(t, c);
xlabel('Time (hours)');
ylabel('Concentration');
title('Concentration vs. Time');

a =  1;
b = 3;
tol = 1e-8;

for k = 1:1000000
    if t1 ==0
        break
    elseif sign(t1) == sign(da)
        a = t;
    else
        b = t;
    end
    t = (a+b) / 2;
    t1 = 1.3* (-1/11* exp(-t/11)+ 4/3 *exp(-4 * t /3));
    da = 1.3* (-1/11* exp(-a/11)+ 4/3 *exp(-4 * a /3));
    db = 1.3* (-1/11* exp(-b/11)+ 4/3 *exp(-4 * b /3));
    A1 = t;
    A2 =  1.3* (exp(-t/11)-exp(-4 * t /3));
    A3 = abs( 1.3* (-1/11* exp(-t/11)+ 4/3 *exp(-4 * t /3)));
end





%%% Problem 2
%%% Implement Newton's method as we did in the Week 5 Coding Lecture
%%
x0 = 2;
iterations = 1;
x = x0;
tolerance = 1e-8;

while abs(2*x) > tolerance
        x = x - 2*x / 2;
        iterations = iterations + 1;
end

 A4 = iterations;
 A5 = x;

x0 = 2;
iterations = 1;
x = x0;
tolerance = 1e-8;

while abs(500*x^499) > tolerance
        x = x - 500*x^499 / (249500*x^498);
        iterations = iterations + 1;
end

A6 = iterations;
A7 = x;

 x0 = 2;
iterations = 1;
x = x0;
tolerance = 1e-8;

while abs(1000*x^999) > tolerance
        x = x - (1000*x^999)/ (999000*x^998);
        iterations = iterations + 1;
end

 A8 = iterations;
 A9 = x;
 


