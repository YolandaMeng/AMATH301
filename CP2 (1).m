%%% Problem 1

%%% Initialize A as a 20 by 21 matrix of zeros (Week 2 Lecture 1)
%%% To fill in the matrix create a nested for loop (Week 2 Lecture 1)
%%% Save the matrix A as the variable A1

for i = 1:20
    for j = 1:21
        A(i,j) = 1/(i*j);
    end
end
A1 = A;




%%% Let B equal A, and set the entire 15th row as zero (Week 2 Lecture 1)
%%% Do the same thing for the entire 16th column
%%% Save the matrix B as A2
B = A;
B(15, 1:21) = zeros;
for j = 16
    for i = 1:20
        B(i,j) = 0;
    end
end
A2 = B;
    



%%% For A3, since we want the last few columns/rows you want to use the end
%%% command (Week 2 Lecture 1)
A3 = B(end-2:end, end-4:end);
A3



%%% Set A4 as the 10th column of B (Week 2 Lecture 1)
A4 = B(1:20, 10);
A4



%%% Problem 2

%%% For A5 and A6 it's exactly like Week 2 Theory lecture.
sum = 0;
for n = 1:20
    sum = 1/n + sum;
end
A5 = sum;

sum = 0;
for n = 1:200
    sum = 1/n +sum;
end
A6 = sum;




%%% For A7 through A10 you're still doing a Sum as you did for A5 and A6
%%% but now you want to break out of the loop when the sum surpasses 10
%%% for A7 and A8, and 20 for A9 and A10
%%% (very similar to Week 2 Lecture 2 Fibonacci)
sum = 0;
n = 0;
while sum < 10
    n = n+1; 
    sum = 1/n +sum;
end
A7 = n;
A8 = sum;

sum = 0;
n = 0;
while sum < 20
    n = n+1; 
    sum = 1/n +sum;
end
A9 = n;
A10 = sum;





%%% Problem 3
%%% First go to the bottom of this document to create a function because in
%%% MATLAB functions go at the end of the m-file, then come back here.

%%% After you have made a function at the end of the m-file you can start
%%% assigning your variables here.

%%% Set N and x0 according to the assignment file





%%% For each of the next three set its respective r value.  Then set the
%%% vector x as the iterates of the logistic map using the output of the
%%% function you created.

%%% Write the code for A11 and A12 here
a = 2.75;
x = zeros(1,100);
x(1) = 0.2;

for n = 1:99
    x(n+1)= a* x(n)* (1-x(n));
end
A11 = x;

behavior = 1;
x_s = std(x);
A12 =[x_s, 1];

%%% Do the same as above except for A13 and A14

for n = 1:99
    x(n+1) = a* x(n)* (1 - x(n));
end
A13 = x;

behavior = 2;
x_s = std(x);
A14 =[x_s, 2];


%%% Do the same as above except for A15 and A16

for n = 1:99
    x(n+1)= a*x(n)*(1 - x(n));
end
A15 = x;

behavior = 3;
x_s = std(x);
A16 =[x_s, 3];



%%% Create a function here (Week 2 Lecture 3).  The function will take r, 
%%% x0, and N as inputs and output a vector x of all N iterates starting at
%%% x0.  Inside the function create a loop that calculates the value of the
%%% logistic map and saves it in its respective entries of x.


