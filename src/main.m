% Testing script for the CGNR + ILQ implementation
% The data folder contains the sparse representations of the coefficient matrix A
% and vector b for a Poisson problem on WinkelStructured geometry

disp("\nReading the data")

A = dlmread("data/linsys_a.dat");
b = dlmread("data/linsys_b.dat");

n_elems = size(b, 1);

A(:, [1, 2]) = A(:, [2, 1]);
A = spconvert(A);

b = [b(:, 1), ones(n_elems, 1), b(:, 2)];
b = spconvert(b);


% Solve the problem with cgnr_ilq
param = 1;
disp(["\nSolving the problem with CGNR + ILQ {", num2str(param), "} ... \n"])

x0 = sparse(zeros(size(b)));

tic
x = cgnr_ilq(A, x0, b, 1000, 0.000001, param, false);
toc

disp(["Residual norm: ", num2str(norm(b - A * x))])


% Solve the problem with regular cg

% disp("\nSolving the problem with CG ... \n")

% x0 = zeros(size(b));

% tic
% x = cg(A, x0, b, 10000, 0.000001);
% toc

% disp(["Residual norm: ", num2str(norm(b - A * x))])


% Solve the problem with cgnr

disp("\nSolving the problem with CGNR ... \n")

x0 = sparse(zeros(size(b)));

tic
x = cgnr(A, x0, b, 10000, 0.000001);
toc

disp(["Residual norm: ", num2str(norm(b - A * x))])



