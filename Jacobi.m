function [x, k, rho] = Jacobi(A, b, x0, tol)
% Computes the solution to Ax = b using Jacobi iteration
% Use tol as the tolerance/accuracy level

D = diag(diag(A));
Dinv = inv(D);

k = 0;
r = b - A*x0;
%re = Dinv*(b - A*x0);
n = length(b);
while norm(r, 1) > tol
%while norm(re, 1) > tol
    x = x0 + Dinv*(b - A*x0);
    k = k + 1;
    x0 = x;
    %re = Dinv*(b - A*x0);
    r = b - A*x;
    if k > 100
        T = eye(n) - Dinv*A;
        rho = max(abs(eig(T)));
        return
    end
end

T = eye(n) - Dinv*A;
rho = max(abs(eig(T)));