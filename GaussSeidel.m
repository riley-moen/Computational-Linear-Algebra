function [x, k, rho] = GaussSeidel(A, b, x0, tol)
% COmputes the solution to Ax = b using Gauss-Seidel numerical method
% Use tol as the tolerance/accuracy level

n = length(b);
E = tril(A);
Einv = inv(E);

k = 0;
r = b - A*x0;
re = Einv*(b - A*x0);

while norm(r, 1) > tol
% while norm(re, 1) > tol
    x = x0 + Einv*(b - A*x0);
    k = k + 1;
    x0 = x;
    re = Einv*(b - A*x0);
    %r = b - A*x;
    if k > 100
        T = eye(n) - Einv*A;
        rho = max(abs(eig(T)));
        return
    end
end

T = eye(n) - Einv*A;
rho = max(abs(eig(T)));