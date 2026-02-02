function [l, v, k] = powMethod(A, x0, tol)
% Computes the dominant eigenvalue of A using the power method

x0 = x0/norm(x0);
l0 = x0'*A*x0;
v = A*x0/norm(x0);
l = v'*A*v;
k = 1;
fprintf('l = %4.4f, %4.0f iterations\n', [l, k])

while abs(l0 - l) > tol
    x0 = v;
    l0 = l;
    v = (A*x0)/norm(A*x0);
    l = v'*A*v;
    k = k + 1;
    fprintf('l = %4.4f, %4.0f iterations\n', [l, k])
end