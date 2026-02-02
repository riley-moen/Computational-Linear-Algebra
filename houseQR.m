function [V, R] = houseQR(A)
% Computes the Householder QR factorization of an (m,n)-matrix A.
% R is an upper triangular matrix
% V is a lower triangular matrix of reflection vectors
% Use the function Q = formQ(V) to compute the orthogonal matrix Q.

% Get m and n from A 
[m, n] = size(A);
% Fill V and R
for k=1:n
    x = A(k:m, k);
    e1 = zeros(length(x), 1); e1(1) = 1;
    V(k:m, k) = x + sign(x(1))*norm(x, 2)*e1;
    V(:, k) = V(:, k) / norm(V(:, k), 2);
    A(k:m, k:n) = A(k:m, k:n) - 2*V(k:m, k)*(V(k:m, k)'*A(k:m, k:n));
end
R = A;