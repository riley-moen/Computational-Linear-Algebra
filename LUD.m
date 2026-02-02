function [L, U, x] = LUD(A,b)
% Computes the LU decomposition of matrix A to solve Ax = b

n = length(b);
L = eye(n);

%Make L a lower triangular matrix while also making M upper triangular
for j = 1:n
    for i = j+1:n
        M = A(i,j)./A(j,j);
        A(i,:) = A(i,:) - M.*A(j,:);
        L(i, j) = M;
    end
end

%Set U = A since A is upper triangular
U = A;
sum = 0;

%Solve for y by doing forwards substitution
y = zeros(n, 1);
for i = 1:n
    for j = 1:i
        if i == 1
            sum = sum + 0;
        else
            sum = sum + L(i, j).*y(j);
        end
    end
    y(i) = (b(i) - sum);
    sum = 0;
end

x = zeros(n,1);
sum = 0;
%Find x by doing backwards substitution
for i = n:-1:1
    for j = n:-1:i
        if i == n
            x(i) = U(i, n)/y(i);
        else
            sum = sum + U(i, j).*x(j);
        end
    end
    x(i) = (y(i) - sum)./U(i, i);
    sum = 0;
end
