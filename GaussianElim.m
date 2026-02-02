function [x] = GaussianElim(A,b)
% Computes the solution to Ax = b using Gaussian Elimination

%Assume A is size nxn
n = length(b);
Ab = [A b];

%Reduce matrix into lower triangular
for j = 1:n
    for i = j+1:n
        M = Ab(i,j)./Ab(j,j);
        Ab(i,:) = Ab(i,:) - M.*Ab(j,:);
    end
end

x = zeros(n,1);
sum = 0;

%Perform backwards substitution for the reduced matrix to get x vector
for i = n:-1:1
    for j = n:-1:i
        if i == n
            x(i) = Ab(i, n+1)/Ab(i, i);
        else
            sum = sum + Ab(i, j).*x(j);
        end
    end
    x(i) = (Ab(i, n+1) - sum)/Ab(i, i);
    sum = 0;
end