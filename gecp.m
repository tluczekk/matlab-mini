function [xfin, U, L, P, Q] = gecp(A, B)
% Usage of the gecp function to calculate X from the equation AX = B:
% Function takes two arguments - matrices A and B from the equation AX = B 
% The output may vary, depending which results we want to extract from the
% function. We can call it for example:
% X = gecp(A, B)
% which will return just the calculated X matrix, or we can call it even:
% [X, U, L, P, Q] = gecp(A, B)
% which returns X, as well as all matrices from PAQ = LU equation
% error is printed on the output as well
[n, n2] = size(A) % extracts A dimensions
[n1, m] = size(B) % extracts B dimensions
if n ~= n2 || n1 ~= n || m > n % checks whether dimensions are valid
    disp('Dimensions not correct'); % if not returns the function
    return;
end
p = 1:n; % for P and Q evaluation
q = 1:n;
A_aug = [A, B]; % creating augmented matrix
col_per = zeros(n, 1) % vector for keeping track of column swaps
for i = 1 : 1 : n
    col_per(i, 1) = i;
end
for k = 1 : n-1
    [maxc, rowindices] = max(abs(A_aug(k:n, k:n))); % looking for pivot
    [maxm, colindex] = max(maxc);
    row = rowindices(colindex) + k - 1; col = colindex + k - 1;
    A_aug([k, row], :) = A_aug([row, k], :); % row swap
    A_aug(:, [k, col]) = A_aug(:, [col, k]); % column swap
    temp_col = col_per(k, 1); % column swap tracking
    col_per(k, 1) = col_per(col, 1);
    col_per(col, 1) = temp_col;
    p([k, row]) = p([row, k]); % P and Q updates
    q([k, col]) = q([col, k]);
    if A_aug(k, k) == 0 % if pivot is 0 there is no sense to add rows
        break
    end
    A_aug(k+1:n, k) = A_aug(k+1:n, k)/A_aug(k, k); % adding rows
    i = k+1:n;
    j = k+1:n+m; % we don't zero values below pivot to save them for L
    A_aug(i, j) = A_aug(i, j) - A_aug(i, k) * A_aug(k, j);
end

A = A_aug(1:n, 1:n); % extracting left side (A) of augmented matrix
B = A_aug(1:n, n+1:n+m); % extracting right side (B) of an augmented matrix
L = tril(A, -1) + eye(n); % extracting L from changed A
U = triu(A); % extracting U from changed A
P = eye(n);
P = P(p, :); % create P matrix based on 'tracking' vector
Q = eye(n);
Q = (:, q); % the same as above

% calculating x with entries not on the right places:
x = zeros(n, m);
for c = 1 : 1 : m
    res = B(:, c);
    for j = n : -1 : 1
        if(U(j, j) == 0)
            error('singular matrix')
        end
        x(j, c) = res(j) / U(j, j);
        res(1:j-1) = res(1:j-1) - U(1:j-1, j) * x(j, c);
    end
end

xfin = zeros(n, m); % correcting it thanks to col_per vector
for i = 1 : 1 : n
    for j = 1 : 1 : m 
        xfin(col_per(i), j) = x(i, j);
    end
end