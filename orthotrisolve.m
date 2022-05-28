function x = orthotrisolve(a, b, c, r)
% Jakub Tłuczek
% Function orthotrisolve takes four vectors as an input - b and r are
% n-element vectors, while a and c have n-1 elements each. Vector x is an
% n-element vector which is a solution to the equation A*x=r, where A is a
% tridiagonal matrix with a, b and c being its diagonals, going upwards.

% Getting length of the longest vector
n = length(b);

% preparing memory for d-vector which represent another diagonal vector
% created by Given's rotation, and x result vector
d = zeros(n-2, 1);
x = zeros(n, 1);
% Loop of given's rotations, applied to two rows at a time
for i = 1 : n-1
    % Getting Given's matrix which when applied to first nonzero column
    % of two rows in question, will zero the element under the diagonal
    % (i.e. from "a" vector), and set b(i) element to
    % sqrt(a(i)^2 + b(i)^2)
    [giv_c, giv_s, res] = givens_rotation(a(i, 1), b(i, 1));
    givens_matrix = [giv_c -giv_s; giv_s giv_c];
    b(i, 1) = res;
    
    % Applying Given's matrix to second nonzero column, i.e.
    % [c(i);b(i+1)]
    temp_bc = [c(i, 1); b(i+1, 1)];
    res_bc = givens_matrix * temp_bc;
    c(i, 1) = res_bc(1, 1);
    b(i+1, 1) = res_bc(2, 1);
    
    % Applying Given's matrix to the third non-zero column, if there's
    % any, i.e. if we're not considering last and penultimate rows.
    % We're considering then [d(i); c(i+1)] with d(i) being initialy
    % equal to 0
    if i < n-1
        temp_cd = [0; c(i+1, 1)];
        res_cd = givens_matrix * temp_cd;
        c(i+1, 1) = res_cd(2, 1);
        d(i, 1) = res_cd(1, 1);
    end
    
    % Applying Given's matrix to the r vector, i.e. to [r(i); r(i+1)]
    temp_rr = [r(i, 1); r(i+1, 1)];
    res_rr = givens_matrix * temp_rr;
    r(i, 1) = res_rr(1, 1);
    r(i+1, 1) = res_rr(2, 1);
end

% Solving for x using backwards substitution, iterating from n to 1
for i = n : -1 : 1
    if i == n
        x(i, 1) = r(i, 1) / b(i, 1);
    elseif i == n-1
        x(i, 1) = (r(i, 1) - c(i, 1) * x(i+1, 1)) / b(i, 1);
    else
        x(i, 1) = (r(i, 1) - c(i, 1) * x(i+1, 1) - d(i, 1) * x(i+2, 1)) / b(i, 1);
    end
end
end