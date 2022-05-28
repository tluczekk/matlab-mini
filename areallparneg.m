function out = areallparneg(A)
% Jakub Tłuczek
%
% Function areallparneg(A) checks, if in matrix A each row and column
% contains at least one negative value. If it is the case, function returns
% logical 1, and 0 if it has at least one row or column without negative
% value, or dimensions of provided matrix are not correct

    % Getting dimensions of matrix A
    [rows, cols] = size(A);
    % Selecting the bigger value for upper bound of for loop
    iters = max(rows, cols);
    % Flags that show, whether non partial-negative row/column was found
    % and if dimensions of A are alright
    non_par_neg_found = 0;
    dims_correct = 0;
    % Condition checking if dimensions of A are correct
    if (rows > 0) && (cols > 0)
        % Setting the flag if condition is cleared
        dims_correct = 1;
        % Main loop, which checks i-th row and i-th column at the same time
        % Conditions are required if dimensions are irregular, and thus
        % prevent program from checking non-existen row or columns
        % If the minimal value of row or column is negative, it means it
        % does have at least one negative element, making it
        % partial-negative
        % In case the minimal value is greater or equal to 0, appropriate
        % flag is set
        for i = 1 : iters
            if i <= rows
                temp_r = A(i, :);
                if min(temp_r >= 0)
                    non_par_neg_found = 1;
                end
            end
            
            if i <= cols
                temp_c = A(:, i);
                if min(temp_c >= 0)
                    non_par_neg_found = 1;
                end
            end
        end
    end
    % Program returns 0 if dimensions are incorrect, or at least one row
    % or column is not partial-negative, and 1 if each of rows and columns
    % fulfills the requirement of having at least one negative element
    if (dims_correct == 0) || (non_par_neg_found == 1)
        out = 0;
    else
        out = 1;
    end
end