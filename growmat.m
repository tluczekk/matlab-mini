function mtr = growmat(x)
% Jakub Tłuczek
%
% Function growmat(x) creates matrix with x rows and x columns, in which
% each column is an arithmetic sequence of length x, starting with 1 and
% with common difference of c, where c denotes number of the column
    
    % Creating an empty matrix
    mtr = [];
    % Inserting column one by one
    % First as a row, and then transposing the matrix in the end
    for c = 1 : x
        % Creating column vector starting with 1, with step c, and ending
        % when last element is going to be bigger than x*c - it will be 
        % exactly 1 + x*c
        temp = [1 : c : x*c];
        % Adding new vector transposed to the matrix (temporarily as a row)
        mtr = [mtr; temp];
    end
    % Transposing the matrix, so that sequences are represented by columns
    mtr = transpose(mtr);
end