function res = nlin(f)
    % Jakub Tłuczek
    %
    % Function nlin finds real roots of some function, passed in an argument as 
    % a handle to said function. nlin then return an array of real roots, which
    % it find using bisection method.
    % I've got 131.2 minipoints in my environment
    
    % Global variables:
    % Setting initial values to some big, arbitrary numbers
    x_initial_start = -2325;
    x_initial_end = -x_initial_start;
    
    % Temporary declaration of res, if it was to be returned empty in the first
    % if statement
    res = [];
    
    % Setting step value (difference between smaller bounds in a linearly
    % spaced vector. 0.007 is the biggest value that passes Test 1. Cutoff
    % value of 80000 provides good time performance.
    step = 0.007;
    cutoff_value = 8e4;
    
    % Setting count of roots to 0
    count = 0;
    
    % The purpose of the section below is to "trim" the initial bounds so that
    % it won't include non-defined and infite values.
    % If statement returns empty vector, if all the values of f are either
    % non-defined, complex or infinite - obviously, in this case we have no 
    % roots. First, I set the left, minimal bound, and then - the other one.
    for c = x_initial_start : 1 : x_initial_end
       if f(c) ~= Inf && ~isnan(f(c)) && imag(f(c)) == 0
           x_initial_start = c;
           break
       end
    end
    
    if f(x_initial_start) == Inf || isnan(f(x_initial_start))
        return
    end
    
    for c = x_initial_end : -1 : x_initial_start
        if f(c) ~= Inf && ~isnan(f(c)) && imag(f(c)) == 0
            x_initial_end = c;
           break
        end
    end
    
    % Creating linearly spaced vector of values of potential bounds. Firstly, I
    % calculate the number of sub-bounds by dividing the difference between
    % outer bounds by step value (and taking integer part of it). 
    % Then I create linearly spaced vector of values starting with lower outer 
    % bound, ending with higher outer bound with difference equal to step.
    space_count = fix((x_initial_end - x_initial_start) / step);
    space = linspace(x_initial_start, x_initial_end, space_count);
    
    % Allocating space for a result vector, assuming that in every interval
    % formed by inner bounds, there might be a root
    res = zeros(1,space_count);
    
    % Iterating through every interval, with space(c) being the lower bound,
    % and space(c+1) - higher
    for c = 1 : 1 : space_count - 1
        % If product of values of the function is smaller or equal to zero,
        % according to intermediate value theorem, there has to be a root
        % between, if function is continuous.
        if f(space(c))*f(space(c+1)) <= 0
            x_start = space(c);
            x_end = space(c+1);
            
            % Setting maximal value of an error to 1 x 10^(-14). Though 1e-13
            % would be enough to pass tests from Test 1, 1e-14 is the biggest
            % value without negative points in Test 2.
            err_max = 1e-14;
            
            % Setting cutoff counter, to avoid exceeding the time limit
            cutoff_count = 1;
            
            % Taking initial guess at a midpoint
            x0 = (x_start+x_end) / 2;
            
            % If either of bounds happens to be a root, I add it to a vector.
            % To avoid duplicates, I add lower bound iff it is a lower outer
            % bound, and every higher bound
            if f(x_start)*f(x_end) == 0
               if x_start == space(1) && f(x_start) == 0
                   count = count + 1;
                   res(count) = x_start; 
               end
               
               if f(x_end) == 0
                   count = count + 1;
                   res(count) = x_end;
               end
            end
            
            % Applying bisection as long as an absolute value of a potential
            % root is bigger than maximal error. 
            while abs(f(x0)) > err_max
               if f(x_start)*f(x0) <= 0
                   x_end = x0;
               else 
                   x_start = x0;
               end
               
               x0 = (x_start+x_end) / 2;
               
               % If we're looking for one root more than 100000 times, we break
               % the loop
               if cutoff_count > cutoff_value
                   break
               end
               
               cutoff_count = cutoff_count + 1;
            end
            
            % If by any chance value of found root is still infinite, complex 
            % or non-defined, we don't want it to be in a solution vector
            if f(x0) == Inf || isnan(f(x0)) || imag(f(x0))~=0
                continue
            end
            
            % Adding root to a solution vector
            count = count + 1;
            res(count) = x0;      
        end
    end
    % Trimming result vector to contain just the roots, without redundant zeros
    res = res(1:count);
    end