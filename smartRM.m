function [Q, err, K] = smartRM(f, a, b, tol, m0, Kmin, Kmax)
% Jakub Tłuczek
%
% Function smartRM calculates an integral of the function f on the interval
% [a, b] using the Romberg method of numerical integration. What makes it
% different however, is that it can stop when certain condition is smaller
% than a tolerance tol, specified in arguments. Function smartRM will do at
% least Kmin interpolation steps, and Kmax at maximum. Output variables are
% the result of the integration Q, err which is equal either to the value
% of aforementioned condition or to tol value, and K which indicates number
% of extrapolation steps performed.

% Predefining Romberg table's memory
RM = zeros(Kmax + 1);

% Iterating through each row of Romberg table
for i = 1 : Kmax + 1
    % Calculating number of trapezoids required to calculate an integral,
    % and height of each one respectively. T_0,0 is calculated using m0
    % trapezoids, with each trapezoid being divided in 2 in the next step.
    num_steps = m0 * 2^i;
    step = (b - a)/(num_steps);
    % Iterating through each column of Romberg's table
    for j = 1 : i
        % If in first column, integral must be calculated by summing areas
        % of the trapezoids by hand
        if j == 1
            % it represents first point for value of function to be
            % evaluated, the other being it + step. it_num is the counter
            % of steps performed. Areas of trapezoids are calculated until
            % it reaches num_steps
            it = a;
            it_num = 1;
            zero_sum = 0;
            while it_num <= num_steps
                % Calculating area of each trapezoid and adding it to the
                % sum
                step_area = step*(f(it)+f(it+step))/2;
                % Skipping trapezoid if area is infinite or undefined.
                % Every trapezoid is divided in two in the next step, so
                % we're going to get a better approximation then in case of
                % non-continuous functions
                if step_area ~= Inf && ~isnan(step_area)
                    zero_sum = zero_sum + step_area;
                end
                it = it + step;
                it_num = it_num + 1;
            end
            % Setting T_i-1,j-1 value to the sum of the areas of the
            % trapezoids
            RM(i, j) = zero_sum;
        else
            % If not in the first column, T_i-1,j-1 is calculated using the
            % formula below
            RM(i, j) = (4^(j-1) * RM(i, j-1) - RM(i-1, j-1))/(4^(j-1) - 1);
        end
    end
    if i>=Kmin+1
        % When number of extrapolation steps reaches Kmin (in other words,
        % i reaches Kmin + 1, since we haven't performed extrapolation for
        % i=1), we might start to check if the condition specified in the
        % task description is fulfilled, i.e. if
        % abs(RM(i,i) - RM(i-1, i-1))/max([1 abs(RM(i,i))]) < tol
        % for i in {K-1, K}, with K in this case being i-1.
        % When accessing Romberg table, 1 has to be added to index because
        % of the MATLAB's convention of enumeration starting with 1, not 0.
        % As stated in  the task description, it is assumed that:
        % Kmin >= 2 && Kmax >= Kmin + 2
        condition = abs(RM(i,i)-RM(i-1,i-1))/max([1 abs(RM(i,i))]);
        condition2 = abs(RM(i-1,i-1)-RM(i-2,i-2))/max([1 abs(RM(i-1,i-1))]);
        if condition < tol && condition2 < tol
            % If conditions are fulfilled, program returns
            err = tol;
            K = i - 1;
            Q = RM(i,i);
            return;
        end
    end
end
% If Kmax steps of interpolation are reached, program returns
err = abs(RM(Kmax+1,Kmax+1)-RM(Kmax, Kmax))/max(1, abs(RM(Kmax+1,Kmax+1)));
K = Kmax;
Q = RM(Kmax+1, Kmax+1);
end