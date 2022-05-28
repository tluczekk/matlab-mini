function [p, val] = Steffensen(n, a, p0, tol)
% Function takes following arguments:
% n - from the formula of the function
% a - series from the formula
% p0 - initial guess
% tol - tolerance
% It produces the output consisting of found root, which we can later 
% substitute as an argument of our function, and array of approximations
format compact
format long
val = [p0];
for i=1:1000    %1000 is an arbitrary max, in which we expect function to
                % converge
    f = nasza(n, a, p0);
    g = (nasza(n, a, p0+f) - f)/f;
    if g == 0      %display result if g hits 0
        p
        break
    end
    p = p0 - f/g;   % calculating the next term
    val = [val, p];
    if abs(nasza(n,a,p)) < tol % display value if within tolerance
        p
        break
    end
    p0 = p;         % assign new value and proceed to next iteration
end
if abs(nasza(n,a,p)) > tol  % If failed to converge, display message
    'Failed to converge in 1000 iterations';
end