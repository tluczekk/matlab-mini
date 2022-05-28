function[c, s, r] = givens_rotation(a, b)
% Jakub Tłuczek
%
% Source:
% Discontinuous Plane Rotations and the Symmetric Eigenvalue Problem
% E. Anderson, Lockheed Martin Services, 2000
% available at: http://www.netlib.org/lapack/lawnspdf/lawn150.pdf
% Function looks for c and s components of Given's rotation 2x2 matrix,
% which are the solution to equation [c -s; s c] * [a; b] = [r; 0] and r
% is equal to sqrt(a^2 + b^2)
if b == 0
    c = sign(a);
    if (c == 0)
        c = 1.0;
    end
    s = 0;
    r = abs(a);
elseif a == 0
    c = 0;
    s = sign(b);
    r = abs(b);
elseif abs(a) > abs(b)
    t = b / a;
    u = sign(a) * sqrt(1 + t * t);
    c = 1 / u;
    s = c * t;
    r = a * u;
else
    t = a / b;
    u = sign(b) * sqrt(1 + t * t);
    s = 1 / u;
    c = s * t;
    r = b * u;
end
end