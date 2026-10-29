function [f,g] = getSystem(eg,vars)
% GETSYSTEM  True system x+ = f(x) + g*u used to generate data.
%   [f, g] = getSystem(eg, vars) evaluates f at vars (symbolic, sdpvar or
%   numeric).
%   eg = 1: Example 1 (linear), eg = 2: Example 3 (polynomial),
%   eg = 4: Example 2 (3D polynomial).

switch eg
    case 1
        % linear
        A = [0, 1; -2, -5];
        f = A*vars;
        g = [0; 1];
    case 2
        % poly
        x1 = vars(1);
        x2 = vars(2);
        f = [x2; -x1 + 1/3*x1^2 - x2];
        g = [0; 1];
    case 4
        % 3d poly
        x1 = vars(1);
        x2 = vars(2);
        x3 = vars(3);
        f = [x2^2+x3; x1*x2^2 + x3; 0];
        g = [2;-1;1];
end

end
