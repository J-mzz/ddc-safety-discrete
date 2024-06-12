function [f,g] = getSystem(eg,vars)

switch eg
    case 1
        % linear
        A = [0, 1; -2, -5]; 
        f = A*vars; % random linear system
        g = [0; 1];
    case 2
        % poly
        x1 = vars(1);
        x2 = vars(2);
        f = [x2; -x1 + 1/3*x1^2 - x2]; 
%         f = [x2; -x1 + 1/3*x1^3 - x2]; 
        g = [0; 1];
    case 3
        % linear scalar
        f = -2*vars;
        g = 1;
    case 4
        % 3d poly
        x1 = vars(1);
        x2 = vars(2);
        x3 = vars(3);
        f = [x2^2+x3; x1*x2^2 + x3; 0];
        g = [2;-1;1];
end

 
end