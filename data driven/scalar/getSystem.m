function [f,g] = getSystem(eg,vars)

switch eg
    case 1
        % linear
        A = [0, 1; -2, -5]; 
        f = A*vars; % random linear system
        g = [0; 1];
    case 2
        % ranzter's "On Analysis and Synthesis of Safe Control Laws"
        x1 = vars(1);
        x2 = vars(2);
        f = [x2; -x1 + 1/3*x1^3 - x2]; 
        g = [0; 1];
    case 3
        % linear scalar
        f = -2*vars;
        g = 1;
end

 
end