function [f,g] = getSystem(eg,vars)

x1 = vars(1);
x2 = vars(2);

switch eg
    case 1
        % linear
        A = [0, 1; -2, -5]; 
        f = A*vars; % random linear system
    case 2
        % ranzter's "On Analysis and Synthesis of Safe Control Laws"
        f = [x2; -x1 + 1/3*x1^3 - x2]; 
end

g = [0; 1];
 
end