function Cons = setCons(eg)
% SETCONS  Problem dimensions for example eg (see getSystem.m).
%   Cons.n   number of states        Cons.T   number of samples
%   Cons.df  degree of f             Cons.dg  degree of g

switch eg
    case 1
        Cons.n 	  = 2;	    % # of state
        Cons.T    = 40;	    % # of sample
        Cons.df   = 1;      % degree of f
        Cons.dg   = 0;      % degree of g
    case 2
        Cons.n 	  = 2;	    % # of state
        Cons.T    = 10;     % # of sample
        Cons.df   = 2;      % degree of f
        Cons.dg   = 0;      % degree of g
    case 4
        Cons.n 	  = 3;	    % # of state
        Cons.T    = 8;	    % # of sample
        Cons.df   = 3;      % degree of f
        Cons.dg   = 0;      % degree of g
end

end
