function Cons = setCons(eg)

switch eg
    case 1
        Cons.n 	  = 2;	    % # of state
        Cons.T    = 40;	    % # of sample
        Cons.drho = 4;	    % degree of rho
        Cons.dpsi = 4;      % degree of psi
        Cons.df   = 1;      % degree of f
        Cons.dg   = 0;      % degree of g
    case 2
        Cons.n 	  = 2;	    % # of state
        Cons.T    = 200;	% # of sample
        Cons.drho = 4;	    % degree of rho
        Cons.dpsi = 4;      % degree of psi
        Cons.df   = 3;      % degree of f
        Cons.dg   = 0;      % degree of g
    case 3
        Cons.n 	  = 1;	    % # of state
        Cons.T    = 40;	    % # of sample
        Cons.drho = 4;	    % degree of rho
        Cons.dpsi = 4;      % degree of psi
        Cons.df   = 1;      % degree of f
        Cons.dg   = 0;      % degree of g
end

end