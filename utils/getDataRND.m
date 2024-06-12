function x = getDataRND(eg,eps,N)
% Randomly sample points on the plane, generate xdot with input and noise. 
% System dynamics: x_dot = f(x) + g(x)u(x) + noise
% u \in [-1,1]

box = 5;
switch eg
    case {1,2}
        syms x1 x2
        vars = [x1;x2];
        [f,g] = getSystem(eg,vars);
        
        x = [];
        for i = 1:N
            
            u = 2*rand(1,1)-1;
        
	        noise = eps*(2*rand(2,1)-1);
            
            x0 = (2*rand(2,1)-1) * box;
            
            x1 = x0(1);
            x2 = x0(2);
            
            xdot = subs(f+g*u+noise);
            
            curr = [u,i,x1,x2,xdot(1),xdot(2)];
            
            x = [x;curr];
        end

    case 3
        syms x1
        vars = x1;
        [f,g] = getSystem(eg,vars);
        
        x = [];
        for i = 1:N
            
            u = 2*rand(1,1)-1;
        
	        noise = eps*(2*rand(1,1)-1);
            
            x0 = (2*rand(1,1)-1) * box;
            
            x1 = x0;
            
            xdot = subs(f+g*u+noise);
            
            curr = [u,i,x1,x1,xdot,xdot];
            
            x = [x;curr];
        end

   case 4
        syms x1 x2 x3
        vars = [x1;x2;x3];
        [f,g] = getSystem(eg,vars);
        
        x = [];
        for i = 1:N
            
            u = 2*rand(1,1)-1;
        
	        noise = eps*(2*rand(3,1)-1);
            
            x0 = (2*rand(3,1)-1) * box;
            
            x1 = x0(1);
            x2 = x0(2);
            x3 = x0(3);

            xdot = subs(f+g*u+noise);
            
            curr = [u,i,x1,x2,x3,xdot(1),xdot(2),xdot(3)];
            
            x = [x;curr];
        end
end

x = double(x);

end