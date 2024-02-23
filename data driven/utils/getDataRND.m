function x = getDataRND(eg,eps,N)
% Randomly sample points on the plane, generate xdot with input and noise. 
% System dynamics: x_dot = f(x) + g(x)u(x) + noise
% u \in [-1,1]

box = 5;

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

x = double(x);
end