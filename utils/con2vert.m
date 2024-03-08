function [V,nr] = con2vert(A,b,F,d, DEBUG)
% CON2VERT - convert a convex set of constraint inequalities into the set
%            of vertices at the intersections of those inequalities;i.e.,
%            solve the "vertex enumeration" problem. Additionally,
%            identify redundant entries in the list of inequalities.
% 
% V = con2vert(A,b)
% [V,nr] = con2vert(A,b)
% 
% Converts the polytope (convex polygon, polyhedron, etc.) defined by the
% system of inequalities A*x <= b into a list of vertices V. Each ROW
% of V is a vertex. For n variables:
% A = m x n matrix, where m >= n (m constraints, n variables)
% b = m x 1 vector (m constraints)
% V = p x n matrix (p vertices, n variables)
% nr = list of the rows in A which are NOT redundant constraints
%
% ... = con2vert(A,b,F,d)
%   also adds the equality constraints F*x = d in addition to A*x <= b
%   New in ver 1.2
%
% ... = con2vert(A,b,F,d,DEBUG)
%   will give more verbose output.  Set F=[] and d=[] if you do not
%   want to include the equality constraints.
%   New in ver 1.2
% 
% NOTES: (1) This program employs a primal-dual polytope method.
%        (2) In dimensions higher than 2, redundant vertices can
%            appear using this method. This program detects redundancies
%            at up to 6 digits of precision, then returns the
%            unique vertices.
%        (3) Non-bounding constraints give erroneous results; therefore,
%            the program detects non-bounding constraints and returns
%            an error. You may wish to implement large "box" constraints
%            on your variables if you need to induce bounding. For example,
%            if x is a person's height in feet, the box constraint
%            -1 <= x <= 1000 would be a reasonable choice to induce
%            boundedness, since no possible solution for x would be
%            prohibited by the bounding box.
%        (4) This program requires that the feasible region have some
%            finite extent in all dimensions. For example, the feasible
%            region cannot be a line segment in 2-D space, or a plane
%            in 3-D space.
%                NOTE: the updated version can handle affine constraints
%                   as they are explicit in the Fx == d constraints 
%        (5) At least two dimensions are required.
%        (6) See companion function VERT2CON.
%        (7) ver 1.0: initial version, June 2005
%        (8) ver 1.1: enhanced redundancy checks, July 2005
%        (9) Written by Michael Kleder
%        (10) ver 1.2 handles Fx == d constraints and uses linear
%         programming rather than fminsearch to find strictly feasible
%         point. 
%        Explicitly handling equality constraints is necessary since if you
%         tried to implicitly encode them using Fx <= d and -Fx <= -d
%         then this code would fail since the code requires finding a *strictly*
%         feasible point (e.g., an interior point). By explicitly handling
%         the equality constraints, we only need to find a point in the
%         relative interior.
%        Almost all the significant computing time is due to convhulln,
%         which is Matlab's wrapper to qhull (see www.qhull.org).
%         As dimensions increase above 10, this rapidly gets slow.
%        Ver 1.2 programmed by Stephen Becker, May 1 2021,
%         Stephen.Becker@Colorado.edu
%
% LICENSE ISSUE:
%   For ver 1.2, which is the 2021 update by a different author,
%    we are assuming the original version (which doesn't have an explicit
%    license) is now covered by the default BSD license for File Exchange
%     https://www.mathworks.com/matlabcentral/FX_transition_faq.html
%   This is an old code and the author hasn't responded to any
%    of the comments on Mathworks File Exchange for over 15 years,
%    so we think this is OK, but will of course revise if we hear
%    otherwise from the original author Michael Kleder
%   This new version is released under the BSD license as well.
%
% EXAMPLES:
%
% % FIXED CONSTRAINTS:
% A=[ 0 2; 2 0; 0.5 -0.5; -0.5 -0.5; -1 0];
% b=[4 4 0.5 -0.5 0]';
%
% figure(1); clf;   hold on
% grd = linspace(-3,5,500);  [x,y] = meshgrid( grd );
% p=[x(:) y(:)]';   p=double(all((A*p <= b)));
% p =reshape(p,size(x));
% imagesc(grd,grd,p); axis image; colormap summer; xlabel('x'); ylabel('y');
% V=con2vert(A,b);
% plot(V(:,1),V(:,2),'rs','MarkerSize',12,'MarkerFaceColor','r')
% 
% % RANDOM CONSTRAINTS:
% A=rand(30,2)*2-1;
% b=ones(30,1);
%
% figure(2); clf;   hold on
% grd = linspace(-3,5,500);  [x,y] = meshgrid( grd );
% p=[x(:) y(:)]';   p=double(all((A*p <= b)));
% p =reshape(p,size(x));
% imagesc(grd,grd,p); axis image; colormap summer; xlabel('x'); ylabel('y');
% V=con2vert(A,b);
% plot(V(:,1),V(:,2),'rs','MarkerSize',12,'MarkerFaceColor','r')
% 
% % HIGHER DIMENSIONS:
% A=rand(15,5)*1000-500;
% b=rand(15,1)*100;
% V=con2vert(A,b)
% 
% % NON-BOUNDING CONSTRAINTS (ERROR):
% A=[0 1;1 0;1 1];
% b=[1 1 1]';
% figure('renderer','zbuffer')
% hold on
% [x,y]=ndgrid(-3:.01:3);
% p=[x(:) y(:)]';
% p=(A*p <= repmat(b,[1 length(p)]));
% p = double(all(p));
% p=reshape(p,size(x));
% h=pcolor(x,y,p);
% set(h,'edgecolor','none')
% set(h,'zdata',get(h,'zdata')-1) % keep in back
% axis equal
% set(gca,'color','none')
% V=con2vert(A,b); % should return error
% 
% % NON-BOUNDING CONSTRAINTS WITH BOUNDING BOX ADDED:
% A=[0 1;1 0;1 1];
% b=[1 1 1]';
% A=[A;0 -1;0 1;-1 0;1 0];
% b=[b;4;1000;4;1000]; % bound variables within (-1,1000)
% figure('renderer','zbuffer')
% hold on
% [x,y]=ndgrid(-3:.01:3);
% p=[x(:) y(:)]';
% p=(A*p <= repmat(b,[1 length(p)]));
% p = double(all(p));
% p=reshape(p,size(x));
% h=pcolor(x,y,p);
% set(h,'edgecolor','none')
% set(h,'zdata',get(h,'zdata')-1) % keep in back
% axis equal
% set(gca,'color','none')
% V=con2vert(A,b);
% plot(V(:,1),V(:,2),'y.','markersize',20)
%
% % JUST FOR FUN (3D visualization):
% A=rand(80,3)*2-1;
% n=sqrt(sum(A.^2,2));
% A=A./repmat(n,[1 size(A,2)]);
% b=ones(80,1);
% V=con2vert(A,b);
% k=convhulln(V);
% figure
% hold on
% for i=1:length(k)
%     patch(V(k(i,:),1),V(k(i,:),2),V(k(i,:),3),'w','edgecolor','none')
% end
% axis equal
% axis vis3d
% axis off
% h=camlight(0,90);
% h(2)=camlight(0,-17);
% h(3)=camlight(107,-17);
% h(4)=camlight(214,-17);
% set(h(1),'color',[1 0 0]);
% set(h(2),'color',[0 1 0]);
% set(h(3),'color',[0 0 1]);
% set(h(4),'color',[1 1 0]);
% material metal
% for x=0:5:720
%     view(x,0)
%     drawnow
% end

numDigits           = 10;
constraintTolerance = 1e-8;


% New debugging options, ver 1.2
if nargin < 5 || isempty(DEBUG)
    DEBUG = false;
end

if DEBUG
    printf = @fprintf;
else
    printf = @(varargin) 0; % no output
end


% ==== Handle equality constraints via change-of-variables =====
% ( new in Ver 1.2 )
if nargin >= 4 && ~isempty(F) && ~isempty(d)
    % Fx == d constraints have been specified
    % { x : Fx=d } can be written as x = xp + H*y where H=null(F)
    % so then Ax <= d becomes A*xp + A*H*y <= d
    
    % First, find a "particular solution" to shift by
    xp = F\d;
    % Now find a basis for the null space of F
    H  = null(F);
    
    % call this function again, but without equality constraints
    [VV,nr]=con2vert(A*H,b - A*xp, [], [], DEBUG);
    % Finally, postprocess output to undo the change-of-variables
    V = VV*H' + repmat(xp',size(VV,1),1);
    
    return
end



% ==== From here on, there are no explicit equality constraints =====

% ver 1.2: add check for dim >= 2, as require by convhull
if size(A,2) < 2
    error('con2vert:dim','con2vert requires at least 2 dimensions, since it depends on convhulln');
end

c = A\b; % try a easy-to-compute guess and hope it is an interior point
if ~all(A*c < b)
    printf('Finding strictly feasible point\n');
    
    NEW_SOLVER = true;
    if NEW_SOLVER
        % New behavior, ver 1.2
        % Solve the linear program:  min_{x,t} 0*x+1*t s.t. A*x - b <= t*ones
        %   i.e., Find x such that A*x < b as much as possible (with
        %     respect to the worst-case entry).
        % This is 1000% more robust than using fminsearch, so we recommend
        % you use this option.
        dim = size(A,2);
        m   = size(A,1);
        if DEBUG, display = 'final'; else, display = 'off'; end
        opts = optimoptions('linprog','Display',display,'ConstraintTolerance',constraintTolerance);
        [x_and_t,fc,ef] = linprog( [zeros(dim,1);1], [A,-ones(m,1)], b, [],[],[],[],opts );
        c = x_and_t(1:dim);
        t = x_and_t(end); % ought to be negative
        if t > 0
            error('con2vert:NotStrictlyFeasible','Unable to locate a point within the interior of a feasible region.')
        end
    else
        % This is the old behavior
        fcn = @(c) obj(c,{A,b});
        [c,f,ef] = fminsearch(fcn,c); % 2021 modification of original code, same slow derivative-free algorithm
        %[c,f,ef] = fminsearch(@obj,c,'params',{A,b}); % 2005 original code, incompatiable with recent Matlab
        t = -f;
    end
    if ef ~= 1
        error('con2vert:NotStrictlyFeasible','Unable to locate a point within the interior of a feasible region.')
    else
        printf('... found strictly feasible point (feasible by amount %g)\n', -t);
    end
end
b = b - A*c;        % change-of-variables: shift so that the point "c" is now at the origin

% In the following step, because 0 is strictly feasible (A*0 < b ), this means the new b
%   is all positive, so we can safely dvide by b without changing
%   the direction of the inequality:
%D = A ./ repmat(b,[1 size(A,2)]);   % implicitly, b is now all ones
D = bsxfun( @times, A, 1./b );      % ver 1.2, more memory efficient

% -- This is the expensive step --
printf('First call to convhulln\n');
[~,volume2] = convhulln([D;zeros(1,size(D,2))]);
printf('Second call to convhulln\n');
[k,volume1] = convhulln(D);
printf('... done with all convhulln\n');
if volume2 > volume1
    error('Non-bounding constraints detected. (Consider box constraints on variables.)')
end

nr = unique(k(:));
G  = zeros(size(k,1),size(D,2));
o  = ones( size(k,2), 1 );  % this is the RHS, which was scaled to be all ones
if ~DEBUG
    warnState = warning('off','MATLAB:singularMatrix'); % don't freak out user
end
badIndices = [];
for ix = 1:size(k,1)
    F = D(k(ix,:),:);
    G(ix,:)=F\o;
    [~, LASTID] = lastwarn;
    if any( isnan(G(ix,:)) ) || strcmpi(LASTID,'MATLAB:singularMatrix')
        % This comes up via the example from Matt J in 18 Feb 2012
        % https://www.mathworks.com/matlabcentral/fileexchange/7894-con2vert-constraints-to-vertices
        % Ver 1.2 tries to handle this
        %       We could solve (e.g., pinv(F)*o) if o is in range of F,
        %       but there are an infinite number of solutions. Rather,
        %       this condition implies that these are actually NOT
        %       defining a vertex, so we shoud skip this
        badIndices(end+1) = ix;
        printf(' Vertex %d: system not invertible, so removing vertex\n', ix);
    end
    lastwarn(''); % reset
end
G(badIndices,:) = []; % Ver 1.2: remove the fake verticies
if ~DEBUG
    warning( warnState ); % reset back to how it was
end

% V = G + repmat(c',[size(G,1),1]);   % undo the change-of-variables / shift

V1 = round( G*randn(size(G,2),1),  numDigits, 'significant' );
V2 = round( G*randn(size(G,2),1),  numDigits, 'significant' );
[~,I1]=unique(V1,'stable'); % ver 1.2: cheaper way to do it. Also, don't sort
[~,I2]=unique(V2,'stable'); % do it twice to be extra careful with numerics
I   = union( I1, I2 );
if length(I) > max( [length(I1), length(I2)] )
    printf('Scheme for finding unique rows may not have worked\n');
end
% [~,I]=unique(num2str(V,6),'rows'); % Old behavior, slow
G=G(I,:);

V = bsxfun( @plus, G, c' );         % ver 1.2: more memory efficient
%V = G + repmat(c',[size(G,1),1]);  % undo the change-of-variables / shift 
                                    % (*after* finding unique rows, probably 
                                    %  slightly more stable this way )
end

% Only needed if NEW_SOLVER = false (which isn't recommended)
function d = obj(c,params)
    A=params{1};
    b=params{2};
    d = A*c-b;
    k=(d>=-1e-15);
    d(k)=d(k)+1;
    d = max([0;d]);
end


