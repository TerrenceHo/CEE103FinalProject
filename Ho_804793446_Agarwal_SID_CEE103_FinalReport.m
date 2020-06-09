%% Solving for First Concentration
D = 5000;
U = 100;
k = 2;
L = 100;
cin = 100;
dx = 5;

D2 = 2500;
U2 = 50;

% default parameters
[displacement1, concentration1] = reactor_solver(D, U, k, L, cin, dx);

% varying parameters
[displacement2, concentration2] = reactor_solver(D2, U, k, L, cin, dx);
[displacement3, concentration3] = reactor_solver(D, U2, k, L, cin, dx);

%% Plotting Concentration at each position
% set this to false if you do not want to generate images
save_figures = false;

fig1 = figure(1);
plot(displacement1, concentration1,'-o','LineWidth',2); hold on;
xlabel("Position along reactor (m)");
ylabel("Concentration (mol/L");
title("Concentration at each position in the tank (Baseline)");
set(gca,'FontSize',17);
hold off;
if save_figures == true
    saveas(fig1,'images/fig1.png')
end

fig2 = figure(2);
plot(displacement2, concentration2,'-o','LineWidth',2); hold on;
xlabel("Position along reactor (m)");
ylabel("Concentration (mol/L");
title("Concentration at each position in the tank (D=2500 m^s/hr)");
set(gca,'FontSize',17);
hold off;
if save_figures == true
    saveas(fig2,'images/fig2.png')
end

fig3 = figure(3);
plot(displacement3, concentration3,'-o','LineWidth',2); hold on;
xlabel("Position along reactor (m)");
ylabel("Concentration (mol/L");
title("Concentration at each position in the tank (U=50 m/hr)");
set(gca,'FontSize',17);
hold off;
if save_figures == true
    saveas(fig3,'images/fig3.png')
end

%% Reactor Solver Function
function [xx, yy] = reactor_solver(D, U, k, L, cin, dx)
% A function to solve the steady state reactor problem
% Inputs:
%   D: dispersion coefficient 
%   U: mean velocity through tank
%   k: reaction rate
%   L: length of the tank
%   cin: concentration in the inflow
%   dx: discretization of the length of the reactor
% Outputs:
%   xx: xx is an (N+1) column vector of the node points
%   yy: yy is an (N+1) column vector of the solution values

    a = 0;
    b = L;
    N = L/dx;
    
    function p = reactor_p_func(x)
    p = U./D.*(ones(size(x)));
    end
    function q = reactor_q_func(x)
        q = k./D.*(ones(size(x)));
    end
    function r = reactor_r_func(x)
        r = 0.*x;
    end
    function v = ghost_left_func1(x, h)
        v = (-U/D)*2*h.*(ones(size(x)));
    end
    function v = ghost_left_func2(x, h)
        v = 1.*(ones(size(x)));
    end
    function v = ghost_left_func_rhs(x, h)
        v = (U/D)*2*h*cin.*(ones(size(x)));
    end
    function v = ghost_right_func1(x, h)
        v = 1.*(ones(size(x)));
    end
    function v = ghost_right_func2(x, h)
        v = 0.*(ones(size(x)));
    end
    function v = ghost_right_func_rhs(x, h)
        v = 0.*(ones(size(x)));
    end
    z = ODEBVP_boundary(@reactor_p_func, @reactor_q_func, @reactor_r_func,... 
        @ghost_left_func1, @ghost_left_func2, @ghost_left_func_rhs,...
        @ghost_right_func1, @ghost_right_func2, @ghost_right_func_rhs,...
        a,b,N);
    
    xx = z(:,1);
    yy = z(:,2);
end

%% Custom ODEBVP Solver Function
function z = ODEBVP_boundary(p_func,q_func,r_func,...
    ghost_left_func1, ghost_left_func2, ghost_left_func_rhs,...
    ghost_right_func1, ghost_right_func2, ghost_right_func_rhs,...
    a,b,N,varargin)
% A program to solve the two point boundary value problem
%   y''=p(x)y'+q(x)y+r(x),  a<x<b
%   y(a)=g1,  y(b)=g2
% 
% Rather than take a boundary value, the boundary condition may be
% formulated as an equation instead. When trying to apply the governing 
% equation above to boundary nodes, "ghost" nodes are produced. You can 
% eliminate the ghost nodes with the boundary condition function. This adds
% extra values to the other nodes for that boundary computation. For
% example the following is the general governing equation:
%   -(1+h/2p(x=0))C_-1 + (2+h^2q(x=0))C_0 + (h/2p(x=0) - 1)C_1 = -h^2r(x=0)
% And take the following as the boundary condition function
%   c'(x = 0) = 0
% With the finite difference, 
%   C_1 - C_-1 = 0, C_1 = C_-1.
% We can rewrite as C_-1 = C_1 + 0(C_0) + 0. Now we can replace the ghost
% node C_-1 as a function of the other nodes C_1, C_0, and some constant
% value. This results in three functions for the left boundary, and three
% more functions for the right boundary. In this example, ghost_left_func2
% = 1, ghost_left_func1 = ghost_left_func_rhs = 0, but they can be
% arbitrary functions. 
% 
% Lastly, ghost functions can be dependent on both x (input) and h, as well
% as any extra variables you desire (varargin).
% 
% Input
%   p, q, r: coefficient functions 
%   ghost_left_func1: left boundary function addition for the smaller node
%   ghost_left_func2: left boundary function addition for the bigger node
%   ghost_left_func_rhs: left boundary function for a constant value
%   ghost_right_func1: right boundary function addition for the smaller node
%   ghost_right_func2: right boundary function addition for the bigger node
%   ghost_right_func_rhs: right boundary function for a constant value
%   a, b: the end-points of the interval
%   N: number of sub-intervals
% Output
%   z = [ xx yy ]: xx is an (N+1) column vector of the node points
%                yy is an (N+1) column vector of the solution values
% The user must provide m-files to define the functions p, q and r.
% 
% Other MATLAB program called: tridiag.m
%
% Initialization
N1 = N+1;
h = (b-a)/N;
h2 = h*h;
xx = linspace(a,b,N1)';
yy = zeros(N1,1);

% Boundary conditions (at two endpoints)
% yy(1) = ga;
% yy(N1) = gb;

% Define the sub-diagonal avec, main diagonal bvec, superdiagonal cvec for
% the tridiagonal system
% Key point: Some of these functions may be zero-valued, but they still 
% need to be passed as vectors rather than scalars.
% pp(2:N) = p_func(xx(2:N),varargin{:});
% avec(2:N-1) = -1-(h/2)*pp(3:N);
% bvec(1:N-1) = 2+h2*q_func(xx(2:N),varargin{:});
% cvec(1:N-2) = -1+(h/2)*pp(2:N-1);

pp = p_func(xx, varargin{:});
avec(2:N+1) = -1-(h/2)*pp(2:N+1);
bvec(1:N+1) = 2+h2*q_func(xx, varargin{:});
cvec(1:N)   = -1+(h/2)*pp(1:N);

bvec(1) = bvec(1) + -(1+(h/2)*pp(1)).*ghost_left_func1(xx(1), h, varargin{:});
cvec(1) = cvec(1) + -(1+(h/2)*pp(1)).*ghost_left_func2(xx(1), h, varargin{:});

avec(N+1) = avec(N+1) + ((h/2)*pp(N+1) - 1).*ghost_right_func1(xx(N+1), h, varargin{:});
bvec(N+1) = bvec(N+1) + ((h/2)*pp(N+1) - 1).*ghost_right_func2(xx(N+1), h, varargin{:});

% Define the right hand side vector fvec
fvec = -h2*r_func(xx,varargin{:});
fvec(1)   = fvec(1)   + (1+(h/2)*pp(1)).*ghost_left_func_rhs(xx(1), h, varargin{:});
fvec(N+1) = fvec(N+1) + -((h/2)*pp(N+1) - 1).*ghost_right_func_rhs(xx(N+1), h, varargin{:});

% Solve the tridiagonal system
yy = tridiag(avec,bvec,cvec,fvec,N+1,0);
% Store output (independent variable xx and dependent variable yy)
z = [xx'; yy']';

end % end of function

%% Embedded Functions
function [x, lambda, delta, ier] = tridiag(l,d,u,b,n,iflag)

% function [x, lambda, delta, ier] = tridiag(l,d,u,b,n,iflag)
%
% Solve a tridiagonal linear system A*x=b
%
% SM: Note: This is a modified version of the code from the Atkinson
% textbook where the notation has been modified to avoid using variables 
% in conflict with those already defined for Ax = b.
%
% INPUT:
% The order of the linear system is given in n.
% The subdiagonal, diagonal, and superdiagonal of A are given
% by the arrays l,d,u, respectively.  More precisely,
%     A(i,i-1) = l(i), i=2,...,n (lower diagnonal)
%       Note: The lower diagnonal vector is assumed to be populated in
%       elements 2 through n (i.e. with the first element left as a zero).
%       So the input vector should be of size nx1 with the first element
%       zero (empty).
%     A(i,i)   = d(i), i=1,...,n (main diagonal)
%       Note: The main diagnonal vector should be of size nx1. 
%     A(i,i+1) = u(i), i=1,...,n-1 (upper diagonal)
%        Note: The upper diagnonal vector should be of size (n-1)x1. 
%
% iflag=0 means that the original matrix A is given as specified above.
% iflag=1 means that the LU factorization of A is already known and is
% stored in l,d,u.  This will have been accomplished by a previous call
% to this routine.  In that case, the vectors lambda and delta should 
% have been substituted for a and b in the calling sequence.
%
% OUTPUT:
% Upon exit, the LU factorization of A is already known and is stored
% in lambda,delta,u.  The solution x is given as well.
% ier=0 means the program was completed satisfactorily.
% ier=1 means that a zero pivot element was encountered and the 
% solution process was abandoned.

l(1) = 0;
if iflag == 0
   % Compute LU factorization of matrix A.
   for j=2:n
      if d(j-1) == 0
         ier = 1;
         return
      end
      l(j) = l(j)/d(j-1);
      d(j) = d(j) - l(j)*u(j-1);
   end
   if d(n) == 0 
      ier = 1;
      return
   end
end

% Compute solution x to A*x = b.
% Do forward substitution to solve lower triangular system.
for j=2:n
   b(j) = b(j) - l(j)*b(j-1);
end

% Do backward substitution to solve upper triangular system.
b(n) = b(n)/d(n);
for j=n-1:-1:1
   b(j) = (b(j) - u(j)*b(j+1))/d(j);
end

% Set output variables.
ier = 0;
x = b;
lambda = l; delta = d;
end % end of function
