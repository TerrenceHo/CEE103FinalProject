
function v = reactor_solver(D, U, k, L, cin, dx)
% A function to solve the steady state reactor problem
% Inputs:
%   D: dispersion coefficient 
%   U: mean velocity through tank
%   k: reaction rate
%   L: length of the tank
%   cin: concentration in the inflow
%   dx: discretization of the length of the reactor
% Outputs:
%   z: solution output

    a = 0;
    b = L;
    function p = reactor_p_func(x)
    p = U./D;
    end
    function q = reactor_q_func(x)
        q = k./D;
    end
    function r = reactor_r_func(x)
        r = 0.*x;
    end

    z = ODEBVP(@reactor_p_func, @reactor_q_func, @reactor_r_func, a, b,...
            ga, gb, N);
end



%% Embedded Functions
function z = ODEBVP(p_func,q_func,r_func,a,b,ga,gb,N,varargin)
% A program to solve the two point boundary value problem
%   y''=p(x)y'+q(x)y+r(x),  a<x<b
%   y(a)=g1,  y(b)=g2
% Input
%   p, q, r: coefficient functions 
%   a, b: the end-points of the interval
%   ga, gb: the prescribed function values at the end-points
%   N: number of sub-intervals
% Output
%   z = [ xx yy ]: xx is an (N+1) column vector of the node points
%                yy is an (N+1) column vector of the solution values
% A sample call would be
%   z=ODEBVP('p','q','r',a,b,ga,gb,100)
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
yy(1) = ga;
yy(N1) = gb;

% Define the sub-diagonal avec, main diagonal bvec, superdiagonal cvec for
% the tridiagonal system
% Key point: Some of these functions may be zero-valued, but they still 
% need to be passed as vectors rather than scalars.
pp(2:N) = p_func(xx(2:N),varargin{:});

avec(2:N-1) = -1-(h/2)*pp(3:N);
bvec(1:N-1) = 2+h2*q_func(xx(2:N),varargin{:});
cvec(1:N-2) = -1+(h/2)*pp(2:N-1);

% Define the right hand side vector fvec
fvec(1:N-1) = -h2*r_func(xx(2:N),varargin{:});

fvec(1) = fvec(1)+(1+h*pp(2)/2)*ga;
fvec(N-1) = fvec(N-1)+(1-h*pp(N)/2)*gb;

% Solve the tridiagonal system
yy(2:N) = tridiag(avec,bvec,cvec,fvec,N-1,0);
% Store output (independent variable xx and dependent variable yy)
z = [xx'; yy']';

end % end of function
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
