function [X,Y] = ode4xy(odefun,tspan,x0,varargin)
%ODE4XY Solves differential equations with a non-adaptive method of order 4.
%
% Based on the original ode4 function available at:
% https://www.mathworks.com/matlabcentral/answers/uploaded_files/5693/ODE_Solvers.zip
% https://www.mathworks.com/matlabcentral/answers/98293-is-there-a-fixed-step-ordinary-differential-equation-ode-solver-in-matlab-8-0-r2012b
%
% Modified by: Antonio Bernardo Guimarães Neto (ITA, São José dos Campos, Brazil)
% Date: 2020-03-10
%   -> The modification of the original ode4 function enables the calculation of
%      an output matrix Y containing outputs of the ODEFUN function
%
%   [X,Y] = ODE4XY(ODEFUN,TSPAN,X0) with TSPAN = [T1, T2, T3, ... TN] integrates 
%   the system of differential equations x' = f(t,x), y = g(t,x) by stepping from 
%   T0 to T1 to TN. Function ODEFUN(T,X) must return f(t,x) and g(t,x) in column 
%   vectors. The vector X0 contains the initial conditions at T0. Each row in the
%   solution arrays X and Y corresponds to a time specified in TSPAN.
%
%   [X,Y] = ODE4XY(ODEFUN,TSPAN,X0,P1,P2...) passes the additional parameters 
%   P1,P2... to the derivative function as ODEFUN(T,X,P1,P2...). 
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN
%   but the derivative function ODEFUN is evaluated multiple times per step.
%
%   The solver implements the classical Runge-Kutta method of order 4.   
%

if ~isnumeric(tspan)
  error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(x0)
  error('X0 should be a vector of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
  error('Entries of TSPAN are not in order.') 
end  

try
  [f0,y0] = feval(odefun,tspan(1),x0,varargin{:});
catch
  msg = ['Unable to evaluate the ODEFUN at t0,x0. ',lasterr];
  error(msg);  
end  

x0 = x0(:);   % Make a column vector.
if ~isequal(size(x0),size(f0))
  error('Inconsistent sizes of x0 and f(t0,x0).');
end  

neq = length(x0);
nout = length(y0);
N = length(tspan);
X = zeros(neq,N);
Y = zeros(nout,N);

X(:,1) = x0;
Y(:,1) = y0;
for i = 2:N
  ti = tspan(i-1);
  hi = h(i-1);
  xi = X(:,i-1);
  [f1,y1] = feval(odefun, ti,        xi,           varargin{:});
  [f2,y2] = feval(odefun, ti+0.5*hi, xi+0.5*hi*f1, varargin{:});
  [f3,y3] = feval(odefun, ti+0.5*hi, xi+0.5*hi*f2, varargin{:});  
  [f4,y4] = feval(odefun, tspan(i),  xi+hi*f3,     varargin{:});
  X(:,i) = xi + (hi/6)*(f1 + 2*f2 + 2*f3 + f4);
  Y(:,i) = (1/6)*(y1 + 2*y2 + 2*y3 + y4);
end
X = X.';
Y = Y.';
