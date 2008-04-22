## Copyright (C) 2006,2007,2008  Carlo de Falco            
##
## This file is part of:
## OCS - A Circuit Simulator for Octave
##
## OCS is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program (see the file LICENSE); if not,
## see <http://www.gnu.org/licenses/>.
##
## author: carlo.defalco@gmail.com 

## -*- texinfo -*-
##
## @deftypefn{Function File} {[@var{out}] =}  TSTdaspk @
## (@var{cirstruct}, @var{x},@var{t},  @var{tol}, @var{maxit},@
## @var{dmp}, @var{pltvars},@var{verbosity}, @var{daskopts},@var{dae_fun});
##
## Performs a transient simulation of the system described by
## @var{cirstruct}  over the time interval @var{t} using daspk.
##
## An initial value for the state vector is computed by solving
## a steady state problem at @var{t}(1).
## The initial guess for the state vector is set to  @var{x}.   
## @var{tol} @var{maxit} and @var{dmp} are parameters passed
## to NLSnewtonraphson.
##
## The output @var{out} will contain the value of the state vector
## at each point of @var{t}
##
## Extra options for daspkcan be passed as name/value pairs in
## the cellarray @var{daskopts}.
##
## The optional parameter  @var{verbosity} controls the amount of
## output produced: 
## 
## @itemize @minus 
## @item if verbosity(1) != 0, information on the progress
## of the algorithm are output at runtime
## @item if verbosity(2) != 0, the plot of the variables whose names
## are listed in @var{pltvars} is
## produced after the computation
## @end itemize
##
## For special purposes one may need to pass modified jacobian
## and residual functions, this can be done
## via the cell array of function handles @var{dae_fun}, such
## functions should have the same input and output
## parameter list as the default functions
## TSTBWEFUNJAC0,TSTBWEFUNRES0, TSTDASPKFUNJAC,TSTDASPKFUNRES
## @seealso{TSTbweuler, NLSnewtonraphson, daspk}
##
## @end deftypefn

function out = TSTdaspk (outstruct, x, t, tol, maxit, 
			 dmp, pltvars, verbosity, daspkopts, dae_fun)
  if ~exist("verbosity")
    verbosity = [0,0];
  elseif length(verbosity)<2
    verbosity(2) =0;
  endif

  if(verbosity(1))
    fprintf(1,"initial value:\n");
  endif
  
  daspk_options ("print initial condition info",1);
  daspk_options("maximum order",2);
  daspk_options("initial step size",t(2)-t(1));
  daspk_options("relative tolerance",1e-3);

  if ( nargin > 8 )
    for ii = 1:2:length(daspkopts)
      daspk_options (daspkopts{ii},daspkopts{ii+1});
    endfor
  endif
  

  [A0,B,C,outstruct] = ASMinitsystem(outstruct,x,t(1));

  if nargin > 9
    JAC = @(x) dae_fun{1}(outstruct,x,t(1),B);
    RES = @(x) dae_fun{2}(outstruct,x,t(1),B,C);
  else
    JAC = @(x) TSTBWEFUNJAC0(outstruct,x,t(1),B);
    RES = @(x) TSTBWEFUNRES0(outstruct,x,t(1),B,C);
  endif

  [x,ii,resnrm] = NLSnewtonraphson(x, RES, JAC, tol, maxit,verbosity);
  nrm(1) = resnrm(ii);
  
  if nargin > 9
    JAC = @(x) dae_fun{3}(outstruct,x,t(1),B);
    RES = @(x) dae_fun{4}(outstruct,x,t(1),B,C);
  else  
    JAC = @(x,xdot,t,c) TSTDASPKFUNJAC(outstruct,x,xdot,A0,B,t,c);
    RES = @(x,xdot,t) TSTDASPKFUNRES(outstruct,x,xdot,A0,B,C,t);
  endif

  [out, xdot, istate, msg] = daspk ({RES,JAC}, x, zeros(size(x)), t);

  out = out';
  if verbosity(2)
    UTLplotbyname(t,out,outstruct,pltvars)
  endif

endfunction