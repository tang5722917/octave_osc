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
## author: Carlo de Falco <cdf _AT_ users.sourceforge.net> 

## -*- texinfo -*-
##
## @deftypefn{Function File} {[@var{out}, @var{niter}] =}  @
## TSTbweuler(@var{cirstruct}, @var{x},@var{t}, @var{tol},@
##            @var{maxit},@var{dmp}, @var{pltvars},@
##            @var{verbosity} ,@var{dae_fun});
##
## Performs a transient simulation of the system described by
## @var{cirstruct} over the time interval @var{t} using the backward
## Euler algorithm
##
## An initial value for the state vector is computed by solving a
## steady  state problem at @var{t}(1)
## the initial guess for the state vector is set to  @var{x} 
## @var{tol} @var{maxit} and @var{dmp} are parameters passed to
## NLSnewtonraphson.
##
## The output @var{out} will contain the value of the state vector at
## each  point of @var{t}.
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
## For special purposes one may need to pass modified jacobian and
## residual functions, this can be done
## via the cell array of function handles @var{dae_fun}, such
## functions  should have the same input and output
## parameter list as the default functions
## TSTBWEFUNJAC0,TSTBWEFUNRES0, TSTBWEFUNJAC,TSTBWEFUNRES
##
## The optional output @var{niter} returns the number of Newton iterations
## needed to reach convergence.
## 
## 
## @seealso{TSTdaspk, NLSnewtonraphson, TSTthetamethod}
##
## @end deftypefn

function [out, varargout] = TSTbweuler(outstruct,x,t,tol,maxit,dmp,pltvars,verbosity,dae_fun)

  if ~exist("verbosity")
    verbosity = [0,0];
  elseif length(verbosity)<2
    verbosity(2) =0;
  endif
  
  out      = zeros(rows(x),columns(t));
  out(:,1) = x;
  
  if nargout > 1
    niter = zeros(length(t),1);
  endif

  if (verbosity(1))
    fprintf(1,"initial value:\n");
  endif
  
  
  [A0,B,C,outstruct] = ASMinitsystem(outstruct,x,t(1));
  
  if (nargin > 8)
    JAC = @(x) dae_fun{1}(outstruct,x,t(1),B);
    RES = @(x) dae_fun{2}(outstruct,x,t(1),B,C);
  else
    JAC = @(x) TSTBWEFUNJAC0(outstruct,x,t(1),B);
    RES = @(x) TSTBWEFUNRES0(outstruct,x,t(1),B,C);
  endif
  
  %%out = repmat (x, 1, length(t));
  [out(:,1),ii,resnrm] = NLSnewtonraphson(x, RES, JAC, tol, maxit,verbosity(1));
  %%nrm(1) = resnrm(ii);
  
  for it=2:length(t)
    
    if (verbosity)
      fprintf(1,"timestep #%d:\n",it);
    endif

    if nargin > 8
      JAC = @(x) dae_fun{3}(outstruct,x,t(it-1),t(it),A0,B);
      RES = @(x) dae_fun{4}(outstruct,x,out(:,it-1),t(it-1),t(it),A0,B,C);
    else
      JAC = @(x,A1,Jac,res) TSTBWEFUNJAC1(outstruct, x, t(it-1), 
					  t(it), A0, B, A1, Jac, res);
      RES = @(x,A1,Jac,res) TSTBWEFUNRES1(outstruct, x, out(:,it-1), 
					  t(it-1), t(it), A0, B, C, 
					  A1, Jac, res);
      UPDT = @(x) TSTBWEFUNUP1 (outstruct, x, t(it));
    endif

    [out(:,it),ii,resnrm] = NLSnewtonraphson(out(:,it-1), RES, JAC, ...
					     tol,  maxit, verbosity(1), 
					     UPDT);
    %%nrm(it) = resnrm(ii);

    if nargout > 1
      niter(it) = ii;
    endif
    
    if (verbosity(2))
     UTLplotbyname(t(1:it),out(:,1:it),outstruct,pltvars), pause(.01)
    endif
  
    if exist("~/.stop_ocs","file")
      printf("stopping at timestep %d\n",it);
      unix("rm ~/.stop_ocs");
      break
    end

  endfor

  if nargout > 1
    varargout{1} = niter;
  endif

endfunction