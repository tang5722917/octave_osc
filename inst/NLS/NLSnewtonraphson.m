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
##  @deftypefn{Function File}{[@var{y},@var{numit},@var{resnrm}] =} @
## NLSnewtonraphson @
## (@var{y0}, @var{RES}, @var{JAC}, @var{tol}, @
##  @var{maxit},@var{verbosity}, @var{update});
##
##  Solves a non-linear system of equations using the Newton-Raphson
##  method with damping.
##
##  Input:
##  @itemize @minus 
##  @item @var{y0}: initial guess
##  @item @var{RES}: function handle to compute the residual
##  @item @var{JAC} function handle to compute the jacobian
##  @item @var{tol}: tolerance for convergence check
##  @item @var{maxit}: maximum number of iterations
##  @item @var{verbosity}: verbosity level
##  @item @var{update}: update function to run at each timestep
##  @end itemize
##
##  Output:
##  @itemize @minus 
##  @item @var{y}: initial guess
##  @item @var{numit}: number of iterations performed
##  @item @var{resnrm}: residual norm at each step
##  @end itemize
##  @seealso{NLSstationary, TSTbweuler, TSTdaspk}
##  @end deftypefn 

function [y,ii,resnrm] = NLSnewtonraphson (y0,RES,JAC,tol,maxit,verbosity,update);

  if ~exist("verbosity")
    verbosity = 0;
  endif

  if ~exist("update")
    update = @(x) ({});
  endif
  
  jjtot = 0;
  y = y0;
  
  uptodate = update(y);
  res_y = RES(y,uptodate{:});
  resnrm(1) = norm(res_y,inf);
  
  for ii=1:maxit
    
    jac_y      = JAC(y,uptodate{:}); 
    ynew       = jac_y\(-res_y+jac_y*y);
    uptodate   = update(ynew);
    res_y      = RES(ynew,uptodate{:}); 
    resnrm(ii+1) = norm(res_y,inf);
    
    jj=0;
    while (resnrm(ii+1)>resnrm(ii))&(jj<10)
      jj++;
      damp=2^(-jj);
      ynew = y*(1-damp)+ynew*damp;
      uptodate = update(ynew);
      res_y  = RES(ynew,uptodate{:});
      resnrm(ii+1) = norm(res_y,inf);
    endwhile
    
    jjtot+=jj;
    y = ynew;
    
    if resnrm(ii+1)<tol 
      if (verbosity)
	fprintf(1,"converged in %d newton iterations and ",ii);
	fprintf(1,"%d damping iterations\n",jjtot);
      endif
      break
    elseif ii==maxit
      if(verbosity)
	fprintf(1,"not converged, nrm=%g\n",resnrm(maxit))
      endif
      break
    endif
  endfor
  
endfunction