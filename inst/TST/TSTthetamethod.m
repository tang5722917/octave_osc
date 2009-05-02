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
## @deftypefn{Function File} {[@var{out}] =}  TSTthetamethod @
## (@var{cirstruct}, @var{x},@var{t}, @var{tol},@
## @var{maxit},@var{dmp}, @var{theta}, @var{pltvars},@
## @var{verbosity});
##
## Performs a transient simulation of the system described by
## @var{cirstruct} over the time interval @var{t} using the
## theta-method with parameter @var{theta}
##
## An initial value for the state vector is computed by solving a
## steady  state problem at @var{t}(1)
## the initial guess for the state vector is set to  @var{x} 
## @var{tol} @var{maxit} and @var{dmp} are parameters passed to
## NLSnewtonraphson.
##
## The output @var{out} will contain the value of the state vector at
## each  point of @var{t}
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
## 
## @seealso{TSTbweuler, NLSnewtonraphson, daspk}
##
## @end deftypefn

function out = TSTthetamethod(outstruct, x, t, tol, maxit, dmp, theta,
			      pltvars, verbosity)

  if ~exist("verbosity")
    verbosity = [0,0];
  elseif length(verbosity)<2
    verbosity(2) = 0;
  endif

  out=zeros(rows(x),columns(t));
  out(:,1) = x;
  
  if (verbosity(1))
    fprintf(1,"initial value:\n");
  endif
  
  
  [A0,B,C,outstruct] = ASMinitsystem(outstruct,x,t(1));
  
  JAC = @(x) TSTBWEFUNJAC0(outstruct,x,t(1),B);
  RES = @(x) TSTBWEFUNRES0(outstruct,x,t(1),B,C);
  
  [out(:,1),ii,resnrm] = NLSnewtonraphson(x, RES, JAC, tol, maxit,verbosity(1));
  nrm(1) = resnrm(ii);
  
  for it=2:length(t)
    
    if(verbosity)
      fprintf(1,"timestep #%d:\n",it);
    endif


    [A1old,Jacold,resold] = ASMbuildsystem(outstruct, out(:,it-1), t(it-1));
    
    JAC = @(x,A1,Jac,res) TSTTHETAFUNJAC1(outstruct, x, t(it-1), 
					  t(it), A0, B, theta, 
					  A1, Jac, res);
    RES = @(x,A1,Jac,res) TSTTHETAFUNRES1(outstruct, x, out(:,it-1), 
					  t(it-1), t(it), A0, B, C, 
					  resold, theta, A1, Jac, res);
    UPDT = @(x) TSTBWEFUNUP1 (outstruct, x, t(it));
    
    [out(:,it),ii,resnrm] = NLSnewtonraphson(out(:,it-1), RES, JAC, 
					     tol,  maxit, verbosity(1), 
					     UPDT);
    nrm(it) = resnrm(ii);
    
    if (verbosity(2))
      UTLplotbyname(t(1:it),out(:,1:it),outstruct,pltvars), pause(.1)
    endif
    
    if exist("~/.stop_ocs","file")
      break
    end
  endfor

endfunction