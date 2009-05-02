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
## TSTbweulernr(@var{cirstruct}, @var{x},@var{t}, @var{tol},@
##            @var{maxit},@var{dmp}, @var{pltvars},@
##            @var{verbosity} ,@var{dae_fun});
##
## TSTbweulernr is the same as TSTbweuler except that no steady state
## simulation is computed to initialize the first timestep.
## 
## @seealso{TSTbweulernr, TSTdaspk, NLSnewtonraphson, TSTthetamethod}
##
## @end deftypefn

function [out, varargout] = TSTbweulernr(outstruct, out, t, tol, maxit, dmp, pltvars, verbosity, istart, dae_fun)

  if ~exist("verbosity")
    verbosity = [0,0];
  elseif length(verbosity)<2
    verbosity(2) =0;
  endif
  
  if nargout > 1
    niter = zeros(length(t),1);
  endif
  
  x = out(:,istart);
  [A0,B,C,outstruct] = ASMinitsystem(outstruct,x,t(1));
  
  if (nargin > 9)
    JAC = @(x) dae_fun{1}(outstruct,x,t(1),B);
    RES = @(x) dae_fun{2}(outstruct,x,t(1),B,C);
  else
    JAC = @(x) TSTBWEFUNJAC0(outstruct,x,t(1),B);
    RES = @(x) TSTBWEFUNRES0(outstruct,x,t(1),B,C);
  endif
    
  for it=istart+1:length(t)
    
    if (verbosity)
      fprintf(1,"timestep #%d:\n",it);
    endif

    if nargin > 9
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
					     tol,  maxit, verbosity(1), UPDT);

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