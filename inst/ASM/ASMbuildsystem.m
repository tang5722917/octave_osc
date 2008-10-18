## Copyright (C) 2006,2007,2008  Carlo de Falco, Culpo Massimiliano
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
## author: culpo@math.uni-wuppertal.de

## -*- texinfo -*-
## @deftypefn{Function File} [@var{A}, @var{Jac}, @var{res} ] = @
## ASMbuildsystem (@var{instruct}, @var{x}, @var{t})
##
## Cycles through the circuit description structure @var{instruct} 
## to build the system matrices @var{A}, @var{Jac}, @var{res} for
## the current step of the Newton method.
##
## @itemize @minus
## @item @var{x} is the current value of the state variables
## @item @var{t} is the current time point
## @end itemize
##
## see the @cite{IFF file format specifications} for details about 
## the output matrices.
## 
## @seealso{ASMinitsystem}
##
## @end deftypefn

function  [A,Jac,res] = ASMbuildsystem(instruct,x,t);

  n = instruct.totextvar+instruct.totintvar;
  A = spalloc(n,n,0);
  Jac = spalloc(n,n,0);
  res = spalloc(n,1,0);
  
  
  ## NLC section
  nblocks = length(instruct.NLC);

  for ibl = 1:nblocks
    for iel = 1:instruct.NLC(ibl).nrows

      ## evaluate element
      if instruct.NLC(ibl).nintvar(iel)    
	intvars = instruct.totextvar+instruct.NLC(ibl).osintvar(iel)+...
	    [1:instruct.NLC(ibl).nintvar(iel)]';
      else
	intvars=[];
      endif
    
      il = instruct.NLC(ibl).vnmatrix(iel,:)';
      nzil = find(il!=0);
      
      y = zeros(size(il));
      y(nzil)=x(il(nzil));
      
      z = x(intvars);
      
      [a,b,c] = feval(instruct.NLC(ibl).func,...
		      instruct.NLC(ibl).section,...
		      instruct.NLC(ibl).pvmatrix(iel,:),...
		      instruct.NLC(ibl).parnames,...
		      y,z,t);
      
      ## ASSEMBLE MATRICES
      
      ## global indexing
      vars    = [il(nzil);intvars];
      ## local indexing
      lclvars = [nzil; instruct.NLC(ibl).nextvar + (1:length(intvars))' ];
      ## reshaping sparse stamps
      [atmp,btmp,ctmp,mm,nn] = shape(a,b,c,vars,lclvars);
      ## stamping
      A   += sparse(nn,mm,atmp,n,n);
      Jac += sparse(nn,mm,btmp,n,n);
      res += sparse(vars,1,ctmp,n,1);
      
    endfor	
  endfor
  
  
endfunction

function [aout,bout,cout,mm,nn] = shape(ain,bin,cin,vars,lclvars)
## reshapes input in a form suitable for "sparse" function
 
  [mm, nn] = meshgrid(vars);
  nn       = reshape(nn,[],1);
  mm       = reshape(mm,[],1);
  aout     = reshape(ain(lclvars,lclvars),[],1);
  bout     = reshape(bin(lclvars,lclvars),[],1);
  cout     = reshape(cin(lclvars),[],1);
  
endfunction