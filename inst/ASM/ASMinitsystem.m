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
## author: Carlo de Falco <cdf _AT_ users.sourceforge.net> 
## author: culpo@math.uni-wuppertal.de

## -*- texinfo -*-
##
## @deftypefn{Function File} {[@var{A}, @var{B}, @var{C},@
## @var{struct}] =} ASMinitsystem (@var{instruct}, @var{x}) 
##
## Cycles through the circuit description structure @var{instruct} 
## to build the system matrices @var{A}, @var{B}, @var{C} for
## the linear and time-invariant part of the system
##
## @itemize @minus
## @item @var{x} is the current value of the state variables
## @end itemize
##
## see the @cite{IFF file format specifications} for details about 
## the output matrices.
## 
## @seealso{ASMbuildsystem}
##
## @end deftypefn

function  [varargout] = ASMinitsystem(instruct,x);

  if (~isfield(instruct,"totintvar"))
    ## Check number of internal variables
    intvar = 0;

    ## NLC section
    nblocks = length(instruct.NLC);
    
    for ibl = 1:nblocks
      for iel = 1:instruct.NLC(ibl).nrows

	## evaluate element
	il = instruct.NLC(ibl).vnmatrix(iel,:)';
	nzil = find(il!=0);
      
	y = zeros(size(il));
	y(nzil)=x(il(nzil));
      
	[a,b,c] = feval(instruct.NLC(ibl).func,...
			instruct.NLC(ibl).section,...
			instruct.NLC(ibl).pvmatrix(iel,:),...
			instruct.NLC(ibl).parnames,...
			y,[],0);

	instruct.NLC(ibl).nintvar(iel) = columns(a)-instruct.NLC(ibl).nextvar;
	instruct.NLC(ibl).osintvar(iel) = intvar;
	intvar += instruct.NLC(ibl).nintvar(iel);

      endfor
    endfor


    ## LCR section
    nblocks = length(instruct.LCR);

    for ibl = 1:nblocks
      for iel = 1:instruct.LCR(ibl).nrows
	
	## evaluate element
	il = instruct.LCR(ibl).vnmatrix(iel,:)';
	nzil = find(il!=0);
	
	y = zeros(size(il));
	y(nzil)=x(il(nzil));
	
      
	[a,b,c] = feval(instruct.LCR(ibl).func,...
			instruct.LCR(ibl).section,...
			instruct.LCR(ibl).pvmatrix(iel,:),...
			instruct.LCR(ibl).parnames,...
			y,[],0);

	instruct.LCR(ibl).nintvar(iel) = columns(a)-instruct.LCR(ibl).nextvar;
	instruct.LCR(ibl).osintvar(iel) = intvar;
	intvar += instruct.LCR(ibl).nintvar(iel);
	
      endfor
    endfor
    instruct.totintvar = intvar;
  endif

  ## Build linear part of the system
  n = instruct.totextvar+instruct.totintvar;
  lx = length(x);
  if lx < n
    x(lx+1:n) = 0;
  endif
  A = spalloc(n,n,0);
  
  ## LCR section
  B = spalloc(n,n,0);
  C = spalloc(n,1,0);

  nblocks = length(instruct.LCR);

  for ibl = 1:nblocks
    for iel = 1:instruct.LCR(ibl).nrows
      
      ## evaluate element
      if instruct.LCR(ibl).nintvar(iel)
	intvars = instruct.totextvar+instruct.LCR(ibl).osintvar(iel)+...
	    [1:instruct.LCR(ibl).nintvar(iel)]';
      else
	intvars=[];
      endif

      il = instruct.LCR(ibl).vnmatrix(iel,:)';
      nzil = find(il!=0);
      
      y = zeros(size(il));
      y(nzil)=x(il(nzil));
      z = x(intvars);
      
      [a,b,c] = feval(instruct.LCR(ibl).func,...
		      instruct.LCR(ibl).section,...
		      instruct.LCR(ibl).pvmatrix(iel,:),...
		      instruct.LCR(ibl).parnames,...
		      y,z,0);
      
      ## ASSEMBLE MATRICES
      
      ## global indexing
      vars    = [il(nzil);intvars];
      ## local indexing
      lclvars = [nzil; instruct.LCR(ibl).nextvar + (1:length(intvars))' ];
      ## reshaping sparse stamps
      a = a(lclvars,lclvars);
      b = b(lclvars,lclvars);
      c = reshape(c(lclvars),[],1);
      
      [na,ma,va] = find(a);
      [nb,mb,vb] = find(b);
      [nc,mc,vc] = find(c);

      ## stamping
      A   += sparse(vars(na),vars(ma),va,n,n);
      B   += sparse(vars(nb),vars(mb),vb,n,n);
      C   += sparse(vars(nc),1,vc,n,1);

    endfor
  endfor

  varargout{1}=A;
  varargout{2}=B;
  varargout{3}=C;
  varargout{4}=instruct;
endfunction