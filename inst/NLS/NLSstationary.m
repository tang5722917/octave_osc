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
##  @deftypefn{Function File} @var{out} = NLSstationary @
##  (@var{instruct},@var{x},@var{tol},@var{maxit},@var{dmp})
##  Computes a stationary state for the system described by
##  @var{sysstruct}.
##
##  Input:
##  @itemize @minus
##  @item @var{instruct}: system description structure
##  @item @var{x}: initial guess 
##  @item @var{tol, maxit, dmp}: parameters to be passed to NLSnewtonraphson
##  @end itemize
##  @end deftypefn

function out = NLSstationary(outstruct,x,tol,maxit,dmp)

  [A0,B,C,outstruct] = ASMinitsystem(outstruct,x,0);
  JAC = @(x) TSTBWEFUNJAC0(outstruct,x,0,B);
  RES = @(x) TSTBWEFUNRES0(outstruct,x,0,B,C);
  [out,ii,resnrm] = NLSnewtonraphson(x, RES, JAC, tol, maxit);

endfunction