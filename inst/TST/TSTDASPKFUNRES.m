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
## @deftypefn{Function File} @
## {[@var{rhs}]=} TSTDASPKFUNJAC@
## (@var{outstruct}, @var{x}, @var{xdot}, @var{A0}, @var{B}, @var{C}, @var{t})
##
## INTERNAL FUNCTION:
##
## NOT SUPPOSED TO BE CALLED DIRECTLY BY USERS
## @end deftypefn

function rhs = TSTDASPKFUNRES(outstruct,x,xdot,A0,B,C,t)
  
  [A1,Jac,res] = ASMbuildsystem(outstruct,x,t);
  rhs = (A0+A1)*xdot + B*x + C + res; 

endfunction

