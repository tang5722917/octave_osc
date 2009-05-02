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
## @deftypefn{Function File} @
## {[@var{lhs}]=} TSTDASPKFUNJAC@
## (@var{outstruct}, @var{x}, @var{xdot}, @var{A0}, @var{B}, @var{t},
## @var{c})
##
## INTERNAL FUNCTION:
##
## NOT SUPPOSED TO BE CALLED DIRECTLY BY USERS
## @end deftypefn

function lhs = TSTDASPKFUNJAC(outstruct,x,xdot,A0,B,t,c)

  [A1,Jac,res] = ASMbuildsystem(outstruct,x,t);
  lhs = ( c*(A0+A1) + B + Jac); 

endfunction