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
## {[@var{a},@var{b},@var{c}]=}Mcapacitors@
## (@var{string}, @var{parameters}, @var{parameternames}, @
## @var{extvar},@var{intvar},@var{t})
##
## SBN file implementing models for capacitors.
##
## @var{string} is used to select among models. Parameters are listed
## as inner items. Possible models are:
##
## @enumerate
## @item @var{string} = "LIN"  (Linear Capacitor)
## @itemize @minus
## @item C -> capacitance value
## @end itemize
## @item @var{string} = "MULTICAP" (Multipole Capacitor)
## @itemize @minus
## @item C -> capacitance values
## @end itemize
## @item @var{string} = "PDE_NMOS" (Drift-Diffusion PDE NMOS capacitor)
## @itemize @minus
## @item tbulk  -> bulk thickness
## @item tox    -> oxide thickness
## @item Nnodes -> number of nodes of 1D grid 
## @item Na     -> bulk doping
## @item toll   -> absolute tolerance
## @item maxit  -> max iterations number
## @item Area   -> device area
## @end itemize
## @end enumerate
##
## @seealso{ PRSiffparse, ASMinitsystem, ASMbuildsystem, the IFF file
## format  specifications }
## @end deftypefn

function [a,b,c] = Mcapacitors(string,parameters,parameternames,extvar,intvar,t)
  
  if isempty(intvar)
    intvar = 0;
  endif

  switch string 
      ##LCR part
    case "LIN"
      for ii=1:length(parameternames)
	eval([parameternames{ii} "=" num2str(parameters(ii)) ";"])	
      endfor
            
      a = [0 0 1; 0 0 -1; 0 0 0];
      b = [0 0 0;0 0 0;-C C 1];
      c = [0 0 0]';
      break

    case "MULTICAP"
      
      n = length(extvar);
      C = reshape(parameters,n,n);
      
      a = [zeros(n) eye(n); zeros(n) zeros(n)];
      b = [zeros(n) zeros(n); -C eye(n)];
      c = [zeros(2*n,1)]';
      
      break  

      ##NLC part
    case "PDE_NMOS"
      
      constants
      
      tbulk =  1.5e-6;
      tox   =  90e-9;
      len = tbulk + tox;
      Nnodes = 300;
      Na=1e21;
      toll  = 1e-10;
      maxit = 1000;
      Area = 1e-12;

      for ii=1:length(parameternames)
	eval([parameternames{ii} "=" num2str(parameters(ii)) ";"])	
      endfor
      
      Vg = extvar(1) - extvar(2);
      q  = intvar(1);

      [Q,C]=Mnmoscap(tbulk,tox,Area,Vg,Na,Nnodes,toll,maxit);
      
      a = [0 0 1; 0 0 -1; 0 0 0];
      b = [0 0 0;0 0 0;C -C -1];
      c = [0 0 Q-q]';
      break  

    otherwise
      error (["unknown section:" string])
  endswitch

endfunction