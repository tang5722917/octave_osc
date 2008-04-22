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
## {[@var{Q},@var{C}]=}Mnmoscap@
## (@var{tbulk}, @var{tox}, @var{Area}, @
## @var{Vg},@var{Na},@var{Nnodes},@var{toll},@var{maxit})
##
## INTERNAL FUNCTION:
##
## NOT SUPPOSED TO BE CALLED DIRECTLY BY USERS
## @end deftypefn

function [Q,C]=Mnmoscap(tbulk,tox,Area,Vg,Na,Nnodes,toll,maxit);
  
  constants 

  Nelements = Nnodes - 1;
  len = tox+tbulk;
  x = linspace(0,len,Nnodes)';
  sinodes = find(x<=tbulk);
  Nsinodes = length(sinodes);
  NelementsSi = Nsinodes-1;
  D = - Na* ones(Nsinodes,1);
  pp = Na ;
  p = pp* ones(Nsinodes,1);
  n = (ni^2)./p;
  Fn = 0*n;
  Fp = 0*n;
  

  V = -Phims + Vg * ones(Nnodes,1);
  V(sinodes) = Fn + Vth*log(n/ni);
  
  ## Scaling
  xs  = len;
  ns  = norm(D,inf);
  Din = D/ns;
  Vs  = Vth;
  xin   = x/xs;
  nin   = n/ns;
  pin   = p/ns;
  Vin   = V/Vs;
  Fnin  = (Fn - Vs * log(ni/ns))/Vs;
  Fpin  = (Fp + Vs * log(ni/ns))/Vs;
  
  l2    = (Vs*esio2)/(q*ns*xs^2)* ones(Nelements,1);
  l2(1:NelementsSi)    = (Vs*esi)/(q*ns*xs^2);
  
  ## Solution of Nonlinear Poisson equation
  [V,nout,pout,res,niter] = DDGnlpoisson (xin,sinodes,Vin,nin,...
				       pin,Fnin,Fpin,Din,l2,...
				       toll,maxit,0);
    
  L = Ucomplap(xin,Nnodes,[],Nelements,l2);
  C22 = L(end,end);
  C12 = L(2:end-1,end);
  C11 = L(2:end-1,2:end-1);

  drdv  = zeros(Nnodes,1);    drdv(sinodes) = nout + pout;
  coeff = zeros(Nelements,1); coeff(1:NelementsSi) = 1;
  M     = Ucompmass(xin,Nnodes,[],[],drdv,coeff);
  C     = C22 - C12'*((C11+M(2:end-1,2:end-1))\C12);
  Q     =(C12'*V(2:end-1)+C22*V(end));

  ## Descaling
  C = Area*C*(q*ns*xs/Vs);
  Q = Area*Q*(q*ns*xs);

endfunction