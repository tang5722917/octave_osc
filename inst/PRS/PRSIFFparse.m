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
## @deftypefn{Function File} @var{outstruct} = PRSIFFparse(@var{name})
## Parses a netlist in IFF format and produces the system description
## structure @var{outstruct}.
## @var{name} is the basename of the CIR and NMS files to
## be parsed.
##
## See the @cite{IFF file format specifications} for more details
## @end deftypefn

function outstruct = PRSIFFparse(name)
  
  ## init
  version ="0.1b1";
  outstruct = struct("NLC",[],...
                     "LCR",[],...
                     "totextvar",0);
  
  ## open cir file
  filename = [name ".cir"];
  if isempty(file_in_path(".",filename))
    error([".cir file not found:" filename]);
  endif
  fid = fopen(filename,"r");

  ## Check version
  line = fgetl(fid);
  
  if line(1)~="%"
    error(["missing version number in file " filename]);
  endif
  
  if ~strcmp(version,sscanf(line(2:end),"%s"));
    error(["conflicting version number in file " filename]);
  endif

  ## NLC section
  NLCcount = 0;
  while ~strcmp(line,"END")

    ## skip  comments
    while line(1)=="%"
      line = fgetl(fid);
    endwhile

    if strcmp(line,"END")
      break
    else
      NLCcount++;
    endif
    
    ## parse block header
    [outstruct]=parseNLCblockheader(fid,line,outstruct,NLCcount);

    ## parse block par-value matrix
    [outstruct.NLC(NLCcount).pvmatrix]=...
	fscanf(fid,"%g",[outstruct.NLC(NLCcount).npar,...
			 outstruct.NLC(NLCcount).nrows])';

    ## parse block var-number matrix
    [outstruct.NLC(NLCcount).vnmatrix]=...
	fscanf(fid,"%g",[outstruct.NLC(NLCcount).nextvar,...
			 outstruct.NLC(NLCcount).nrows])';
    
    outstruct.totextvar = max([max(outstruct.NLC(NLCcount).vnmatrix(:)) 
                               outstruct.totextvar]);

    ## skip the newline char after the matrix
    line = fgetl(fid);
    
    ## proceed to next line
    line = fgetl(fid);

  endwhile

  ## LCR section
  LCRcount = 0;
  line = fgetl(fid);

  while (~strcmp(line,"END"))

    ## skip  comments
    while line(1)=="%"
      line = fgetl(fid);
    endwhile

    if strcmp(line,"END")
      break
    else
      LCRcount++;
    endif
    
    ## parse block header
    [outstruct]=parseLCRblockheader(fid,line,outstruct,LCRcount);
    
    ## parse block par-value matrix
    [outstruct.LCR(LCRcount).pvmatrix]=...
	fscanf(fid,"%g",[outstruct.LCR(LCRcount).npar,...
			 outstruct.LCR(LCRcount).nrows])';
    
    ## parse block var-number matrix
    [outstruct.LCR(LCRcount).vnmatrix]=...
	fscanf(fid,"%g",[outstruct.LCR(LCRcount).nextvar,...
			 outstruct.LCR(LCRcount).nrows])';

    outstruct.totextvar = max([max(outstruct.LCR(LCRcount).vnmatrix(:)) 
                               outstruct.totextvar]);
    
    ## skip the newline char after the matrix
    line = fgetl(fid);
    
    ## proceed to next line
    line = fgetl(fid);

  endwhile

  ## fclose cir file
  fclose(fid); 

  ## open nms file
  filename = [name ".nms"];
  if isempty(file_in_path(".",filename))
    error([".nms file not found:" filename]);
  endif
  fid = fopen(filename,"r");

  ## Check version
  line = fgetl(fid);
  
  if line(1)~="%"
    error(["missing version number in file " filename]);
  endif
  
  if ~strcmp(version,sscanf(line(2:end),"%s"));
    error(["conflicting version number in file " filename]);
  endif

  ## Init
  cnt = 1;
  outstruct.namesn = [];
  outstruct.namess = {};
  nnames = 0;
  
  while cnt
    [nn,cnt] = fscanf(fid,"%d","C");
    [ns,cnt] = fscanf(fid,"%s","C");
    if cnt
      outstruct.namesn(++nnames)=nn;
      outstruct.namess(nnames)=ns;
    endif
  endwhile
  
  ## fclose nms file
  fclose(fid);

endfunction


##############################################
function [outstruct]=parseNLCblockheader(fid,line,outstruct,NLCcount);

  [func,section,nextvar,npar]=sscanf(line,"%s %s %g %g","C");
  outstruct.NLC(NLCcount).func = func;
  outstruct.NLC(NLCcount).section = section;
  outstruct.NLC(NLCcount).nextvar = nextvar;
  outstruct.NLC(NLCcount).npar = npar;
  [nrows,nparnames]=fscanf(fid,"%g %g","C");
  outstruct.NLC(NLCcount).nrows = nrows;
  outstruct.NLC(NLCcount).nparnames = nparnames;
  outstruct.NLC(NLCcount).parnames = {};
  for ii=1:nparnames
    outstruct.NLC(NLCcount).parnames{ii}=fscanf(fid,"%s","C");
  endfor

endfunction

##############################################
function     [outstruct]=parseLCRblockheader(fid,line,outstruct,LCRcount);

  [func,section,nextvar,npar]=sscanf(line,"%s %s %g %g","C");
  outstruct.LCR(LCRcount).func = func;
  outstruct.LCR(LCRcount).section = section;
  outstruct.LCR(LCRcount).nextvar = nextvar;
  outstruct.LCR(LCRcount).npar = npar;
  [nrows,nparnames]=fscanf(fid,"%g %g","C");
  outstruct.LCR(LCRcount).nrows = nrows;
  outstruct.LCR(LCRcount).nparnames = nparnames;
  outstruct.LCR(LCRcount).parnames = {};
  for ii=1:nparnames
    outstruct.LCR(LCRcount).parnames{ii}=fscanf(fid,"%s","C");
  endfor

endfunction
