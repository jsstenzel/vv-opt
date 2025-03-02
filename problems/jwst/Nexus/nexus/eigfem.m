function [phi,frq]=eigfem(K,M,varargin)
%eigfem		Symmetric solution of finite element eigenproblem
%
%	[phi,frq]=eigfem(K,M,nmodes)
%
%where,		K	is the stiffness matrix.
%		M	is the system mass matrix.
%               nmodes  the number of modes to be solved for (optional).
%		phi	are the mass normalized modeshapes.
%		frq	are the modal frequencies (rad/s).
%
%	The eigenproblem is solved by making the eigenproblem
%symmetric.  This is done by the following transformation:
%
%		phibar=M^(1/2)*phi
%		Kbar=M^(-1/2)'*K*M^(-1/2);
%		
%	where,	M^(1/2)		is the cholesky decomposition
%		eig(Kbar)	is the symmetric eigenproblem
%
%	The symmetric eigenproblem solution is more numerically
%robust.  This is important for large systems (>250 dofs) with 
%rigid body modes.
%	Guyan reduction is used to remove any dofs without mass.
%

% 
% Copyright 1992.  National Aeronautics and Space Administration,
% all rights reserved.
% Copyright 1993-2000.  California Institute of Technology, National
% Aeronautics and Space Administration.  All rights reserved.  US Government
% Sponsorship acknowledged.
% 
% THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
% INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTIBILITY
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
% CALTECH BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR
% ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
% WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTUOUS
% ACTION, ARISING OUT OF OR IN CONNECTION WITH THE ACCESS, USE OF
% PERFORMANCE OF THIS SOFTWARE.
% 

%History
%   ~Jun92  eldred:   conceived the symmetric solution to the 
% 		      eigenproblem
%   ~Jun92  jmelody:  used in script file format
%   1Oct93  jmelody:  written into function
%   8Feb94  jmelody:  added full() statement to eig()
%   9Feb94  jmelody:  fooled around w/ order of fprintf statements
%  28Aug94  jmelody:  include automatic guyan reduction for dofs
%		      w/o mass
%  14Apr95  jmelody:  added number of modes to be solved for as optional input
%  19Apr95  marie:    fixed indexing error associated with number of modes
%  23Apr96  lneedels: fixed bug with the number of modes index when guyan 
%                     reduction is automatically induced (jmelody)
%  17Dec97  lneedels  changed max(size( to length for matlab5
%  12May98  lneedels: changed call to complement to the matlab5 setdiff 
%                     command.  Order of input arguments had to be swapped.
%  14May98 lneedels:  modified for varargin for matlab5

% fprintf('\n');
%check for dofs w/o mass
  oset=find(diag(M) == 0);
  if (length(oset) > 0)
% disp('  EIGFEM:  Reducing out dofs w/o mass');
    aset=setdiff(1:length(M),oset);
    [K,M,Toa]=guyan(K,M,oset,aset);
  end
% disp('  EIGFEM:  Solve eigenproblem')
  M=chol(M);      %calculate M^(1/2)
  M=inv(M);        %M^(-1/2)
% disp('  EIGFEM:  Cholesky decomposition performed')
  K=M'*K*M;
              %K should be symmetric!! (unlike M\K)
  K=tril(K)+tril(K,-1)'; %this ensures K symmetric
% disp('  EIGFEM:  Kbar calculated')
  [phi,frq]=eig(full(K));
% disp('  EIGFEM:  phibar and eig found')
  [frq,in]=sort(sqrt(diag(frq)));   %sort the reduced model
  phi=phi(:,in);
  if (nargin == 2)      %give back all of the modes
    phi=M*phi;
  else                  %only the first nmodes
    nmodes=[varargin{1}];
    [row,col]=size(M);
    if (nmodes < col)
      phi=M*phi(:,1:nmodes);
    else
      phi=M*phi;
    end
  end
  if (length(oset) > 0)
% disp('  EIGFEM:  Expanding Modes');
    phi=unguyan(Toa,phi,aset,oset);
  end
% fprintf('\n');
