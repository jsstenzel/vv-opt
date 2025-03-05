function [xyzcg]=cg_calc(m,xyz,bc,varargin)
%CG_CALC	finds the center of mass of a finite element model.
%
%  [xyzcg]=cg_calc(m,xyz,bc,ti,tf)
%	
%  m is the system mass matrix.
%  xyz is the nodal coordinates array.
%  bc is the array of dofs.
%  ti is the  #nodes x 1 index vector referencing transform location (optional).
%  tf is the #coord syst x 9 matrix for storing 3x3 transf matrices (optional).
%  xyzcg is the location of the center of mass.
%
%  xyzcg is the average of the two values obtained for each of the 
%  coordinates.  If the difference between the two values obtained
%  for a coordinate cm location is greater than 1e-6, a warning is issued
%  and the all of calculated values are printed.  This difference 
%  occurs when "non-real" mass is included in a model.  As an example,
%  a constrained model using the (consistent mass) beam may trigger this 
%  warning.


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
%  14Aug92 jmelody:  copied from HCB
%  14Aug92 jmelody:  arbitrary node numbering capabilities added
%   8Apr94 jmelody:  changed name from cg.m to cg_calc.m
%  12Feb97 lneedels: modified for compilation
%  09Jun98 akissil:  substantial rewrite using wtcgmex function 
%  25Jun98 lneedels: modified for varargin
%  15Jul98 lneedels: added error checking to make sure calculate cm values
%                    are the same
%  06Aug98 lneedels: added flag which will print out the 3x3 cm matrix 
%                    if any of the cm values don't match
%

flag=0;

[mxyz,nxyz]=size(xyz);

if nxyz == 4		%non-sequencial nodal numbering
  xyz=xyz(:,2:4);	%throw out nodal number column
end

Rpt=[0,0,0];        %define reference point at the origin

if nargin==3,
  ti=[];
  tf=[];
elseif nargin==4,
  error('Must include both ti & tf if you want to use local coordinates')
elseif nargin==5,
  ti=[varargin{1}];
  tf=[varargin{2}];
end

dg=wtcgmex(bc,xyz,Rpt,ti,tf);

mpt=dg'*m*dg;

% Compute the offset of the cm from the origin
Xy=mpt(2,6)/mpt(2,2);
Xz=-mpt(3,5)/mpt(3,3);
Yx=-mpt(1,6)/mpt(1,1);
Yz=mpt(3,4)/mpt(3,3);
Zx=mpt(1,5)/mpt(1,1);
Zy=-mpt(2,4)/mpt(2,2);

% average calculated values
xyzcg(1,1)=0.5*(Xy+Xz);
xyzcg(1,2)=0.5*(Yx+Yz);
xyzcg(1,3)=0.5*(Zx+Zy);

%  check to see if any of the cm points are not the same
tol=1e-6;
if abs(xyzcg(1,1))>tol,
  val=abs(Xy-Xz)/xyzcg(1,1);
else,
  val=abs(Xy-Xz);
end
if val>tol
  str=['Difference of X cm values is ', num2str(val*xyzcg(1,1))];
  warning(str);
  flag=1;
end

if abs(xyzcg(1,2))>tol,
  val=abs(Yx-Yz)/xyzcg(1,2);
else,
  val=abs(Yx-Yz);
end
if val>tol
  str=['Difference of Y cm values is ', num2str(val*xyzcg(1,2))];
  warning(str);
  flag=1;
end

if abs(xyzcg(1,3))>tol,
  val=abs(Zx-Zy)/xyzcg(1,3);
else,
  val=abs(Zx-Zy);
end
if val>tol
  str=['Difference of Z cm values is ', num2str(val*xyzcg(1,3))];
  warning(str);
  flag=1;
end

%  print whole matrix if cm values don't match
if flag==1,
  [ 0 Yx Zx
   Xy  0 Zy
   Xz Yz 0]
end
