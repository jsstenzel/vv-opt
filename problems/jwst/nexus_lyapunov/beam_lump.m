function [ko,mo]=beam_lump(ni,xyz,prop,varargin)
%BEAM	computes and assembles stiffness and mass for beam 
%	elements with pin flags, axial, torsional, shear, and, 
%      2 transverse stiffnesses, using a lumped mass matrix 
%      formulation (mass on trans dofs only).
%
%		[Ko,Mo]=beam_lump(ni,xyz,prop,mat,bc,K,M,ti,tf) 
%
%  where,  
%          ni    is the (:,5), (:,7), or (:,10) nodal connectivity array, 
%                each row describes an element:
%                ni(:,1) is the element number
%                ni(:,2) is the node number of one end of the element
%                ni(:,3) is the node number of the other end of the element
%                ni(:,4) is the node number of a third node defining the
%			 plane for in-plane and out-of-plane bending
%                ni(:,5) is the property number, or row index into prop
%                ni(:,6) is the optional pin flag for the first end, ni(:,2).
%                ni(:,7) is the optional pin flag for the second end, ni(:,3).
%		 ni(:,8) is the optional first component of the orient vector
%                ni(:,9) is the optional second component of the orient vector 
%                ni(:,10) is the optional third component of the orient vector 
%                The orientation vector is used if the third node is assigned
%                a value of zero.  The vector is assumed to be in the output 
%                coordinate system of the first nodal point
%          xyz   is the array of nodal coordinates
%          prop  is a (:,10) array of element properties
%                prop(:,1) is property id
%                prop(:,2) is material id
%                prop(:,3) is cross-section area
%                prop(:,4) is in-plane bending moment of inertia, I22
%                prop(:,5) is out-of-plane bending moment of inertia, I33
%                prop(:,6) is torsion moment of inertia, J
%                prop(:,7) is the nonstructural mass on the element (OLD IMOS)
%                prop(:,8) is the nonstructural mass per length (NASTRAN-style)
%                prop(:,9) is the factor for multiplying the cross-section
%                          area to obtain the shear area in the beam local 2
%                          direction (toward the reference node).
%                prop(:,10) is the factor for multiplying the cross-section
%                          area to obtain the shear area in the beam local 3
%                          direction (cross-product of nodes 1-2 and 1-3). 
%       mat     the materials matrix, beam elements must use mat1 entries
%       bc      is the array of degree of freedom numbers
%       K,M     are the input system stiffness and mass matrices
%       ti      #-of-nodes by 1 and contains the index for local
%               coordinate transformations.  (optional)
%       tf      #-of-transforms by 9 and holds the 3x3 transformation
%               matrices stored by rows.  (optional)
%       Ko,Mo   are the returned system mass and stiffness matrices
%
%	The third node of ni defines the plane for in-plane and out-of-plane 
%bending.  In-plane bending refers to bending where the displacemnts are 
%in the plane defined by the third node.
%	The pin flag removes connections between the grid points and 
%selected degrees of the freedom of the bar (using the local element 
%coordinate system).  Up to five of the unique integers 1 through 6 
%with no imbedded blanks may be entered for each pin flag.  Note that 
%if pin flags are used on any of the beam elements, they must be 
%present (perhaps as 0) for all elements.

 
% Copyright 1992.  National Aeronautics and Space Administration,
% all rights reserved.
% Copyright 1993-1998.  California Institute of Technology, National
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
 
 
%HISTORY
%  14Apr92 hcb:       for sign error in local-to-global coord xform
%  17Jul92 BNorton:   use NASTRAN style consistent mass matrix
%                     and correct an error in the sign of the 2,6
%                     and 3,5 elements.
%  13Aug92 BNorton:   add shear flexibility
%  26Aug92 jmelody:   add nodal numbering
%  17Jun94 jmelody:   removed element type from calls to nifix.m
%  24Jun96 BNorton:   added pin flag capability (allows release of selected
%                     degrees of freedom at either end of the beam element)
%  16Jan96 lneedels:  major remoding, to make it suitable for matlab
%                     compiler
%  14Feb97 jmelody:   added cols input to nifix.m in order to fix bug 
%                     with simultaneous use of pin flags and arbitrary node 
%                     numbers.  Bug found by prapacz.
%  21Feb97 lneedels:  fixed bug with multiple loopings over i
%  26Mar97 lneedels:  removed local transformations
%  03Mar98 lneedels:  modified to use NASTRAN-like properties
%  14Mar98 lneedels:  modified to include local cs & orient vector ala akissil
%  27Apr98 lneedels:  added ability to handle zero cross-section area properties
%  14May98 lneedels:  modified for varargin for matlab5
%  15May98 lneedels:  changed "disp" to "warning" where appropriate
%  15Jun98 lneedels:  added nsm/length for NASTRAN compatibility
%  24Sep98 akissil:   modified the beam function to create a lumped mass version
%  12Apr99 akissil:   cleaned up for formal inclusion

%  old style, 6 arguments
%beam(ni,xyz,prop,bc,k,m)
%  old style, 6 arguments + local coordinates (8 arguments)
%beam(ni,xyz,prop,bc,k,m,ti,tf)
%  v3.1 style (7 arguments)
%beam(ni,xyz,prop,mat,bc,k,m)
%  v3.1 style + local coordinates (9 arguments)
%beam(ni,xyz,prop,mat,bc,k,m,ti,tf)

%  map arguments in IMOS v3.1 into varargin
mat=[varargin{1}];
bc=[varargin{2}];
k=[varargin{3}];
%  match arguments to Matlab5+ varargin
if nargin>6,
  m=[varargin{4}];
end
if nargin>7,
  ti=[varargin{5}];
end
if nargin>8,
  tf=[varargin{6}];
end

if nargin<8,
  ti=[];
  tf=[];
end

% backward compatibility pre-v3.1 prop/mat combination
if nargin==6 | nargin==8,  % old style or old style with local coord
  warning('beam elements changed under IMOS v3.1 - "help beam" for more info')
%  now mix & match to get correct matrices in new places
%  matrices must be moved over one position 
  if nargin==8,
    tf=ti;
    ti=m;
  end
  m=k;
  k=bc;
  bc=mat;
%  must check to make sure prop is "old prop", not "old old prop" (or earlier)
  [temp,propsz]=size(prop);
  if propsz==9, prop(:,10)=zeros(temp,1); end
  if propsz==8 prop(:,9:10)=zeros(temp,2); end
  if propsz==7 prop(:,8:10)=zeros(temp,3); end 

% must create mat matrix & shuffle prop matrix
% make a mat entry for every prop entry 
  [mprop,nprop]=size(prop);
%  check for an area of zero
  tmp=find(prop(:,1)==0);
  lntmp=length(tmp);
  if lntmp==0, % no properties with zeros cross-sectional area
%            mid           1           E        G           nu                rho!!             alpha
    mat=[ [1:mprop]' ones(mprop,1) prop(:,5) prop(:,6) zeros(mprop,1)  prop(:,7)./prop(:,1)  zeros(mprop,1)];
  else,  % have zero area element, make rho equal to zero for this
    warning('beam property with zero cross-sectional area has been detected')
    warning('rho (mass/volume) for this property is set equal to zero')
    mat=[ [tmp] ones(lntmp,1) prop(tmp,5) prop(tmp,6) zeros(lntmp,1)     zeros(lntmp,1)      zeros(lntmp,1)];
%  now do "normal" ones
    tmp=find(prop(:,1));
    lntmp=length(tmp);
    mat=[mat;[tmp ones(lntmp,1) prop(tmp,5) prop(tmp,6) zeros(lntmp,1)  prop(tmp,7)./prop(tmp,1)  zeros(lntmp,1)]];
%  now sort so they are ordered by number
    [y,idx]=sort(mat(:,1));
    mat=mat(idx,:);
  end

%  now rearrange prop
%            pid        mid        A         I1        I2        J     nsm/element   nsm/length      K1       K2
  prop=[ [1:mprop]'  [1:mprop]' prop(:,1) prop(:,3) prop(:,4) prop(:,2) prop(:,8)  zeros(mprop,1) prop(:,9) prop(:,10) ];

end

%  now check to make sure nsm/length exists
[mprop,nprop]=size(prop);
if nprop==9,
  prop=[prop(:,1:7)  zeros(mprop,1)  prop(:,8:9)];
end

%  arbitrary node numbering
[temp1,temp2]=size(xyz);
if temp2==4,
  index=nodesort(xyz);
  xyz=xyz(:,2:4);
  ni=nifix(index,ni,[2:4]);
end

%  get size of "active" ni
activ=find(ni(:,1)>0);
ni=ni(activ,:);  %  "throw away inactive ones"

%  get size of k
[mdof,ndof]=size(k);

[ivec,jvec,kvec,mvec]=beammex_lump(ni,xyz,prop,mat,bc,ti,tf);

ko=k+sparse(ivec,jvec,kvec,mdof,ndof);
mo=m+sparse(ivec,jvec,mvec,mdof,ndof);
