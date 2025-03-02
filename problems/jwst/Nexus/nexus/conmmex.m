function [ivec,jvec,mvec]=conmmex(em,nodes,bc)
 
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

%History
%   5Oct93 jmelody:  expanded error checking on size of nodes (make sure
%		     nodes is a matrix)
%   5Oct93 jmelody:  made more robust for vector and matrrix em input
%  10Oct93 jmelody:  arbitrary node numbering capabilities added
%  29Jan97 lneedels: modified for compiler
%   4Feb97 lneedels: modified so that warning messages are printed once
%                    per call, instead of once per node be removing indexing
%		     over nodes
%  17Dec97 lneedels: changed max(size( to length for matlab5
%  12May98 lneedels:  modified find's with row vector input to be Matlab5.1+
%                     compliant
%

[nem,mem]=size(em);
[nnodes,mnodes] = size(nodes);  % conm makes sure its a row vector

% preallocate maximum possible length of vectors
if (nem~=1 & mem~=1),  % em's a matrix
  ivec=zeros(mnodes*36,1);
else,  % em's a scalar or vector
  ivec=zeros(mnodes*6,1);
end

jvec=ivec;
kvec=ivec;
mvec=ivec;
cntr=0;

%
err=[sprintf('\n')];	%this string vector is for handling warnings

if nem == 1	%assume that em is translational mass
  [idx,jdx,dof] = find(bc(nodes,1:3));
  sz=length(dof);
  if sz>0,
    if mnodes==1,  % only one node sent in
      ivec(cntr+1:cntr+sz,1)=dof';
      jvec(cntr+1:cntr+sz,1)=dof';
    else
      ivec(cntr+1:cntr+sz,1)=dof;
      jvec(cntr+1:cntr+sz,1)=dof;
    end
    mvec(cntr+1:cntr+sz,1)=em*ones(sz,1);
    cntr=cntr+sz;
  end
  if sz<3*length(nodes),
    err=[err sprintf('Note:    Translational mass has not been added to inactive\n')];
    err=[err sprintf('         translational dofs.\n\n')];
  end

elseif (nem ~= 1)&(mem == 1)	%assume the vector, em, is the diagonal
				%of the full em matrix
  [idx,jdx,dof]=find(bc(nodes,:));  %active dofs for the node (1 thru 6)
  sz=length(dof);           %number of node dofs that are active

  if sz < nem*length(nodes),   %em is larger than active dofs
    if nem  ~= 6,   %if em isn't a six vector, I don't know where to put it
      err=sprintf('Since the vector em is not a six-vector and its size does not');
      err=[err sprintf(' correspond\nto the number of dofs for the nodes, ')];
      err=[err sprintf('I dont know where to put it!')];
      error(err);
    else,
      err=[err sprintf('Warning: you have tried to add mass to dofs that dont exist\n')];
      err=[err sprintf('         This mass has been ignored.\n\n')];
    end
  elseif sz > nem*length(nodes)  %em is smaller than active dofs
    error('The vector em is smaller than the number of active dofs of the nodes!');
  end

  if mnodes==1,  % only one node sent in
    ivec(cntr+1:cntr+sz,1)=dof';
    jvec(cntr+1:cntr+sz,1)=dof';
  else
    ivec(cntr+1:cntr+sz,1)=dof;
    jvec(cntr+1:cntr+sz,1)=dof;
  end

  mvec(cntr+1:cntr+sz,1)=em(jdx);
  cntr=cntr+sz;

elseif (nem~=1)&(mem~=1)          %em is assumed full
  for i=nodes                     %loop over nodes

    [idx,jdx,dof]=find(bc(i,:));  %active dofs for the node (1 thru 6)
    sz=length(dof);               %number of node dofs that are active

    if sz > nem,    %em is smaller than active dofs
      error('The matrix em is smaller than the number of active dofs of the nodes!');
    elseif sz < nem     %active dofs is smaller than em
      if nem ~= 6,      %if em isn't a six vector, I don't know where to put
        err=sprintf('Since the matrix em is not 6-by-6 and its size does not');
        err=[err sprintf(' correspond\nto the number of dofs for the nodes, ')];
        err=[err sprintf('I dont know where to put it!')];
        error(err);
      else,
        if length(err)==1,
          err=[err sprintf('Warning: you have tried to add mass to dofs that dont exist\n')];
          err=[err sprintf('         This mass has been ignored.\n\n')];
        end
      end
    end

    for idx=1:sz,
      ivec(cntr+(idx-1)*sz+1:cntr+(idx-1)*sz+sz,1)=dof';
      jvec(cntr+(idx-1)*sz+1:cntr+(idx-1)*sz+sz,1)=dof(idx)*ones(sz,1);
    end
    mvec(cntr+1:cntr+sz*sz,1)=reshape(em(jdx,jdx),sz*sz,1);
    cntr=cntr+sz*sz;
  end    
end
 
if length(err) > 1	%there is a warning message
  fprintf(err)
end

%  trim output vectors
ivec=ivec(1:cntr);
jvec=jvec(1:cntr);
kvec=kvec(1:cntr);
mvec=mvec(1:cntr);
