function [kr,mr,toa] = guyan(k,m,oset,aset)
% GUYAN    [kr,mr,toa]=guyan (k,m,oset,aset)    computes the guyan reduction
%       transformation matrix and reduces the system mass and stiffness
%       matrices.
%
%   k is the system stiffness matrix.
%
%   m is the system mass matrix.
%
%   oset is a vector of equation numbers that represent the omitted
%       dofs.
%
%   aset is a vector of equation numbers that represents the retained or
%       analysis set of dofs.
%
%   toa is the transformation between the retained and omitted dofs.
%
%       Note that usually the union of oset and aset accounts for all the
%       dofs but this is not enforced.
%
%       Note that if k(oset,oset) has rigid body modes in it then Matlab
%       doesn't complain but this should be avoided

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

%       Let the unknowns x be [xa; xo] so that k is partitioned
%           k= [kaa, kao; koa, koo] and the applied forces f
%           are partitioned f = [fa ; fo].
%       For a static problem f = k x is written
%           |fa| =  | kaa   kao |  |xa|
%           |fo|    | koa   koo |  |xo|
%
%       Solve for xo = koo^-1 ( fo - koa xa ) and
%           kaa xa = fa - kao xo = fa - kao koo^-1 (fo -koa xa) 
%           (kaa - kao koo^-1 koa) xa = fa - kao koo^-1 fo.
%       The reduced stiffness matrix is kaa - kao koo^-1 koa.
%       The transformation is xo = -koo^-1 koa xa  so
%           toa = -koo^-1 koa.
%       The reduced stiffness matrix can be found to be
%       | Iaa toa' | | kaa kao| | Iaa | = | Iaa toa'| | kaa+kao toa | 
%                    | koa koo| | toa |               | koa+koo toa |
%
%       = kaa+kao toa + toa'(koa+koo toa)
%       = kaa - kao koo^-1 koa -kao koo^-1 ( koa - koo koo^-1 koa)
%       = kaa - kao koo^-1 koa !
% 
%      The reduced mass matrix is
%         maa+mao toa + toa'(moa+moo toa)
%       = maa+mao toa + toa'moa + toa'moo toa
%
 
toa = -k(oset,oset)\k(oset,aset);
kr=k(aset,aset) + k(aset,oset)*toa;
mr=m(aset,aset) + m(aset,oset)*toa + toa'*m(oset,aset) + toa'*m(oset,oset)*toa;
