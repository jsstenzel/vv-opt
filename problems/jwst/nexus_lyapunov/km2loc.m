function [k,m]=km2loc(k,m,ti,tf)
%KM2LOC   [k,m]=km2loc(k,m,ti,tf)
%   Transforms the stiffness and mass matrices from basic to local
%   output coordinate systems.
%   k,m= ndof x ndof stiffness and mass matrices
%        where ndof= #nodes x 6
%   ti= #nodes x 1 index referencing transform matrix location
%   tf= #coord sys x 9 matrix storing 3x3 tranf matrices by row
%

% Version 2.0 1997, Copyright California Institute of Technology, National
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

%   HISTORY   APR  1997  akissil   created

%disp('Transforming mass and stiffness to local coordinates')
T=gentloc(ti,tf);
k=T*k*T';
m=T*m*T';
