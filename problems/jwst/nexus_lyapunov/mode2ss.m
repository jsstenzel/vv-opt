function    [a,b,c,d]=mode2ss(omega,phi,zeta,lm,lb,lc)
% MODE2SS    [a,b,c,d]=mode2ss(omega,phi,zeta,lm,lb,lc)
%     forms a state space model from a finite element
%     modal model.  No order reduction except modal
%     selection is done.  
%
%         Modes are presumed normalized to unity modal mass.  If this is
%         not already done, then do this:
%
%               phi = chol(m) \ phi
%
%         where m is the mass matrix.
%  omega is vector of modal frequencies.
%  phi is array of mode shapes in columns.
%  zeta is a vector of modal damping values.
%      Note zeta and omega must be the same shape, i.e. both row vectors
%      or both column vectors.
%  lm is vector of mode numbers to include in ss model.
%  lb is vector of dofs for force input.
%  lc is vector of dofs for velocity sensor output.
%
%       For position sensors swap the left and right halves of c as follows:
%
%          c(:,[nc+1:2*nc,1:nc]) = c;
%
%       where nc is the number of columns of C (also the size of the state).
%
%  State space model is 
%    xdot = a x + b u
%       y = c x + d u
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
%   6Oct93 jmelody:  omega and zeta don't have to both be rows or
%		     columns
%  18Feb86 lneedels: modified for compilation

[a,b,c,d]=mode2ssmex(omega,phi,zeta,lm,lb,lc);
