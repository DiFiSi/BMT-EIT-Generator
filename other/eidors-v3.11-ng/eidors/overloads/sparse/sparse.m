function S = sparse(varargin)
%SPARSE Create sparse matrix (EIDORS overload).
%   S = SPARSE(X) converts a sparse or full matrix to sparse form by
%   squeezing out any zero elements.
%
%  See also SPARFUN/SPARSE
%
% Sparse doesn't work for uint* inputs, so we need to preconvert to double

% (C) 2011 Bartlomiej Grychtol. License: GPL v2 or v3.
% $Id: sparse.m 6413 2022-11-26 20:29:19Z bgrychtol $

for i= 1:nargin
  if ~isa(varargin{i},'double')
     varargin{i} = double(varargin{i});
  end
end
S = builtin('sparse',varargin{:});
