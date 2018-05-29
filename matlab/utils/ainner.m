% AINNER Array inner product along given dimensions
%
% Usage
%    c = ainner(A, B, calc_dims);
%
% Input
%    A, B: Two arrays of arbitrary size and shape. The shapes of both arrays
%       have to coincide.
%    calc_dims: A subset of `1:ndims(A)` along which to compute the inner
%       product. If empty (the default), the product is calculated along all
%       dimensions.
%
% Output
%    c: The inner product of `A` and `B` along the dimensions specified by
%       `calc_dims`.
%
% See also
%    asum, anorm, acorr

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function c = ainner(A, B, calc_dims)
    if nargin < 3
        calc_dims = [];
    end

    if ndims(A) ~= ndims(B) || any(size(A) ~= size(B))
        error('Input arrays must have the same shape');
    end

    c = asum(A.*conj(B), calc_dims);
end
