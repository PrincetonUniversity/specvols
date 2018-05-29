% ACORR Array correlation along given dimensions
%
% Usage
%    c = acorr(A, B, calc_dims);
%
% Input
%    A, B: Two arrays of arbitrary size and shape. The shapes of both arrays
%       have to coincide.
%    calc_dims: A subset of `1:ndims(A)` along which to compute the correla-
%       tion. If empty (the default), the correlation is calculated along all
%       dimensions.
%
% Output
%    c: The correlation of `A` and `B` along the dimensions specified by
%       `calc_dims`.
%
% See also
%    asum, ainner, anorm

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function c = acorr(A, B, calc_dims)
    if nargin < 3
        calc_dims = [];
    end

    c = ainner(A, B, calc_dims)./(anorm(A, calc_dims).*anorm(B, calc_dims));
end
