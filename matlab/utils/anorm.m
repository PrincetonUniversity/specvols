% ANORM Array norm along given dimensions
%
% Usage
%    c = anorm(A, calc_dims);
%
% Input
%    A: An array of arbitrary size and shape.
%    calc_dims: A subset of `1:ndims(A)` along which to compute the norm. If
%       empty (the default), the norm is calculated along all dimensions.
%
% Output
%    c: The Euclidean (l^2) norm of `A` along the dimensions specified by
%       `calc_dims`.
%
% See also
%    asum, ainner, acorr

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function c = anorm(A, calc_dims)
    if nargin < 2
        calc_dims = [];
    end

    c = sqrt(ainner(A, A, calc_dims));
end
