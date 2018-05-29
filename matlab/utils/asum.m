% ASUM Array sum along given dimensions
%
% Usage
%    c = asum(A, calc_dims);
%
% Input
%    A: An array of arbitrary size and shape.
%    calc_dims: A subset of `1:ndims(A)` along which to compute the sum. If
%       empty (the default), the sum is calculated along all dimensions.
%
% Output
%    c: The sum of `A` along the dimensions specified by `calc_dims`.
%
% See also
%    ainner, anorm, acorr

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function c = asum(A, calc_dims)
    if nargin < 2 || isempty(calc_dims)
        calc_dims = 1:ndims(A);
    end

    calc_dims = calc_dims(ismember(calc_dims, 1:ndims(A)));

    sz = size(A);
    dims = 1:ndims(A);
    idx_dims = dims(~ismember(dims, calc_dims));

    perm = [calc_dims idx_dims];

    A = permute(A, perm);

    A = reshape(A, prod(sz(calc_dims)), []);

    c = sum(A, 1);

    c = reshape(c, [ones(1, numel(calc_dims)) sz(idx_dims)]);
    c = permute(c, invperm(perm));
end
