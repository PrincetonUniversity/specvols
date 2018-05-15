% BASIS_EXPAND Expand array in basis
%
% Usage
%    v = basis_expand(basis, x);
%
% Input
%    basis: A basis object (for example, obtained from `dirac_basis`) in which
%       the array is to be expanded.
%    x: An array whose first few dimensions are to be expanded in `basis`.
%       These dimensions must equal `basis.sz`.
%
% Output
%    v: The coefficients of `x` expanded in `basis`. If more than one array
%       of size `basis.sz` is found in `x`, the second and higher dimensions
%       of `v` correspond to those higher dimensions of `x`.
%
% Description
%    If `v` is a matrix of size `basis.ct`-by-..., `B` is the change-of-basis
%    matrix of `basis`, and `x` is a matrix of size `basis.sz`-by-..., the
%    function calculates
%
%       v = (B' * B)^(-1) * B' * x
%
%    where the rows of `B` and columns of `x` are read as vectorized arrays.
%
% See also
%    basis_evaluate

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function v = basis_expand(basis, x)
    sz_x = size(x);

    if numel(sz_x) < numel(basis.sz) || ...
        ~all(sz_x(1:numel(basis.sz)) == basis.sz)
        error('First dimensions of `x` must match `basis.sz`.');
    end

    if basis.type == basis_type_dirac()
        v = dirac_basis_expand(basis, x);
    elseif basis.type == basis_type_fb()
        v = fb_basis_expand(basis, x);
    else
        error('Invalid basis type.');
    end
end
