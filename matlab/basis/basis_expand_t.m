% BASIS_EXPAND_T Expand array in dual basis
%
% Usage
%    x = basis_expand_t(basis, v);
%
% Input
%    basis: A basis object (for example, obtained from `dirac_basis`) in whose
%       dual basis the array is to be expanded.
%    v: An array whose first dimension is to be expanded in the dual basis.
%       This dimension must be equal to `basis.count`.
%
% Output
%    x: The coefficients of `v` expanded in the dual of `basis`. If more than
%       one vector is supplied in `v`, the higher dimensions of `x` correspond
%       to second and higher dimensions of `v`.
%
% Description
%    If `v` is a matrix of size `basis.ct`-by-..., `B` is the change-of-basis
%    matrix of `basis`, and `x` is a matrix of size `basis.sz`-by-..., the
%    function calculates
%
%       x = (B * B')^(-1) * B * v
%
%    where the rows of `B` and columns of `x` are read as vectorized arrays.
%
% See also
%    basis_expand

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = basis_expand_t(basis, v)
    if size(v, 1) ~= basis.count
        error('First dimension of `v` must be of size `basis.count`.');
    end

    if basis.type == basis_type_dirac()
        x = dirac_basis_evaluate(basis, v);
    elseif basis.type == basis_type_fb()
        x = fb_basis_expand_t(basis, v);
    elseif basis.type == basis_type_matrix()
        x = matrix_basis_expand_t(basis, v);
    elseif basis.type == basis_type_ffb()
        error('Not implemented.');
    else
        error('Invalid basis type.');
    end
end
