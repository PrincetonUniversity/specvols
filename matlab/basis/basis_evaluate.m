% BASIS_EVALUATE Evaluate coefficient vector in basis
%
% Usage
%    x = basis_evaluate(basis, v);
%
% Input
%    basis: A basis object (for example, obtained from `dirac_basis`) with
%       respect to which the coefficients are defined.
%    v: A coefficient vector (or an array of coefficient vectors) to be
%       evaluated. The first dimension must equal `basis.count`.
%
% Output
%    x: The evaluation of the coefficient vector(s) `v` in the basis `basis`.
%       This is an array whose first dimensions equal `basis.sz` and the
%       remaining dimensions correspond to dimensions two and higher of `v`.
%
% Description
%    If `v` is a matrix of size `basis.ct`-by-..., `B` is the change-of-basis
%    matrix of `basis`, and `x` is a matrix of size `basis.sz`-by-..., the
%    function calculates
%
%       x = B * v
%
%    where the rows of `B` and columns of `x` are read as vectorized arrays.
%
% See also
%    basis_expand

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = basis_evaluate(basis, v)
    if size(v, 1) ~= basis.count
        error('First dimension of `v` must be of size `basis.count`.');
    end

    if basis.type == basis_type_dirac()
        x = dirac_basis_evaluate(basis, v);
    elseif basis.type == basis_type_fb()
        x = fb_basis_evaluate(basis, v);
    elseif basis.type == basis_type_matrix()
        x = matrix_basis_evaluate(basis, v);
    else
        error('Invalid basis type.');
    end
end
