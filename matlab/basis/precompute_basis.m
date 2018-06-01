% PRECOMPUTE_BASIS Precompute basis vectors and store in matrix
%
% Usage
%    basis = precompute_basis(original_basis);
%
% Input
%    original_basis: A basis object (for example, obtained from `dirac_basis`
%       or `fb_basis`) whose basis vectors we want to precompute.
%
% Output
%    basis: A matrix basis object corresponding to `original_basis`, but with
%       all basis vectors precomputed and stored in a matrix for fast
%       evaluation.
%
% Description
%    The returned basis object has the exact same properties as the supplied
%    `original_basis`, but is faster to evaluate and expand since all of its
%    basis vectors have already been computed. Note, however, that this comes
%    with an increase in required memory to store the matrix of basis vectors.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function basis = precompute_basis(original_basis)
    if original_basis.type == basis_type_dirac
        basis = original_basis;
        return;
    end

    B = basis_evaluate(original_basis, eye(original_basis.count));

    B = reshape(B, [prod(original_basis.sz) original_basis.count]);

    basis = struct();

    basis.type = basis_type_matrix;

    basis.sz = original_basis.sz;
    basis.count = original_basis.count;

    basis.B = B;
end
