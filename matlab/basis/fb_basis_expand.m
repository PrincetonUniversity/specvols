% FB_BASIS_EXPAND Expand array in Fourier-Bessel basis
%
% Usage
%    v = fb_basis_expand(basis, x);
%
% Input/Output
%    See documentation for `basis_expand`.
%
% See also
%    basis_expand

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function v = fb_basis_expand(basis, x)
    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    b = fb_basis_evaluate_t(basis, x);

    A = @(v)(fb_basis_evaluate_t(basis, fb_basis_evaluate(basis, v)));

    % TODO: Check that this tolerance make sense for multiple columns in x

    cg_opt.max_iter = Inf;
    cg_opt.rel_tolerance = 10*eps(class(x));
    cg_opt.verbose = 0;

    v = conj_grad(A, b, cg_opt);

    v = roll_dim(v, sz_roll);
end
