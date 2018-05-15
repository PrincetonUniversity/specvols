% FB_BASIS_EXPAND_T Expand array in dual Fourier-Bessel basis
%
% Usage
%    x = fb_basis_expand_t(basis, v);
%
% Input/Output
%    See documentation for `basis_expand_t`.
%
% See also
%    basis_expand_t

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = fb_basis_expand_t(basis, v)
    d = numel(basis.sz);

    [v, sz_roll] = unroll_dim(v, 2);

    if d == 2
        b = im_to_vec(fb_basis_evaluate(basis, v));

        A = @(x)(im_to_vec(fb_basis_evaluate(basis, ...
            fb_basis_evaluate_t(basis, vec_to_im(x)))));
    elseif d == 3
        b = vol_to_vec(fb_basis_evaluate(basis, v));

        A = @(x)(vol_to_vec(fb_basis_evaluate(basis, ...
            fb_basis_evaluate_t(basis, vec_to_vol(x)))));
    end

    % TODO: Check that this tolerance make sense for multiple columns in x

    cg_opt.max_iter = Inf;
    cg_opt.rel_tolerance = 10*eps(class(v));
    cg_opt.verbose = 0;

    x = conj_grad(A, b, cg_opt);

    x = roll_dim(x, sz_roll);

    if d == 2
        x = vec_to_im(x);
    elseif d == 3
        x = vec_to_vol(x);
    end
end
