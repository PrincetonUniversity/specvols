% SRC_COVAR_BACKWARD Apply adjoint covariance mapping to source
%
% Usage
%    covar_b_coeff = src_covar_backward(src, basis, mean_vol, noise_var, ...
%       covar_est_opt);
%
% Input
%    src: A source structure containing the images and imaging parameters.
%       This is typically obtained from `star_to_src` or `sim_to_src`.
%    basis: A basis object used for representing the volumes.
%    mean_vol: The (estimated) mean volume of the source as an L-by-L-by-L
%       array. This can be estimated using `estimate_mean`.
%    noise_var: The variance of the noise.
%    covar_est_opt: A struct containing the fields:
%          - 'precision': The precision of the kernel. Either 'double' or
%             'single' (default).
%          - 'batch_size': The size of the batches in which to compute the
%             kernel (default 512).
%          - 'shrinker': Type of shrinker to use. Should be one of 'none'
%             (default), 'frobenius_norm', 'operator_norm', or
%             'soft_threshold'. If set to 'none', it subtracts the expected
%             noise contribution from `covar_b_coeff`. Otherwise,
%             `shrink_covar` is applied using the given shrinker.
%
% Output
%    covar_b_coeff: The sum of the outer products of the mean-subtracted
%       images in `src`, corrected by the expected noise contribution and
%       expressed as coefficients of `basis`.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function covar_b_coeff = src_covar_backward(src, basis, mean_vol, ...
    noise_var, covar_est_opt)

    if nargin < 5 || isempty(covar_est_opt)
        covar_est_opt = struct();
    end

    covar_est_opt = fill_struct(covar_est_opt, ...
        'precision', 'single', ...
        'batch_size', 512, ...
        'shrinker', 'none');

    covar_b = zeros(src.L*ones(1, 6), covar_est_opt.precision);

    for batch = 1:ceil(src.n/covar_est_opt.batch_size)
        batch_s = (batch-1)*covar_est_opt.batch_size+1;
        batch_n = min(batch*covar_est_opt.batch_size, src.n)-batch_s+1;

        im = src_image(src, batch_s, batch_n);

        im_centered = im - vol_forward(src, mean_vol, batch_s, batch_n);

        im_centered_b = zeros([src.L*ones(1, 3) batch_n], covar_est_opt.precision);

        for s = 1:batch_n
            im_centered_b(:,:,:,s) = ...
                im_backward(src, im_centered(:,:,s), batch_s+s-1);
        end

        im_centered_b = vol_to_vec(im_centered_b);

        covar_b = covar_b + 1/src.n*vecmat_to_volmat(im_centered_b*im_centered_b');
    end

    covar_b_coeff = basis_mat_evaluate_t(basis, covar_b);

    mean_kernel_f = src_mean_kernel(src, covar_est_opt);
    An = basis_mat_evaluate_t(basis, kernel_to_toeplitz(mean_kernel_f));

    if strcmp(covar_est_opt.shrinker, 'none')
        covar_b_coeff = covar_b_coeff - noise_var*An;
    else
        An_sqrt = sqrtm(An);

        gamma = size(covar_b_coeff, 1)/src.n;

        covar_b_coeff = An_sqrt\covar_b_coeff/An_sqrt;
        covar_b_coeff = shrink_covar(covar_b_coeff, noise_var, gamma, ...
            covar_est_opt.shrinker);
        covar_b_coeff = An_sqrt*covar_b_coeff*An_sqrt;
    end
end
