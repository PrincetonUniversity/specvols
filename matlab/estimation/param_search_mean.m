% PARAM_SEARCH_MEAN Cross-validate regularizers for mean estimation
%
% Usage
%    residuals = param_search_mean(n_partition, mean_kernels_f, ...
%       mean_b_coeffs, basis, precond_kernels_f, regularizers, mean_est_opt);
%
% Input
%    n_partition: A vector of length n_folds containing the number of images
%       in each fold of the partition. If `partition_idx` is a partition
%       obtained from `create_partition`, this vector can be obtained using
%       `cellfun(@numel, partition_idx)`.
%    mean_kernels_f: An array of size 2*L-by-2*L-by-2*L-by-n_folds containing
%       the mean kernels for the different folds. Typically obtained from
%       `src_partition_mean_kernel`. For more information, see
%       `conj_grad_mean`.
%    mean_b_coeffs: An array of size basis.count-by-n_folds containing
%       backprojected coefficients for the different folds. Typically obtained
%       from `src_partition_mean_backward`. For more information, see
%       `conj_grad_mean`.
%    basis: A basis object used for representing the volumes.
%    precond_kernels_f: The Fourier transforms the preconditioning kernels in
%       the form of a L-by-L-by-L-by-n_folds array. For more information, see
%       `conj_grad_mean`.
%    regularizers: A cell array of regularizers to evaluate. These are of the
%       same type as the `regularizer` field of `mean_est_opt` when supplied
%       to `conj_grad_mean`. If all regularizers are scalars (implying
%       isotropic Tikhonov regularization), the cell array may be in the form
%       of a regular array. For more information on possible regularizers, see
%       `conj_grad_mean`.
%    mean_est_opt: A struct of mean estimation options. For more details, see
%       `conj_grad_mean`.
%
% Output
%    residuals: An array of size n_folds-by-numel(regularizers) containing the
%       residuals for the different folds and choice of regularizer. The
%       residual is defined as
%
%          | im - im_est |^2 - | im |^2,
%
%       where im is the set of held-out images and im_est is the set of
%       projections of the mean estimated from the "held-in" images.
%
% See also
%    src_partition_mean_kernel, src_partition_mean_backward, create_partition

function residuals = param_search_mean(n_partition, mean_kernels_f, ...
    mean_b_coeffs, basis, precond_kernels_f, regularizers, mean_est_opt)

    if nargin < 7 || isempty(mean_est_opt)
        mean_est_opt = struct();
    end

    % We may want to take in an array of regularizers, in which case the
    % standard isotropic Tikhonov regularization is used.
    if isnumeric(regularizers)
        regularizers = num2cell(regularizers);
    end

    n_folds = size(mean_kernels_f, 4);

    % Un-normalize the kernels and coefficients. This will simplify taking
    % sums later.
    mean_kernels_f = ...
        bsxfun(@times, mean_kernels_f, permute(n_partition, [2 3 4 1]));
    mean_b_coeffs = ...
        bsxfun(@times, mean_b_coeffs, permute(n_partition, [2 1]));
    precond_kernels_f = ...
        bsxfun(@times, precond_kernels_f, permute(n_partition, [2 3 4 1]));

    for k = 1:n_folds
        % Compute kernel and backprojected coefficients for the held-out
        % portion of the data.
        kernel_out_f = mean_kernels_f(:,:,:,k);
        b_coeff_out = mean_b_coeffs(:,k);

        % Compute kernels and coefficients for the the remainder (i.e., the
        % "held-in" portion). Note that the kernels and coefficients have to
        % be re-normalized, otherwise their sizes will differ with respect to
        % the regularization term.
        idx = ([1:n_folds] ~= k);
        kernel_in_f = ...
            1/sum(n_partition(idx))*sum(mean_kernels_f(:,:,:,idx), 4);
        b_coeff_in = ...
            1/sum(n_partition(idx))*sum(mean_b_coeffs(:,idx), 2);
        precond_kernel_in_f = ...
            1/sum(n_partition(idx))*sum(precond_kernels_f(:,:,:,idx), 4);

        % For each of the regularizers, calculate the estimated mean on the
        % held-in dataset, then compute the residual with respect to the
        % held-out data.
        for ell = 1:numel(regularizers)
            mean_est_opt_ell = mean_est_opt;
            mean_est_opt.regularizer = regularizers{ell};

            est_coeff_in = conj_grad_mean(kernel_in_f, ...
                b_coeff_in, basis, precond_kernel_in_f, mean_est_opt_ell);

            % Calculate the residual as |im - im_est|^2 - |im|^2, where im is
            % the set of experimental images in the held-out set and im_est
            % consists of the corresponding projections of the mean. This
            % quantity can be efficiently calculated using the mean kernel of
            % the held-out data and the backprojected coefficients.
            est_coeff_in_conv = apply_mean_kernel(est_coeff_in, ...
                kernel_out_f, basis);
            residuals(k,ell) = est_coeff_in'*est_coeff_in_conv + ...
                -2*est_coeff_in'*b_coeff_out;
        end
    end
end
