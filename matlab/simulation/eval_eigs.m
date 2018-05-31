% EVAL_EIGS Evaluate accuracy of estimated volume eigendecomposition
%
% Usage
%    eigs_perf = eval_eigs(eigs_true, lambdas_true, eigs_est, lambdas_est);
%
% Input
%    eigs_true: The true volume eigenvectors in an L-by-L-by-L-by-K1 array.
%    lambdas_true: The true eigenvalues in a K1-by-K1 diagonal matrix (default
%       `diag(ones(K1, 1))`).
%    eigs_est: The estimated volume eigenvectors in an L-by-L-by-L-by-K2 array.
%    lambdas_est: The estimated eigenvalues in a K2-by-K2 diagonal matrix
%       (default `diag(ones(K2, 1))`).
%
% Output
%    eigs_perf: A struct containing the evaluation results. It contains the
%       fields:
%          - rel_err: The relative error of the volume matrices formed by
%             multiplying out the eigendecompositions.
%          - corr: The correlations of the same volume matrices.
%          - cos_principal_angles: The cosines of the principal angles between
%             the subspaces spanned by `eigs_true` and `eigs_est`.
%
% See also
%    eval_volmat

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function eigs_perf = eval_eigs(eigs_true, lambdas_true, eigs_est, lambdas_est)
    sz_true = size(eigs_true);
    sz_est = size(eigs_est);

    if ndims(eigs_true) ~= ndims(eigs_est) || ...
        any(sz_true(1:3) ~= sz_est(1:3))

        error('Volume eigenvectors must be the same size.');
    end

    if isempty(lambdas_true)
        lambdas_true = diag(ones(size(eigs_true, 4), 1));
    end

    if isempty(lambdas_est)
        lambdas_est = diag(ones(size(eigs_est, 4), 1));
    end

    eigs_perf = struct();

    B = vol_to_vec(eigs_est)'*vol_to_vec(eigs_true);

    norm_true = anorm(lambdas_true);
    norm_est = anorm(lambdas_est);

    inner = ainner(B*lambdas_true, lambdas_est*B);

    err = sqrt(norm_true^2 + norm_est^2 - 2*inner);

    rel_err = err/norm_true;
    corr = inner/(norm_true*norm_est);

    cos_principal_angles = cos_principal_angles(vol_to_vec(eigs_true), ...
        vol_to_vec(eigs_est));

    eigs_perf.rel_err = rel_err;
    eigs_perf.corr = corr;
    eigs_perf.cos_principal_angles = cos_principal_angles;
end
