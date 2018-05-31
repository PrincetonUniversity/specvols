% SIM_EVAL_EIGS Evaluate covariance eigendecomposition accuracy
%
% Usage
%    eigs_perf = sim_eval_eigs(sim, eigs_est, lambdas_est);
%
% Input
%    sim: Simulation object from `create_sim`.
%    eigs_est: The estimated volume eigenvectors in an L-by-L-by-L-by-K array.
%    lambdas_est: The estimated eigenvalues in a K-by-K diagonal matrix
%       (default `diag(ones(K, 1))`).
%
% Output
%    eigs_perf: A struct containing the evaluation results. It contains the
%       fields:
%          - rel_err: The relative error of estimated volume matrix formed by
%             multiplying out the eigendecomposition.
%          - corr: The correlations of the estimated volume matrix.
%          - cos_principal_angles: The cosines of the principal angles between
%             the subspaces spanned by `eigs_est` and the true eigenvectors of
%             the simulation's covariance matrix.
%
% See also
%    sim_eval_covar, eval_eigs

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function eigs_perf = sim_eval_eigs(sim, eigs_est, lambdas_est)
    if nargin < 3
        lambdas_est = [];
    end

    [eigs_true, lambdas_true] = sim_eigs(sim);

    eigs_perf = eval_eigs(eigs_true, lambdas_true, eigs_est, lambdas_est);
end
