% SIM_EVAL_COVAR Evaluate covariance estimation accuracy
%
% Usage
%    covar_perf = sim_eval_covar(sim, covar_est);
%
% Input
%    sim: Simulation object from `create_sim`.
%    covar_est: The estimated covariance matrix in the form of an
%       L-by-L-by-L-by-L-by-L-by-L array.
%
% Output
%    covar_perf: A struct containing the evaluation results. It contains the
%       fields:
%          - rel_err: The relative error of the estimate with respect to the
%             true covariance.
%          - corr: The correlation of the estimate with the true covariance.
%
% See also
%    eval_volmat

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function covar_perf = sim_eval_covar(sim, covar_est)
    covar_true = sim_covar(sim);

    covar_perf = eval_volmat(covar_true, covar_est);
end
