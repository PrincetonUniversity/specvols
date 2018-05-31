% SIM_EVAL_MEAN Evaluate mean estimation accuracy
%
% Usage
%    mean_perf = sim_eval_mean(sim, mean_est);
%
% Input
%    sim: Simulation object from `create_sim`.
%    mean_est: The estimated mean volume in the form of an L-by-L-by-L array.
%
% Output
%    mean_perf: A struct containing the evaluation results. It contains the
%       fields:
%          - rel_err: The relative error of the estimate with respect to the
%             true mean.
%          - corr: The correlation of the estimate with the true mean.
%
% See also
%    eval_vol

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function mean_perf = sim_eval_mean(sim, mean_est)
    mean_true = sim_mean(sim);

    mean_perf = eval_vol(mean_true, mean_est);
end
