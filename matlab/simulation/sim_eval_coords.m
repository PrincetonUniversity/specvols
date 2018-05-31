% SIM_EVAL_COORDS Evaluate coordinate estimation
%
% Usage
%    coords_perf = sim_eval_coords(sim, mean_vol, eig_vols, coords_est);
%
% Input
%    sim: Simulation object from `create_sim`.
%    mean_vol: A mean volume in the form of an L-by-L-by-L array.
%    eig_vols: A set of eigenvolumes in an L-by-L-by-L-by-K array.
%    coords_est: The estimated coordinates in the affine space defined
%       centered at `mean_vol` and spanned by `eig_vols`.
%
% Output
%    coords_perf: A struct containing the evaluation results. It contains the
%       fields:
%          - rel_err: The relative error of the volumes obtained by expanding
%             the coordinates in the affine space.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function coords_perf = sim_eval_coords(sim, mean_vol, eig_vols, coords_est)
    [coords_true, residuals] = sim_vol_coords(sim, mean_vol, eig_vols);

    coords_true = coords_true(:,sim.states);
    residuals = residuals(sim.states);

    err_coords = anorm(coords_true - coords_est, 1);
    err = hypot(residuals, err_coords);

    norm_true = hypot(residuals, anorm(coords_true, 1));

    rel_err = err./norm_true;

    coords_perf = struct();

    coords_perf.rel_err = rel_err(:);
end
