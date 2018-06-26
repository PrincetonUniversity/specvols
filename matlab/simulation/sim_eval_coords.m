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
%          - corr: The corrleations of the volumes obtained by expanding the
%             coordinates in the affine space.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function coords_perf = sim_eval_coords(sim, mean_vol, eig_vols, coords_est)
    [coords_true, res_norms, res_inners] = sim_vol_coords(sim, mean_vol, ...
        eig_vols);

    coords_true = coords_true(:,sim.states);
    res_norms = res_norms(sim.states);
    res_inners = res_inners(sim.states);

    mean_eigs_inners = vol_to_vec(mean_vol)'*vol_to_vec(eig_vols);

    coords_err = coords_true - coords_est;

    err = anorm(coords_err, 1);
    err = hypot(res_norms, err);

    norm_true = sqrt(anorm(coords_true, 1).^2 + anorm(mean_vol)^2 + ...
        2*res_inners + 2*mean_eigs_inners*coords_true);
    norm_true = hypot(res_norms, norm_true);

    rel_err = err./norm_true;

    inner = anorm(mean_vol)^2 + ...
        mean_eigs_inners*(coords_true+coords_est) + ...
        sum(coords_true .* coords_est, 1) + res_inners;

    norm_est = sqrt(anorm(coords_est, 1).^2 + anorm(mean_vol)^2 + ...
        2*mean_eigs_inners*coords_est);

    corr = inner./(norm_true.*norm_est);

    coords_perf = struct();

    coords_perf.rel_err = rel_err(:)';
    coords_perf.corr = corr;
end
