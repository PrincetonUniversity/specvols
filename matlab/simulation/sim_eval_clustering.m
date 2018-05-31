% SIM_EVAL_CLUSTERING Evaluate clustering accuracy
%
% Usage
%    clustering_perf = sim_eval_clustering(sim, idx);
%
% Input
%    sim: Simulation object from `create_sim`.
%    idx: The estimated states of each image in the simulation.
%
% Output
%    clustering_perf: A struct containing the evaluation results. It contains
%       the fields:
%          - accuracy: The proportion of correctly assigned images.
%          - state_perm: The permutation applied to `idx` in order to match up
%             with the state labeling in `sim.states`.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function clustering_perf = sim_eval_clustering(sim, idx)
    state_perms = perms(1:max(idx));

    accuracy = 0;
    state_perm = zeros(1, size(state_perms, 2));
    for l = 1:size(state_perms, 1)
        accuracy_l = sum(state_perms(l,idx) == sim.states);

        if accuracy_l > accuracy
            accuracy = accuracy_l;
            state_perm = state_perms(l,:);
        end
    end

    accuracy = accuracy/sim.n;

    clustering_perf = struct();

    clustering_perf.accuracy = accuracy;
    clustering_perf.state_perm = state_perm;
end
