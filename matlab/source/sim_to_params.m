% SIM_TO_PARAMS Create parameters structure from simulation
%
% Usage
%    params = sim_to_params(sim);
%
% Input
%    sim: A simulation object obtained from `create_sim`.
%
% Output
%    params: A parameters object whose entries correspond to those of `sim`.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function params = sim_to_params(sim)
    params = struct();

    params.rots = sim.rots;
    params.filters = sim.filters;
    params.filter_idx = sim.filter_idx;
    params.offsets = sim.offsets;
    params.amplitudes = sim.amplitudes;
end
