% SIM_TO_SRC Create source from simulation
%
% Usage
%    src = sim_to_src(sim);
%
% Input
%    sim: A simulation object obtained from `create_sim`.
%
% Output
%    src: A source object whose images correspond to those of `sim`.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function src = sim_to_src(sim)
    src = struct();

    src.type = src_type_sim();

    src.sim = sim;

    src.L = sim.L;
    src.n = sim.n;
    src.precision = class(sim.vols);

    src.params = sim_to_params(sim);
end
