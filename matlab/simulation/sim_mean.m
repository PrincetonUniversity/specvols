% SIM_MEAN Mean volume of simulation
%
% Usage
%    mean_true = sim_mean(sim);
%
% Input
%    sim: A simulation object obtained from `create_sim`.
%
% Output
%    mean_true: The mean volume of the simulation in the form of an
%       L-by-L-by-L array.
%
% See also
%    sim_covar

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function mean_true = sim_mean(sim)
    mean_true = mean(sim.vols, 4);
end
