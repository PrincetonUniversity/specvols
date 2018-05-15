% SIM_COVAR Volume covariance matrix of simulation
%
% Usage
%    covar_true = sim_covar(sim);
%
% Input
%    sim: A simulation object obtained from `create_sim`.
%
% Output
%    covar_true: The volume covariance matrix of the simulation in the form of
%       an L-by-L-by-L-by-L-by-L-by-L volume matrix.
%
% See also
%    sim_mean

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function covar_true = sim_covar(sim)
    C = size(sim.vols, 4);

    vols_c = bsxfun(@minus, sim.vols, sim_mean(sim));

    p = ones(1, C)/C;

    vols_c = vol_to_vec(vols_c);

    covar_true = vols_c*diag(p)*vols_c';

    covar_true = vecmat_to_volmat(covar_true);
end
