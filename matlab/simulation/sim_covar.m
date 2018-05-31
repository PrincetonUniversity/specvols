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
%    sim_eigs, sim_mean

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function covar_true = sim_covar(sim)
    [eigs_true, lambdas_true] = sim_eigs(sim);

    eigs_true = vol_to_vec(eigs_true);

    covar_true = eigs_true*lambdas_true*eigs_true';

    covar_true = vecmat_to_volmat(covar_true);
end
