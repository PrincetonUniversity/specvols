% SIM_EIGS Eigendecomposition of volume covariance matrix of simulation
%
% Usage
%    [eigs_true, lambdas_true] = sim_eigs(sim);
%
% Input
%    sim: A simulation object obtained from `create_sim`.
%
% Output
%    eigs_true: The eigenvectors of the volume covariance matrix in the form of
%       an L-by-L-by-L-by-(C-1) array, where C is the number of distinct
%       states in `sim`.
%    lambdas_true: The eigenvalues of the covariance matrix in the form of a
%       (C-1)-by-(C-1) diagonal matrix.
%
% See also
%    sim_covar, sim_mean

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function [eigs_true, lambdas_true] = sim_eigs(sim)
    C = size(sim.vols, 4);

    vols_c = bsxfun(@minus, sim.vols, sim_mean(sim));

    p = ones(1, C)/C;

    vols_c = vol_to_vec(vols_c);

    [Q, R] = qr(vols_c, 0);

    % Rank is at most C-1, so remove last vector.
    Q = Q(:,1:end-1);
    R = R(1:end-1,:);

    [V, lambdas_true] = eig(make_symmat(R*diag(p)*R'));

    eigs_true = vec_to_vol(Q*V);

    [eigs_true, lambdas_true] = mdim_sort_eig(eigs_true, lambdas_true);
end
