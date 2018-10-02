% SIM_VOL_COORDS Coordinates of simulation volumes in a given basis
%
% Usage
%    [coords, residuals] = sim_vol_coords(sim, mean_vol, eig_vols);
%
% Input
%    sim: Simulation object from `create_sim`.
%    mean_vol: A mean volume in the form of an L-by-L-by-L array (default
%       `sim_mean(sim)`).
%    eig_vols: A set of eigenvolumes in an L-by-L-by-L-by-K array (default
%       `sim_eigs(sim)`).
%
% Output
%    coords: A vector of size K-by-C containing the coordinates of the
%       simulation volumes of `sim`, where C is equal to `size(sim.vols, 4)`.
%       These coordinates are in the affine space centered in `mean_vol`, with
%       variability along the direction of `eig_vols`.
%    residuals: A vector of length C containing the norm of the residual after
%       projection onto the affine space.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function [coords, residuals] = sim_vol_coords(sim, mean_vol, eig_vols)
    if nargin < 2 || isempty(mean_vol)
        mean_vol = sim_mean(sim);
    end

    if nargin < 3 || isempty(eig_vols)
        [eig_vols, ~] = sim_eigs(sim);
    end

    vols = bsxfun(@minus, sim.vols, mean_vol);

    coords = vol_to_vec(eig_vols)'*vol_to_vec(vols);
    residuals = anorm(vols - vec_to_vol(vol_to_vec(eig_vols)*coords), 1:3);

    residuals = permute(residuals, [5 4 1:3]);
end
