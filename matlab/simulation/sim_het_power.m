% SIM_HET_POWER Heterogeneity power of simulation
%
% Usage
%    het_power = sim_het_power(sim);
%
% Input
%    sim: A simulation object obtained from `create_sim`.
%
% Output
%    het_power: The variance of the mean-subtracted clean images in each pixel
%       of the images from `sim`. The images are the projections of `sim` with
%       the projections of the mean volume subtracted for each image.
%
% See also
%    sim_signal_power, sim_noise_power

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function het_power = sim_het_power(sim)
    mean_vol = sim_mean(sim);

    sim_centered = sim;

    sim_centered.vols = bsxfun(@minus, sim.vols, mean_vol);

    het_power = sim_signal_power(sim_centered);
end
