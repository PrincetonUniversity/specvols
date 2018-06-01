% SIM_HET_SNR Heterogeneity signal-to-noise ratio of simulation
%
% Usage
%    het_snr = sim_het_snr(sim);
%
% Input
%    sim: A simulation object obtained from `create_sim`.
%
% Output
%    het_snr: The heterogeneity signal-to-noise ratio of the images in `sim`.
%
% See also
%    sim_het_power, sim_noise_power

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function het_snr = sim_het_snr(sim)
    het_snr = sim_het_power(sim)/sim_noise_power(sim);
end
