% SIM_SIGNAL_SNR Signal-to-noise ratio of simulation
%
% Usage
%    signal_snr = sim_signal_snr(sim);
%
% Input
%    sim: A simulation object obtained from `create_sim`.
%
% Output
%    signal_snr: The signal-to-noise ratio of the images in `sim`.
%
% See also
%    sim_signal_power, sim_noise_power

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function signal_snr = sim_signal_snr(sim)
    signal_snr = sim_signal_power(sim)/sim_noise_power(sim);
end
