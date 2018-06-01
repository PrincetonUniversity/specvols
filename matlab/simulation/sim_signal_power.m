% SIM_SIGNAL_POWER Signal power of simulation
%
% Usage
%    signal_power = sim_signal_power(sim);
%
% Input
%    sim: A simulation object obtained from `create_sim`.
%
% Output
%    signal_power: The variance of the clean images in each pixel of the
%       images from `sim`.
%
% See also
%    sim_het_power, sim_noise_power

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function signal_power = sim_signal_power(sim)
    signal_power = 0;

    src = sim_to_src(sim);

    for k = 1:size(sim.vols, 4)
        vol = sim.vols(:,:,:,k);

        src_k = subset_src(src, (sim.states == k));

        mean_kernel_f = src_mean_kernel(src_k);

        signal_power_k = ainner(vol, conv_vol(vol, mean_kernel_f))/sim.L^2;

        prop = sum(sim.states == k)/sim.n;
        signal_power = signal_power + signal_power_k*prop;
    end
end
