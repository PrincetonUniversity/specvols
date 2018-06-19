% SIM_NOISE_POWER Noise power of simulation
%
% Usage
%    noise_power = sim_noise_power(sim);
%
% Input
%    sim: A simulation object obtained from `create_sim`.
%
% Output
%    noise_power: The variance of the noise in each pixel of the images from
%       `sim`.
%
% See also
%    sim_signal_power, sim_het_power

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function noise_power = sim_noise_power(sim)
    noise_psd = sim.noise_psd;

    if isnumeric(noise_psd)
        noise_power = noise_psd;
    else
        % TODO: Verify this. Probably not exact, but good enough for now.
        noise_psd = eval_filter_grid(noise_psd, sim.L);
        noise_power = sum(noise_psd(:))/sim.L^2;
    end
end
