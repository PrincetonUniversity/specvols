% ESTIMATE_NOISE_POWER Estimate noise variance from source
%
% Usage
%    noise_var_est = estimate_noise_power(src, noise_est_opt);
%
% Input
%    src: A source structure containing the images whose noise variance is to
%       be estimated.
%    noise_est_opt: An options struct for `estimate_noise`. All fields are
%       passed on to `estimate_noise` as they are, except `noise_type`, which
%       is set to 'white'.
%
% Output
%    noise_var_est: The estimated noise variance of the images.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function noise_var_est = estimate_noise_power(src, noise_est_opt)
    if nargin < 2 || isempty(noise_est_opt)
        noise_est_opt = struct();
    end

    noise_est_opt.noise_type = 'white';

    noise_psd_est = estimate_noise(src, noise_est_opt);

    noise_var_est = eval_filter(noise_psd_est, zeros(2, 1));
end