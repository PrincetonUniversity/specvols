% SIM_IMAGE Extract images from simulation
%
% Usage
%    im = sim_image(sim, start, num);
%
% Input
%    sim: Simulation object from `create_sim`.
%    start: First index of the images to extract.
%    num: The number of images to extract.
%
% Output
%    im: The simulation images from start to start+num-1.
%
% See also
%    sim_clean_image, sim_noise_image, create_sim

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function im = sim_image(sim, start, num)
    im = sim_clean_image(sim, start, num) + sim_noise_image(sim, start, num);
end
