% CREATE_SIM Create cryo-EM simulation
%
% Usage
%    sim = create_sim(sim_params);
%
% Input
%    sim_params: A partial simulation structure whose fields are used to create
%       a full simulation structure (default empty). This can be used to
%       specify certain parameters (such as n, L or vol) while letting the
%       others be automatically generated.
%
% Output
%    sim: A simulation structure containing the fields:
%       - n: the number of images in the problem (default 1024),
%       - L: the resolution, that is L-by-L images and L-by-L-by-L volumes
%          (default 8),
%       - vols: An L-by-L-by-L-by-C array of volumes representing the C
%          volumes in the simulation (default `gaussian_blob_vols`),
%       - states: a 1-by-n array containing the hvolume states for each
%          images (default randomly sampled uniformly between 1 and K),
%       - rots: a 3-by-3-by-n array of rotation matrices corresponding to
%          viewing directions (default generated using 'unif_rand_rots'),
%       - filters: a struct array of filter F objects (see `eval_filter`)
%          (default `identity_filter()`),
%       - filter_idx: a 1-by-n array containing the filter function assigments
%          for each of the images (default randomly sampled uniformly between
%          1 and F),
%       - offsets: a 2-by-n array specifying the shifts of the images (default
%          generated from a Gaussian distribution of standard deviation
%          L/16),
%       - amplitudes: a 1-by-n array specifying the amplitude multipliers of
%          the images (default uniformly sampled between 2/3 and 3/2),
%       - noise_seed: the random seed for generating the noise in the images
%           (default 0), and
%       - noise_psd: the power spectral distribution of the noise in the
%           images in the form of a filter object (default
%           `scalar_filter(1)`).

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function sim = create_sim(sim)
    if nargin < 1
        sim = struct();
    end

    sim = fill_struct(sim, ...
        'n', 1024, ...
        'L', 8,  ...
        'vols', [], ...
        'states', [], ...
        'rots', [], ...
        'filters', identity_filter(), ...
        'filter_idx', [], ...
        'offsets', [], ...
        'amplitudes', [], ...
        'noise_seed', 0, ...
        'noise_psd', scalar_filter(1));

    if isempty(sim.vols)
        rand_push(0);

        sim.vols = gaussian_blob_vols(sim.L);

        rand_pop();
    else
        sim.L = size(sim.vols, 1);
    end

    if isempty(sim.states)
        rand_push(0);

        sim.states = randi(size(sim.vols, 4), 1, sim.n);

        rand_pop();
    end

    if isempty(sim.rots)
        rand_push(0);

        sim.rots = unif_rand_rots(sim.n);

        rand_pop();
    end

    if isempty(sim.filter_idx)
        rand_push(0);

        sim.filter_idx = randi(numel(sim.filters), 1, sim.n);

        rand_pop();
    end

    if isempty(sim.offsets)
        rand_push(0);

        sim.offsets = sim.L/16*randn(2, sim.n);

        rand_pop();
    end

    if isempty(sim.amplitudes)
        rand_push(0);

        mn = 2/3;
        mx = 3/2;

        sim.amplitudes = mn + rand(1, sim.n)*(mx - mn);

        rand_pop();
    end
end
