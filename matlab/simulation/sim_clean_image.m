% SIM_CLEAN_IMAGE Extract clean images from simulation
%
% Usage
%    im = sim_clean_image(sim, start, num);
%
% Input
%    sim: Simulation object from `create_sim`.
%    start: First index of the images to extract.
%    num: The number of images to extract.
%
% Output
%    im: The clean simulation images from start to start+num-1.
%
% See also
%    sim_image, create_sim, vol_project, im_filter, im_translate

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function im = sim_clean_image(sim, start, num)
    im = zeros([sim.L*ones(1, 2) num]);

    all_idx = start:start+num-1;

    unique_states = unique(sim.states(all_idx));

    for k = unique_states
        vol_k = sim.vols(:,:,:,k);

        idx_k = find(sim.states(all_idx) == k);

        im_k = vol_project(vol_k, sim.rots(:,:,all_idx(idx_k)));

        im(:,:,idx_k) = im_k;
    end

    unique_filters = unique(sim.filter_idx(all_idx));

    for k = unique_filters
        idx_k = find(sim.filter_idx(all_idx) == k);

        im(:,:,idx_k) = ...
            im_filter(im(:,:,idx_k), sim.filters(k));
    end

    im = im_translate(im, sim.offsets(:,all_idx));

    im = bsxfun(@times, im, permute(sim.amplitudes(all_idx), [1 3 2]));
end
