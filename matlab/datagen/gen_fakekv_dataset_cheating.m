function coords_true = gen_fakekv_dataset_cheating(L, n, max_angle_1, num_angles_1, max_angle_2, num_angles_2)
    normed_vols = gen_fakekv_volumes_slow(L, max_angle_1, num_angles_1, max_angle_2, num_angles_2);

    %% First, we generate our sim!

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
    %          viewing directions (default generated using 'rand_rots'),
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
    %       - noise_var: the variance of the noise in the images (default 1).

    sim_params = struct();
    sim_params = fill_struct(sim_params, ...
            'n', n, ...
            'L', L,  ...
            'vols', normed_vols, ...
            'states', [], ...
            'rots', [], ...
            'filters', identity_filter(), ...
            'filter_idx', [], ...
            'offsets',zeros(2,n), ...
            'amplitudes', [], ...
            'noise_seed', 0, ...
            'noise_var', 0.1); 

    %fix states for practice!  This makes it easier to debug.
    %Default behavior (achieved by skipping this command) is to draw each
    %volume uniformly at random.  This "cheating" method will probably not work
    %super well as num_angles_1 * num_angles_2 gets big (though it should be
    %fine with 8 of each at n = 511, for example)
%     sim_params.states = ceil(cumsum(ones(1,n))/n * num_angles_1 * num_angles_2);
    
    sim = create_sim(sim_params);

    % Normalize noise_var so that the number we put in is normalized by the
    % amount of noise in each pixel.  Otherwise, 'noise_var' needs to be chosen
    % according to the resolution and scaling of the volumes

    clean_images = sim_clean_image(sim,1,n);
    sim_params.noise_var = sim_params.noise_var * ...
        (( sum(clean_images(:).^2)) / (sim_params.n * sim_params.L^2)) ;
    sim.noise_var = sim_params.noise_var;

    src = sim_to_src(sim);

    %% Cheating params!
    mean_true = sim_mean(sim);
    covar_true = sim_covar(sim);
    [eigs_true, lambdas_true] = sim_eigs(sim);
    [coords_true_states, residuals] = sim_vol_coords(sim, mean_true, eigs_true);
    coords_true = coords_true_states(:,sim.states);

    disp(['Done initializing everything and cheatin, t = ' num2str(toc)]);

    % %% Maybe we want to do it without cheating?
    % 
    % mean_est_opt = struct();
    % mean_est_opt = fill_struct(mean_est_opt ,...
    %         'precision','double',...
    %         'batch_size',512, ...
    %         'preconditioner','none');
    % mean_est_opt.max_iter = 5;
    % mean_est_opt.verbose = 0;
    % mean_est_opt.iter_callback = [];
    % mean_est_opt.preconditioner = [];
    % mean_est_opt.rel_tolerance = 1e-5;
    % mean_est_opt.store_iterates = true;
    % 
    % [mean_est, mean_est_cg_info] = estimate_mean(src);
    % 
    % disp(['Done estimating mean, t = ' num2str(toc)]);
    % 
    % %I'm pretty sure that noise_est doesn't work right, but I also don't think
    % %it matters much
    % noise_est_opt = struct();
    % noise_est_opt = fill_struct(noise_est_opt, ...
    %         'batch_size', 512, ...
    %         'bg_radius', 1, ...
    %         'noise_type', 'white');
    % noise_est = estimate_noise_power(src,noise_est_opt);
    % 
    % 
    % covar_est_opt = struct();
    % covar_est_opt = fill_struct(covar_est_opt ,...
    %         'precision','double',...
    %         'batch_size',512, ...
    %         'preconditioner','none');
    % covar_est_opt.max_iter = 5;
    % covar_est_opt.verbose = 0;
    % covar_est_opt.iter_callback = [];
    % covar_est_opt.preconditioner = [];
    % covar_est_opt.rel_tolerance = 1e-5;
    % covar_est_opt.store_iterates = true;
    % 
    % [covar_est, covar_est_cg_info] = estimate_covar(src, mean_est, noise_est, ...
    %     basis, covar_est_opt);
    % 
    % disp(['Done estimating covar, t = ' num2str(toc)]);
    % 
    % [eigs_est_raw, lambdas_est_raw] = eig(reshape(covar_est,[N^3 N^3]));
    % 
    % [eigs_est_all lambdas_est_all] = mdim_sort_eig(vec_to_vol(eigs_est_raw), lambdas_est_raw);
    % 
    % eigs_est = eigs_est(:,:,:,2:(num_angles_1 * num_angles_2));
    % lambdas_est = lambdas_est(2:(num_angles_1 * num_angles_2),2:(num_angles_1 * num_angles_2));
    % 
    % coords_est = src_wiener_coords(src, mean_est, eigs_est, lambdas_est, noise_est);
    % 
    % disp(['Done estimating coords, t = ' num2str(toc)]);
    % 
    % %% Evaluate!
    % mean_perf = sim_eval_mean(sim, mean_est)
    % covar_perf = sim_eval_covar(sim, covar_est)
    % eigs_perf = sim_eval_eigs(sim, eigs_est, lambdas_est)
    % coords_perf = sim_eval_coords(sim, mean_vol, eig_vols, coords_est)
end
