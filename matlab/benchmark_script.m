%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  __   __       __                 __      
% |__) |_  |\ | /   |__| |\/|  /\  |__) |_/ 
% |__) |__ | \| \__ |  | |  | /--\ | \  | \ 
% 
%          __  __  __     __  ___ 
%         (_  /   |__) | |__)  |  
%         __) \__ | \  | |     |  
% 
% A script to perform benchmarks for reconstructing cryo-EM with
% heterogeneity by the method of diffusion volumes.
%
% Outline of script:
%   1. Initialize variables and settings
%   2. Set up test volumes
%   3. Set up sim object
%   4. Compute mean volume, covariance, covariance coordinates
%   5. Computer diffusion map (laplacian eigenvectors)
%   6. Estimate diffusion volumes
%   7. Check results
%   8. Outputs
%
%   Amit Halevi (ahalevi@princeton.edu)
%   September 29, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize variables and settings
% if this script was a function, these would be the inputs
%
% First, we addpath!

addpath_cov3d;

disp('Let''s get started!');
disp(datetime('now'));
tic

% Now, the parameters relating to the volume
% 
num_angles_1 = 32;
max_angle_1 = pi/2;
num_angles_2 = 1;
max_angle_2 = 0;
basis_type = 'dirac';
basis_precompute = 1;
save_vols = 0;      %save precalculated volumes to save time
save_vol_name = 'pregen_vols';  %will have N and num_angles
                                            %added to the name
load_vols = 0;      %load precalculated volumes if present
load_vol_name = []; %default uses save_vol_name
                                            
% parameters related to the simulation

% N = 54;
% n = 1e5;
% noise_var = 0.1;      %normalized to the average energy of the volumes used
noise_seed = 0;     %for reproducibility
offsets = zeros(2,n); %currently unused
rots = [];          %[] means uniformly drawn from SO(3)
amplitudes = [];    %for image filters
states = ceil(cumsum(ones(1,n))/n * num_angles_1 * num_angles_2);
                    %which volumes to use for each image.  They're placed 
                    %together here for ease of debugging
filters = identity_filter();
filter_idx = [];    %draw which filter to use uniformly at random

% Parameters related to covariance stuff
% Note! You can cheat on any one of the estimations indivually, allowing
% you to e.g. run a real coordinate estimation with cheating covariance
num_cov_coords = 16;    %number of coordinates to extract from the
                        %covariance PCA thing.  These are used to calculate
                        %the adjacency matrix for the diffusion map
mean_cheat = 1;
cov_cheat = 1;
eigs_cheat = 1;
coords_cheat = 1;

% Parameters related to diffusion maps
dmap_t = 0;              %time parameter for diffusion maps
dists_epsilon = 0.2;      %width parameter for kernel used on the distances

% Parameters related to estimation of diffusion volumes
% r = 6;  %also counts for above...
save_results = 1;
save_results_name = 'results';

% Parameters related to result checking

% for right now, we check several

disp(['Finished initializing, t = ' num2str(toc)]);
                    
%% Set up test volumes
uncomputed_basis = dirac_basis(N*ones(1,3),[]);
if basis_precompute
    basis = precompute_basis(uncomputed_basis);
else
    basis = uncomputed_basis;
end

if load_vols
    if exist(load_vol_name)
        full_load_name = [load_vol_name '_N_' num2str(N) '_num_ang_1_' ...
            num2str(num_angles_1) '_num_ang_2_' num2str(num_angles_2) '.mat'];
    else
        full_load_name = [save_vol_name '_N_' num2str(N) '_num_ang_1_' ...
            num2str(num_angles_1) '_num_ang_2_' num2str(num_angles_2) '.mat'];
    end
    if exist(full_load_name,'file');
        load(full_load_name)
    else
        vols = generate_vols(N,num_angles_1,num_angles_2,basis);
    end
end

if save_vols
    full_save_name = [save_vol_name '_N_' num2str(N) '_num_ang_1_' ...
        num2str(num_angles_1) '_num_ang_2_' num2str(num_angles_2) '.mat'];
    if ~exist(save_vol_name,'file')
        save(full_save_name','vols');
    end
end

disp(['Finished with volumes, t = ' num2str(toc)]);
%% Set up sim object

sim_params = struct();
sim_params = fill_struct(sim_params, ...
        'n', n, ...
        'L', N,  ...
        'vols', vols, ...
        'states', states, ...
        'rots', rots, ...
        'filters', filters, ...
        'filter_idx', filter_idx, ...
        'offsets',offsets, ...
        'amplitudes', amplitudes, ...
        'noise_seed', noise_seed, ...
        'noise_var', noise_var); 
sim = create_sim(sim_params);

%Normalize noise
clean_images = sim_clean_image(sim,1,n);
sim_params.noise_var = sim_params.noise_var * ...
    (( sum(clean_images(:).^2)) / (sim_params.n * sim_params.L^2)) ;
sim.noise_var = sim_params.noise_var;
uncached_src = sim_to_src(sim);
src = cache_src(uncached_src);

disp(['Finished setting up sim, t = ' num2str(toc)]);

%% We do the covariance thing

if mean_cheat
    mean_vol = sim_mean(sim);
else
    mean_est_opt = struct();
    mean_est_opt = fill_struct(mean_est_opt ,...
            'precision','single',...
            'batch_size',512, ...
            'preconditioner','none');
    mean_est_opt.max_iter = 500;
    mean_est_opt.verbose = 0;
    mean_est_opt.iter_callback = [];
    mean_est_opt.preconditioner = [];
    mean_est_opt.rel_tolerance = 1e-5;
    mean_est_opt.store_iterates = true;
    
    mean_vol = estimate_mean(src, basis, mean_est_opt);
end

disp(['Finished with mean, t = ' num2str(toc)]);


if cov_cheat
    covar = sim_covar(sim);
else

    noise_var_est = estimate_noise_power(src);
    
    cov_est_opt = struct();
    cov_est_opt = fill_struct(cov_est_opt ,...
            'precision','single',...
            'batch_size',512, ...
            'preconditioner','none');
    cov_est_opt.max_iter = 50;
    cov_est_opt.verbose = 0;
    cov_est_opt.iter_callback = [];
    cov_est_opt.preconditioner = [];
    cov_est_opt.rel_tolerance = 1e-5;
    cov_est_opt.store_iterates = true;
    
    cov = estimate_covar(src, mean_vol, noise_var_est, basis, cov_est_opt);
end

disp(['Finished with covariance, t = ' num2str(toc)]);

if eigs_cheat
    [cov_eigs, lambdas] = sim_eigs(sim);
else
    [cov_eigs, lambdas] = mdim_eigs(covar_est, num_cov_coords, 'la');
end

disp(['Finished with covar eigs, t = ' num2str(toc)]);

if coords_cheat
    [eigs_true, lambdas_true] = sim_eigs(sim);
    [coords_true_states, residuals] = sim_vol_coords(sim, mean_vol, eigs_true);
    coords = coords_true_states(:,sim.states);
else
    coords = src_wiener_coords(src, mean_vol, cov_eigs, lambdas, ...
        noise_var_est);
end

disp(['Finished with covar coords, t = ' num2str(toc)]);

disp(['Finished with all covar stuff, t = ' num2str(toc)]);

%% Diffusion map calculation
%stupid version now, will replace with smarter Amit M code

dmap_coords = coords_to_laplacian_eigs(coords,dists_epsilon,r);

disp(['Finished with dmap coords, t = ' num2str(toc)]);

%% Time to estimate vols!

vols_wt_est_opt = struct();
vols_wt_est_opt = fill_struct(vols_wt_est_opt ,...
        'precision','single',...
        'batch_size',3.4e4, ...
        'preconditioner','none');
vols_wt_est_opt.max_iter = 2000;
vols_wt_est_opt.verbose = 0;
vols_wt_est_opt.iter_callback = [];
vols_wt_est_opt.preconditioner = [];
vols_wt_est_opt.rel_tolerance = 1e-5;
vols_wt_est_opt.store_iterates = true;

[vols_wt_est, cg_info] = estimate_vols_wt(src, basis, dmap_coords, ...
    vols_wt_est_opt);


disp(['Finished estimating vols, t = ' num2str(toc)]);

%% Check results!

mse_err = check_recon_wts(uncached_src,dmap_coords,vols_wt_est,'mse');
corr_err = check_recon_wts(uncached_src,dmap_coords,vols_wt_est,'corr');
fsc_err = check_recon_wts(uncached_src,dmap_coords,vols_wt_est,'fsc');

disp(['Finished calculating error, t = ' num2str(toc)]);

if save_results
    full_save_name = [save_vol_name '_N_' num2str(N) '_n_' ...
        num2str(n) 'r' num2str(r) '.mat'];
    save(full_save_name,'uncached_src','dmap_coords','vols_wt_est');
end

disp('All done!');
disp(datetime('now'));