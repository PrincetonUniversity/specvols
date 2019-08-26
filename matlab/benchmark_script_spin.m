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
%   November 8, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize variables and settings
% if this script was a function, these would be the inputs
%
% First, we addpath!
%
addpath_het3d;
omp_set_num_threads(max(omp_get_max_threads,16));

disp('Let''s get started!');
disp(datetime('now'));
tic

% Now, the parameters relating to the volume
% 
num_angles = 54;
max_angle = pi/2;

num_shifts_x = 1;
min_shifts_x = 0;
max_shifts_x = 0;

num_shifts_y = 1;
min_shifts_y = 0;
max_shifts_y = 0;

shift_mask = ones(num_shifts_x,num_shifts_y);

% parameters related to the simulation

L = 108;
down_L = 16;
n = 1e5; 
noise_var = 30;    %Frankly, not entirely certain how this is normalized
noise_seed = 0;     %for reproducibility
offsets = zeros(2,n); %currently unused
rots = [];          %[] means uniformly drawn from SO(3)
amplitudes = [];    %for image filters
%states = ceil(cumsum(ones(1,n))/n * num_angles_1 * num_angles_2);
                    %which volumes to use for each image.  They're placed 
                    %together here for ease of debugging
states = ceil(rand(1,n) * num_angles * sum(shift_mask(:)));
filters = identity_filter();
filter_idx = [];    %draw which filter to use uniformly at random

% Parameters related to covariance stuff
% Note! You can cheat on any one of the estimations indivually, allowing
% you to e.g. run a real coordinate estimation with cheating covariance
num_cov_coords = 4;    %number of coordinates to extract from the
                        %covariance PCA thing.  These are used to calculate
                        %the adjacency matrix for the diffusion map

                        % Parameters related to diffusion maps
dmap_t = 0;              %time parameter for diffusion maps
%dists_epsilon = 0.2;      %width parameter for kernel used on the distances
knn_graph_k = 10;

% Parameters related to estimation of diffusion volumes
r = 8;  %also counts for above...


% Parameters related to result checking

% for right now, we check several

disp(['Finished initializing, t = ' num2str(toc)]);
                    
%% Set up test volumes
uncomputed_basis = fb_basis(L*ones(1,3),[]);
if basis_precompute
    basis = precompute_basis(uncomputed_basis);
else
    basis = uncomputed_basis;
end
down_basis = fb_basis(down_L * ones(1,3));

angles = 0:(1/num_angles * max_angle):(max_angle * (num_angles-1)/num_angles);

[x_locs y_locs] = find(shift_mask);
shifts_x = x_locs - ceil(num_shifts_x/2);
shifts_y = y_locs - ceil(num_shifts_y/2);

shifts = [shifts_x(:) shifts_y(:)]';

vol_tops = gen_fakekv_top_rotations(L,angles);
vol_bottoms = gen_fakekv_mid_bottom_shifts(L,shifts);
vols = vol_tops + vol_bottoms;
%
down_vol_tops = gen_fakekv_top_rotations(down_L,angles);
down_vol_bottoms = gen_fakekv_mid_bottom_shifts(down_L,shifts * down_L / L);
down_vols = down_vol_tops + down_vol_bottoms;

disp(['Finished with volumes, t = ' num2str(toc)]);
%% Set up sim object

sim_params = struct();
sim_params = fill_struct(sim_params, ...
        'n', n, ...
        'L', L,  ...
        'vols', vols, ...
        'states', states, ...
        'rots', rots, ...
        'filters', filters, ...
        'filter_idx', filter_idx, ...
        'offsets',offsets, ...
        'amplitudes', amplitudes, ...
        'noise_psd', scalar_filter(noise_var / L^3), ...
        'noise_seed', noise_seed);
sim = create_sim(sim_params);
down_sim = sim;
down_sim.L = down_L;
down_sim.vols = down_vols;


uncached_src = sim_to_src(sim);
src = cache_src(uncached_src);

down_src = (downsample_src(uncached_src,down_L));

disp(['Finished setting up sim, t = ' num2str(toc)]);

%% We do the covariance thing

mean_est_opt = struct();
mean_est_opt = fill_struct(mean_est_opt ,...
        'precision','single',...
        'batch_size',2^13, ...
        'preconditioner','none');
mean_est_opt.max_iter = 500;
mean_est_opt.verbose = 0;
mean_est_opt.iter_callback = [];
mean_est_opt.preconditioner = [];
mean_est_opt.rel_tolerance = 1e-3;
mean_est_opt.store_iterates = true;

mean_vol = estimate_mean(down_src, down_basis, mean_est_opt);

disp(['Finished with mean, t = ' num2str(toc)]);

noise_var_est = estimate_noise_power(down_src);

cov_est_opt = struct();
cov_est_opt = fill_struct(cov_est_opt ,...
        'precision','single',...
        'batch_size',2^13, ...
        'preconditioner','none');
cov_est_opt.max_iter = 50;
cov_est_opt.verbose = 0;
cov_est_opt.iter_callback = [];
cov_est_opt.preconditioner = [];
cov_est_opt.rel_tolerance = 1e-3;
cov_est_opt.store_iterates = true;

covar = estimate_covar(down_src, mean_vol, noise_var_est, down_basis, cov_est_opt);

disp(['Finished with covariance, t = ' num2str(toc)]);

cov_flat = reshape(covar,[down_sim.L^3 down_sim.L^3]);
cov_flat_sym = 1/2 * (cov_flat + cov_flat');
cov_sym_unflat = reshape(cov_flat_sym, size(covar));
[cov_eigs, lambdas] = mdim_eigs(cov_sym_unflat, num_cov_coords, 'la');

disp(['Finished with covar eigs, t = ' num2str(toc)]);

coords = src_wiener_coords(down_src, mean_vol, cov_eigs, lambdas, noise_var_est);

disp(['Finished with covar coords, t = ' num2str(toc)]);

disp(['Finished with all covar stuff, t = ' num2str(toc)]);

%% Diffusion map calculation
%stupid version now, will replace with smarter Amit M code

%dmap_coords = coords_to_laplacian_eigs(coords,dists_epsilon,r);
%W = graph_knn(coords, knn_graph_k);
graph_weights = graph_knn(coords', knn_graph_k);
%graph_weights = graph_gaussian_kernel(coords', dists_epsilon);
graph_laplacian = laplacian(graph_weights, 'normalized');
graph_laplacian = (graph_laplacian + graph_laplacian')/2; % Make sure it is exactly symmetric
[dmap_coords, dmap_evals] = eigs(graph_laplacian, r, 'smallestabs');
dmap_coords = dmap_coords';

disp(['Finished with dmap coords, t = ' num2str(toc)]);

%% "Direct fitting" of diffusion volumes!

unique_states = unique(sim.states); 

unique_state_line_numbers = zeros(size(unique_states));

for unique_state_index = 1:length(unique_states)
    state_places = find(sim.states == unique_states(unique_state_index));
    unique_state_line_numbers(unique_state_index) = state_places(1);
end

rs_used_direct = r;

x = reshape(sim.vols(:,:,:,:),[L^3 numel(unique_states)]); 

alphas = dmap_coords(1:rs_used_direct,unique_state_line_numbers);

%diff_vols_direct_flat = linsolve(alphas',x')';

%diff_vols_direct = reshape(diff_vols_direct_flat,[L L L rs_used_direct]);

disp(['Finished fitting diff vols, t = ' num2str(toc)]);

%% Time to estimate vols!

vols_wt_est_opt = struct();
vols_wt_est_opt = fill_struct(vols_wt_est_opt ,...
        'precision','single',...
        'batch_size',3e4, ...
        'preconditioner','none');
vols_wt_est_opt.max_iter = 50;
vols_wt_est_opt.verbose = 0;
vols_wt_est_opt.iter_callback = [];
vols_wt_est_opt.preconditioner = [];
vols_wt_est_opt.rel_tolerance = 1e-5;
vols_wt_est_opt.store_iterates = true;

%[vols_wt_est, cg_info] = estimate_vols_wt(src, basis, dmap_coords, ...
%    vols_wt_est_opt);

%%Start vols_wt_est_unwrapped

src = src;
basis = basis;
wts = dmap_coords;
vols_wt_est_opt = vols_wt_est_opt;

    L = src.L;
    n = src.n;
    vols_wt_est_opt = fill_struct(vols_wt_est_opt, ...
        'preconditioner', 'none', ...
        'precision', 'single');
    if isempty(basis)
        basis = dirac_basis(L*ones(1, 3));
    end

    kermat_f = sqrt(n^2) * src_vols_wt_kermat(src, wts, vols_wt_est_opt);

    disp(['Finished with kermat, t = ' num2str(toc)]);

    precond_kermat_f = [];

    if ischar(vols_wt_est_opt.preconditioner)
        if strcmp(vols_wt_est_opt.preconditioner, 'none')
            precond_kermat_f = [];
        elseif strcmp(vols_wt_est_opt.preconditioner, 'circulant')
            precond_kermat_f = circularize_kernel(kermat_f, 3);
        else
            error('Invalid preconditioner type.');
        end

        % Reset so this is not used by the `conj_grad` function.
        vols_wt_est_opt.preconditioner = @(x)(x);
    end

    vols_wt_b_coeff = sqrt(n) * src_vols_wt_backward(src, basis, wts, vols_wt_est_opt);
    
    disp(['Finished with backproj!, t = ' num2str(toc)]);

    [vols_wt_est_coeff, cg_info] = conj_grad_vols_wt(kermat_f, vols_wt_b_coeff, ...
        basis, precond_kermat_f, vols_wt_est_opt);

    vols_wt_est = basis_evaluate(basis, vols_wt_est_coeff);

%%End vols_wt_est_unwrapped


disp(['Finished estimating vols, t = ' num2str(toc)]);

disp('All done!');
disp(datetime('now'));
