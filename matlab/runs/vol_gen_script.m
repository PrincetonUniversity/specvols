disp('Let''s get started!');
disp(datetime('now'));
tic

max_threads = omp_get_max_threads;
threads_to_use = min(max_threads,16);
omp_set_num_threads(threads_to_use);

%Size, important to both volumes and simulation
L = 108;
down_L = 16;

% Now, the parameters relating to the volumes
%
num_angles_1 = L;
max_angle_1 = pi/2;

num_shifts_x = 1;
min_shifts_x = 0;
max_shifts_x = 0;
num_shifts_y = 1;
min_shifts_y = 0;
max_shifts_y = 0;

% parameters related to the simulation

n = 1e4;
noise_var = 0;    %Frankly, not entirely certain how this is normalized
noise_seed = 0;     %for reproducibility
offsets = zeros(2,n); %currently unused
rots = [];          %[] means uniformly drawn from SO(3)
amplitudes = [];    %for image filters
states = ceil(rand(1,n) * num_angles_1 * num_shifts_x * num_shifts_y);
ctf_params = gen_ctf_params();
filters = ctf_filter(ctf_params);
filter_idx = [];    %draw which filter to use uniformly at random

uncomputed_basis = dirac_basis(L*ones(1,3),[]);
basis = precompute_basis(uncomputed_basis);
down_basis = fb_basis(down_L * ones(1,3));

angles = 0:max_angle_1/num_angles_1:(max_angle_1 - max_angle_1/num_angles_1);

shift_mask = ones(num_shifts_x,num_shifts_y);

[x_locs y_locs] = find(shift_mask);
shifts_x = x_locs - ceil(num_shifts_x/2);
shifts_y = y_locs - ceil(num_shifts_y/2);

shifts = [shifts_x(:) shifts_y(:)]';


vol_tops = gen_fakekv_top_rotations(L,angles);
vol_bottoms = gen_fakekv_mid_bottom_shifts(L,shifts);
vols = vol_tops + vol_bottoms;

down_vol_tops = gen_fakekv_top_rotations(down_L,angles);
down_vol_bottoms = gen_fakekv_mid_bottom_shifts(down_L,shifts * down_L / L);
down_vols = down_vol_tops + down_vol_bottoms;
%
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
