%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop the `altmin_test_002` script for different values of the dataset size
% (`n`), image sizes (`L`), and noise variances (`noise_var`).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ns = 2.^[10:4:22];
Ls = [16 32];
noise_vars = 3.36e-2*2.^[0:2:10];

[~, filename, ~] = fileparts(mfilename('fullpath'));
output_filename = fullfile(pkg_root(), 'output', [filename '.csv']);

for L_ind = 1:numel(Ls)
    for n_ind = 1:numel(ns)
        for noise_var_ind = 1:numel(noise_vars)
            L = Ls(L_ind);
            n = ns(n_ind);
            im_noise_var = noise_vars(noise_var_ind);

            altmin_test_002;

            f = fopen(output_filename, 'a+');

            data = [L*ones(size(corr, 1), 1) ...
                n*ones(size(corr, 1), 1) ...
                im_noise_var*ones(size(corr, 1), 1) ...
                [0:size(corr, 1)-1]' corr];

            fprintf(f, '%d %d %.8e %d %.8e %.8e %.8e\n', data');

            fclose(f);
        end
    end
end
