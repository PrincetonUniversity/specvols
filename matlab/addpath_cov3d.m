% ADDPATH_COV3D Add paths for the package
%
% Usage
%    addpath_cov3d();

path_to_pkg = fileparts(mfilename('fullpath'));

addpath(fullfile(path_to_pkg, 'basis'));
addpath(fullfile(path_to_pkg, 'estimation'));
addpath(fullfile(path_to_pkg, 'imaging'));
addpath(fullfile(path_to_pkg, 'install'));
addpath(fullfile(path_to_pkg, 'io'));
addpath(fullfile(path_to_pkg, 'nufft'));
addpath(fullfile(path_to_pkg, 'simulation'));
addpath(fullfile(path_to_pkg, 'source'));
addpath(fullfile(path_to_pkg, 'utils'));

if exist(fullfile(path_to_pkg, 'extern', 'nufftall-1.33'))
    addpath(fullfile(path_to_pkg, 'extern', 'nufftall-1.33'));
end

if exist(fullfile(path_to_pkg, 'extern', 'nfft'))
    addpath(fullfile(path_to_pkg, 'extern', 'nfft', 'lib'));
    addpath(fullfile(path_to_pkg, 'extern', 'nfft', 'share', 'nfft', ...
        'matlab', 'nfft'));
end

if exist(fullfile(path_to_pkg, 'extern', 'finufft'))
    addpath(fullfile(path_to_pkg, 'extern', 'finufft', 'matlab'));
end

clear path_to_pkg;
