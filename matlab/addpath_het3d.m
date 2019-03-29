% ADDPATH_HET3D Add paths for the package
%
% Usage
%    addpath_HET3d();

path_to_pkg = fileparts(mfilename('fullpath'));

addpath(fullfile(path_to_pkg, 'basis'));
addpath(fullfile(path_to_pkg, 'estimation'));
addpath(fullfile(path_to_pkg, 'examples'));
addpath(fullfile(path_to_pkg, 'imaging'));
addpath(fullfile(path_to_pkg, 'install'));
addpath(fullfile(path_to_pkg, 'io'));
addpath(fullfile(path_to_pkg, 'nufft'));
addpath(fullfile(path_to_pkg, 'runs'));
addpath(fullfile(path_to_pkg, 'simulation'));
addpath(fullfile(path_to_pkg, 'source'));
addpath(fullfile(path_to_pkg, 'utils'));
addpath(fullfile(path_to_pkg, 'graph'));
addpath(fullfile(path_to_pkg, 'datagen'));
addpath(fullfile(path_to_pkg, 'utils', 'mex'));
addpath(fullfile(path_to_pkg, 'wts'));
addpath(fullfile(path_to_pkg, 'figuresgen'));

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
