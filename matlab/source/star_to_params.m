% STAR_TO_PARAMS Create parameters structure from a STAR file
%
% Usage
%    params = star_to_params(star);
%
% Input
%    star: A STAR file struct obtained from `load_star`.
%    star_opt: An options structure with the fields:
%          - pixel_size: the pixel size of the images in angstroms (default 1),
%          - B: the envelope decay of the CTF in inverse square angstrom
%             (default 0),
%       These parameters are not given in the STAR file itself and so need to
%       be supplied.
%
% Output
%    params: A parameters object whose entries correspond to those of the
%       STAR file.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function params = star_to_params(star, star_opt)
    if nargin < 2 || isempty(star_opt)
        star_opt = struct();
    end

    star_opt = fill_struct(star_opt, ...
        'pixel_size', 1, ...
        'B', 0);

    if ~isstruct(star)
        error(['Input `star` must be a STAR file struct obtained from ' ...
            '`load_star`']);
    end

    if isfield(star, 'root')
        star = star.root;
    elseif isfield(star, 'images')
        star = star.images;
    elseif ~isfield(star, 'rlnImageName')
        error('Invalid STAR file format.');
    end

    angles_rot = [star.rlnAngleRot];
    angles_tilt = [star.rlnAngleTilt];
    angles_psi = [star.rlnAnglePsi];

    origin_x = [star.rlnOriginX];
    origin_y = [star.rlnOriginY];

    class = [star.rlnClassNumber];

    n = numel(star);

    angles_rot = angles_rot/180*pi;
    angles_tilt = angles_tilt/180*pi;
    angles_psi = angles_psi/180*pi;

    angles = [angles_rot(:)'; angles_tilt(:)'; angles_psi(:)'];

    [ctf_params, ctf_idx] = star_ctf_params(star, star_opt);

    params = struct();

    params.rots = angles_to_rots(angles);

    params.filters = ctf_filter(ctf_params);
    params.filter_idx = ctf_idx;

    params.offsets = [origin_x(:)'; origin_y(:)'];

    params.amplitudes = ones(1, n);

    params.state = class;
end

function [ctf_params, ctf_idx] = star_ctf_params(star, star_opt)
    is_astigmatic = isfield(star, 'rlnDefocusV');

    if isfield(star, 'rlnDefocusV');
        defocus = [[star.rlnDefocusU]' [star.rlnDefocusV]' ...
            [star.rlnDefocusAngle]'];
    else
        defocus = [[star.rlnDefocusU]' [star.rlnDefocusU]' ...
            zeros(numel(star), 1)];
    end

    params = [[star.rlnVoltage]' defocus ...
        [star.rlnSphericalAberration]' [star.rlnAmplitudeContrast]'];

    [params, ~, ctf_idx] = unique(params, 'rows');

    ctf_params = struct();

    for k = 1:size(params, 1)
        ctf_params(k).pixel_size = star_opt.pixel_size;
        ctf_params(k).lambda = voltage_to_wavelength(params(k,1));

        ctf_params(k).defocus_u = params(k,2);
        ctf_params(k).defocus_v = params(k,3);
        ctf_params(k).defocus_ang = params(k,4)/360*2*pi;

        ctf_params(k).Cs = params(k,5);
        ctf_params(k).alpha = params(k,6);
        ctf_params(k).B = star_opt.B;
    end
end
