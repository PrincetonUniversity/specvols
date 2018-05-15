% CTF_FILTER Create CTF filter
%
% Usage
%    filter = ctf_filter(ctf_params);
%
% Input
%    ctf_params: A structure containing the parameters of the CTF given by
%          - pixel_size: the pixel size in angstroms,
%          - lambda: the electron wavelength in nm,
%          - defocus_u: the defocus depth along the u-axis in angstrom,
%          - defocus_v: the defocus depth along the v-axis in angstrom,
%          - defocus_ang: the angle between the x-axis and the u-axis in
%             radians,
%          - Cs: the spherical aberration constant,
%          - alpha: the amplitude contrast phase in radians,
%          - B: the envelope decay in inverse square angstrom.
%       (default given by `gen_ctf_params()`)
%
% Output
%    filter: A struct array of filters corresponding to the parameters
%       in ctf_params.
%
% See also
%    gen_ctf_params, eval_filter_ctf

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function filter = ctf_filter(ctf_params)
    if nargin < 1 || isempty(ctf_params)
        ctf_params = gen_ctf_params();
    end

    filter.type = filter_type_ctf();

    filter.radial = false;

    filter.dim = 2;

    filter = repmat(filter, [numel(ctf_params), 1]);

    for n = 1:numel(ctf_params)
        filter(n).ctf_params = ctf_params(n);
        filter(n).radial = ...
            (ctf_params(n).defocus_u == ctf_params(n).defocus_v);
    end
end
