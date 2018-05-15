% GEN_CTF_PARAMS Generate a set of CTF parameters
%
% Usage
%    ctf_params = gen_ctf_params(ctf_params, defocus);
%
% Input
%    ctf_params: A pre-filled CTF parameter structure with certain values
%       already specified, excluding the defocus field.
%    defocus: A range of defocus values in angstrom (default a uniformly
%       spaced set of 7 values from 1.5e4 to 2.5e4).
%
% Output
%    ctf_params: An array of CTF parameter structures with defocus values set
%       according to the given defocus array and the remainder set according
%       to the given ctf_params structure or filled in with the default
%       values:
%          - pixel_size: 40 angstrom,
%          - lambda: voltage_to_wavelength(200),
%          - Cs: 2.26,
%          - B: 0,
%          - alpha: 0.07.
%       Note that the generated CTFs are radially symmetric (no astigmatism),
%       so the `defocus_u` and `defocus_v` values are identical and
%       `defocus_ang` is set to zero.
%
% See also
%    ctf_filter

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function ctf_params = gen_ctf_params(ctf_params, defocus)
    if nargin < 1
        ctf_params = [];
    end

    if nargin < 2 || isempty(defocus)
        defocus = linspace(1.5e4, 2.5e4, 7);
    end

    ctf_params = fill_struct(ctf_params, ...
        'pixel_size', 40, ...
        'lambda', voltage_to_wavelength(200), ...
        'defocus_u', [], ...
        'defocus_v', [], ...
        'defocus_ang', [], ...
        'Cs', 2.26, ...
        'B', 0, ...
        'alpha', 0.07);

    ctf_params = repmat(ctf_params, [numel(defocus) 1]);

    defocus = num2cell(defocus);
    [ctf_params.defocus_u] = deal(defocus{:});
    [ctf_params.defocus_v] = deal(defocus{:});
    [ctf_params.defocus_ang] = deal(0);
end
