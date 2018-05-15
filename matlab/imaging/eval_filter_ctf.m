% EVAL_FILTER_CTF Evaluates a filter of the astigmatic CTF type
%
% Usage
%    h = eval_filter_ctf(filter, omega);
%
% Input
%    filter: A filter structure of `filter_type_ctf`.
%    omega: An array of size 2-by-n representing the spatial frequencies at
%       which the filter is to be evaluated. These are normalized so that
%       pi is equal to the Nyquist frequency.
%
% Output
%    h: The value of the filter at the specified frequencies.
%
% See also
%    ctf_filter, eval_filter

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function h = eval_filter_ctf(filter, omega)
    p = filter.ctf_params;

    om_x = omega(1,:)/(2*pi*p.pixel_size);
    om_y = omega(2,:)/(2*pi*p.pixel_size);

    defocus = zeros(size(om_x), class(om_x));

    ind_nz = (abs(om_x) > eps(pi)) | (abs(om_y) > eps(pi));

    angles_nz = atan2(om_y(ind_nz), om_x(ind_nz));
    angles_nz = angles_nz - p.defocus_ang;

    defocus_mean = (p.defocus_u + p.defocus_v)/2;
    defocus_diff = (p.defocus_u - p.defocus_v)/2;

    defocus(ind_nz) = defocus_mean + defocus_diff*cos(2*angles_nz);

    c2 = -pi/2*2*p.lambda*defocus;
    c4 = pi/2*(p.Cs*1e7)*p.lambda^3;

    r2 = om_x.^2 + om_y.^2;
    r4 = r2.^2;

    gamma = c2.*r2+c4*r4;

    h = sqrt(1-p.alpha^2)*sin(gamma) - p.alpha*cos(gamma);

    if p.B ~= 0
        h = h.*exp(-B*r2);
    end
end
