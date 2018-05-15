% WAVELENGTH_TO_VOLTAGE Convert from wavelength to electron voltage
%
% Usage
%    voltage = wavelength_to_voltage(lambda);
%
% Input
%    lambda: The electron wavelength in nanometers.
%
% Output
%    voltage: The electron voltage in kV.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function voltage = wavelength_to_voltage(lambda)
    voltage = ...
        (-1e3 + sqrt(1e6 + 4*12.2643247^2*0.978466/lambda^2))/(2*0.978466);
end
