% VOLTAGE_TO_WAVELENGTH Convert from electron voltage to wavelength
%
% Usage
%    lambda = voltage_to_wavelength(voltage);
%
% Input
%    voltage: The electron voltage in kV.
%
% Output
%    lambda: The electron wavelength in nanometers.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function lambda = voltage_to_wavelength(voltage)
    lambda = 12.2643247/sqrt(voltage*1e3 + 0.978466*voltage^2);
end
