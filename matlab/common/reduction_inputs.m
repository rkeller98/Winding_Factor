function [currents, flux, groups, excitation] = reduction_inputs(data)
% Derive numerical inputs from measured channels and nominal excitation.
currents = [data.id_meas, data.iq_meas, data.ie_meas];
flux = [data.psi_d, data.psi_q];
[excitation, ~, groups] = unique(data.ie_ref);
end
