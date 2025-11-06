function [Voigt_p,Der] = VoigtProfileV3(W_array,param)
% This code generate a voigt profile based on the given parameters using
% error functions.

    % Individual parameters can be any size matrix in a single row
    peaks = param(1,:); % wavelength peaks
    dopwidth = param(2,:); % doppler half width
    lozwidth = param(3,:); % lorentzian half width
    
    % Complex error function equations to get voigt profile
    v = W_array*ones(1,length(peaks)); % wavenumber array for each peak
    v0 = ones(length(W_array),1)*peaks; % Peak wavelengths
    x = (2*sqrt(log(2))).*(v-v0)./(ones(length(W_array),1)*dopwidth);
    y = ones(length(W_array),1)*(lozwidth./dopwidth)*(sqrt(log(2)));
    z = x + 1i*y;
    w_z = fadf(z); % Use faddeeva function code from Sanjar M. Abrarov and Brendan M. Quine
    Voigt_p = 2*real(w_z)*sqrt(log(2))/(sqrt(pi)*dopwidth);
    
    % Finding Derivative for fitting function (This section is for fitting
    % solver only, helps increase solver accuracy)
    K = real(w_z);
    L = imag(w_z);
    dV_dw=2*(x.*y.*K-y.^2.*L)./(ones(length(W_array),1)*lozwidth)...
        .*(ones(length(W_array),1));
    dV_dDop=2*((x.^2-y.^2).*K-2*x.*y.*L+y./sqrt(pi))./(ones(length(W_array),1)*dopwidth)...
        .*(ones(length(W_array),1));
    dV_dLoz=2*(x.*y.*L+y.^2.*K-y./sqrt(pi))./(ones(length(W_array),1)*lozwidth)...
        .*(ones(length(W_array),1));
    
    if nargout > 1
        der = [dV_dw;dV_dDop;dV_dLoz]; % This is also the jacobian matrix
        Der = reshape(der,[length(W_array) numel(param)]);
    end
end