% Spectral Simulation function by Marc Etienne
% This function works similar to spectraplot, although it can be used to 
% simulate other conditions other than in air
% The alpha_total is the absorbance contribution of all lines in range and
% the alpha output is the individual absorbance contribution of each line
function [alpha_total,W_array,alpha,alpha_array] = SpectralSim(Spectral_data,T,P,X,L,v_start,v_end,v_step,SimType,Species)
    % Paramaters for calculation and simulation
    if iscell(Spectral_data)
        Spectral_data = Spectral_data{Species};
    else
        Spectral_data = Spectral_data.data.(Species);
    end
    W_array = v_start:v_step:v_end; W_array = W_array';
    M = Spectral_data(1,25); %Molecular weight of absorbing species [g/mol]
    h = 6.62607015e-34; % [m2*kg/s]
    c = 299792458*100; % [cm/s]
    k = 1.380649e-23; % [m2*kg*s-2*K-1]
    Q_0 = Spectral_data(1,26)*296^6+Spectral_data(2,26)*296^5+Spectral_data(3,26)*296^4+...
        Spectral_data(4,26)*296^3+Spectral_data(5,26)*296^2+Spectral_data(6,26)*296+Spectral_data(7,26);
    Q = Spectral_data(1,26)*T^6+Spectral_data(2,26)*T^5+Spectral_data(3,26)*T^4+...
        Spectral_data(4,26)*T^3+Spectral_data(5,26)*T^2+Spectral_data(6,26)*T+Spectral_data(7,26);

    % Determine required transitions based on user select bound
    selected_lines = Spectral_data(:,1) >= v_start & v_end >= Spectral_data(:,1);
    Transitions = Spectral_data(selected_lines,:);
    
    % Determine transition parameters
    % Transition Information
    V_0 = Transitions(:,1); % Peak location of transition [cm-1]
    S = Transitions(:,2).*(Q_0/Q)*(296/T).*exp((-h*c*Transitions(:,8)/k)*((1/T)-(1/296)))...
        .*(1-exp(-h*c*V_0/(k*T))).*(1-exp(-h*c*V_0/(k*296))).^(-1); %[cm-2/atm]
    g_self_0 = Transitions(:,3); %[cm-1/atm]
    n_self = Transitions(:,5);
    d_air_0 = Transitions(:,6);
    % Determine air bath gas parameters
    if strcmp(SimType,'HITEMP Sim')
        g_air_0 = Transitions(:,4); %[cm-1/atm]
        n_air = Transitions(:,7);

        % Determine parameters at gas condition
        [V,g_self,g_air] = ShiftedParam(V_0,g_self_0,g_air_0,P,X,T,d_air_0,n_air,n_air,n_self);
    
        % Doppler and Collisional (Lorentzian) Full Width at Half Maximum Calculation
        [V_d,V_c] = Dop_and_Col_Terms(V,X,g_self,g_air,P,T,M);
    else
        % If specific mixture is needed
        % User select mixture
        prompt_val = inputdlg({'N2','CO2','H2O','OH','H2','O2','O','H'},'Select Gas Product Concentration');
        Mixture = str2double(prompt_val);
        if isempty(Mixture)
            alpha = nan;
            alpha_array = nan;
            alpha_total = nan;
            return % Cancel operation if no value was selected from user
        elseif isnan(Mixture)
            alpha = nan;
            alpha_array = nan;
            alpha_total = nan;
            return % Cancel operation if no value was selected from user
        end
        select_mix = ~isnan(Mixture);
        X_mix = Mixture(select_mix);

        % Determine the corresponding column for each bath species
        g_options = [9,11,13,15,17,19,21,23]; g_options = g_options(select_mix);
        n_options = [10,12,14,16,18,20,22,24]; n_options = n_options(select_mix);

        % Pre-allocate matrix size to prevent running into size issues
        g_mix_0 = zeros(size(Transitions,1),size(g_options,2));
        n_mix = zeros(size(Transitions,1),size(g_options,2));

        % Determine coefficients for bath mixture
        g_mix_0(:,1) = Transitions(:,4); %[cm-1/atm]
        n_air(:,1) = Transitions(:,7);
        n_mix(:,1) = n_air;

        % Find lines with species parameters other than for air
        isnan_L = ~isnan(Transitions(:,g_options(1)));
        S_L = find(isnan_L);

        for i = S_L
            g_mix_0(i,:) = Transitions(i,g_options);
            n_mix(i,:) = Transitions(i,n_options);
        end

        % Determine parameters at gas condition
        [V,g_self,g_mix] = ShiftedParam(V_0,g_self_0,g_mix_0,P,X,T,d_air_0,n_air,n_mix,n_self);
    
        % Doppler and Collisional (Lorentzian) Full Width at Half Maximum Calculation
        [V_d,V_c] = Dop_and_Col_Terms(V,X,g_self,g_mix,P,T,M,X_mix);
    end
    
    % Define cuttoff value for spectral lines
    alpha_total = zeros(size(W_array,1),1);
    if P <= 2
        cutoff = 2; % [cm-1]
    else
        cutoff = 5; % [cm-1] greater than 2 atm
    end

    % Check if extra outputs are required
    extra_output = (nargout > 2);

    % Start a parallel pool if one is not already running
    pool = gcp('nocreate'); % Get current pool, don't create new one
    if isempty(pool)
        pool = parpool('local'); % Start a pool with default settings
    end

    % Generate Voigt Profile (line shape) of Each Transitions
    parfor i = 1:size(Transitions,1)
        % Find usable index relating to cutoff
        index = abs(W_array-V(i)) <= cutoff;

        % Skip line if no usable index in range was found
        if isempty(index)
            continue; % Skip to next line
        end

        % Create wavenumber array for line being processed.
        L_array = W_array(index);

        param = [V(i);V_d(i);V_c(i)];
        [Voigt_p] = VoigtProfileV3(L_array,param);

        % Absorbance Calculation
        alpha_L = S(i)*Voigt_p*P*X*L;
        alpha_temp = zeros(size(W_array,1),1);
        alpha_temp(index) = alpha_L;

        % Absorbance contribution of each transitions
        if extra_output
            alpha{i} = alpha_L;
            alpha_array{i} = L_array;
        end

        % Spectral feature with all transitions superimposed
        alpha_total = alpha_total + alpha_temp;
    end

    % Function to calculate transition parameters at specified conditions
    % Function that uses reference condition to calculate values at gas condition
    function [V,g_self,g_mix] = ShiftedParam(V_0,g_self_0,g_mix_0,P,X,T,d_air_0,n_air,n_mix,n_self)
        d_air = d_air_0 .* (296/T).^(n_air);
        V = V_0 + P * (1-X) .* d_air; % Transition linecenter shift
        g_self = g_self_0 .*(296/T).^(n_self);
        g_mix = g_mix_0 .* (296/T).^(n_mix);
    end
    
    % Function that take spectral parameters at gas condition to determine
    % doppler and collisional half width
    function [V_d,V_c] = Dop_and_Col_Terms(V,X,g_self,g_mix,P,T,M,X_mix)
        V_d = 7.1623e-7 .* V .* sqrt(T/M);
        for j = 1:size(g_self,1)
            if sum(g_mix(j,:) ~= 0, 2) == 1
                V_c(j) = 2 .* P .* (X.*g_self(j) + (1-X).*g_mix(j));
            else
                V_c(j) = 2 .* P .* (X.*g_self(j) + (g_mix(j,:)*X_mix));
            end
        end
    end

end