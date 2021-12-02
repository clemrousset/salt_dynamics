clear all; close all; format short

%
%==========================================================================
%
% Two layer model for ice salinity.
%
% I use a simplified Semtner model
% Two state variables: Tsu and h
%
% MV, Dec 2021 for Giulia
%
%
%==========================================================================
% 
% User setup
%
%==========================================================================
%

% parameterization choices
i_scheme = 1      % Salinity scheme; '1' - SI3 empirical scheme (exponential-based)

% numerical parameters
dt    = 86400.    % Time step
N_y   = 50.        % Number of years
Nl     = 5;       % Number of layers

h_min = 0.1       % Minimum thickness

Kelvin = 0        % 273.15; % 0 for Celsius; 273.15 for Kelvins, 

% tuning parameters
alb_i = 0.8       % Surface albedo
k = 2.            % Ice thermal conductivity
Q_w = 2.          % Oceanic heat flux

%==========================================================================

%--------------------------------------------------------------------------
%
% physical constants
%
rho       = 917.
L         = 335000.
epsilon   = 0.98;
sigma     = 5.67e-8;
eps_sigma = epsilon * sigma;

phi_new   = 0.75; % New ice liquid fraction

%
% bottom boundary condition
%
T_b        = -1.8 + Kelvin;
Sw         = 32.0096 % Sw value corresponding T_b with the VC19 liquidus

%
% calendar
%
Nts_1d    = 86400. / dt
Nts_1y    = 365 * Nts_1d;       % number of time steps in a year
N_ts      = 365. * Nts_1d * N_y % number of time steps

%--------------------------------------------------------------------------
% Regress sigma T4 into A+BT
% to get a linear function FLW = A + BT
%--------------------------------------------------------------------------

Ttest = [ -40.:1:0. ];
LW    = sigma*(Ttest+273.15).^4;
P     = polyfit(Ttest,LW,1);
A = P(2);
B = P(1);

%==========================================================================
%
% Model calculations
%
%==========================================================================

%--------------------------------------------------------------------------
%
% Heat and mass balance (surface temperature and ice thickness)
%
%--------------------------------------------------------------------------

h_i(1) = h_min         ;
T_su(1) = -5. + Kelvin ;

for i_time = 2:N_ts
    
    doy(i_time) = ( 1 + mod(i_time-1,Nts_1y) ) / Nts_1d ; % --- day of year
    
    day(i_time) = (i_time-1) / Nts_1d;                    % --- day of simulation

    % Temperature calculation
    Q_atm(i_time) = solar_flux(doy(i_time)) * (1. - alb_i ) + non_solar_flux(doy(i_time)); 
    k_h      = k / h_i(i_time-1);

    zp4     = 0 ;
    zp1     = - k_h - epsilon * B;
    zp0     =   Q_atm(i_time) + k_h * T_b - epsilon * ( A - B*Kelvin ) ;
    
    p       = [ zp4 0. 0. zp1 zp0 ]; 

    T_su(i_time)     = min( [ max( real(roots(p))) Kelvin ] );
    
    % Thickness calculation
    Q_net(i_time)    = Q_atm(i_time) - epsilon * ( A + B*(T_su(i_time) - Kelvin) );
    
    dh_dt(i_time)    = - 1/(rho*L) * ( Q_net(i_time) + Q_w );
    h_i(i_time)      = max([ h_i(i_time-1) + dh_dt(i_time) * dt h_min ]);
    
    % Diags
    Q_c(i_time)      = k_h * ( T_su(i_time) - T_b ); % diag
    dhbo_dt(i_time)  = - 1 / (rho*L) * ( Q_w + Q_c(i_time) );
    dhsu_dt(i_time)  = dh_dt(i_time) - dhbo_dt(i_time);
    
end

%--------------------------------------------------------------------------
%
% Temperature and brine salinity profile
%
%--------------------------------------------------------------------------

for i_time = 2:N_ts
    
   dT_dz = ( T_b - T_su(1,i_time) ) / h_i(1,i_time-1) ;
   dz    = h_i(1,i_time-1) / Nl ;
   zdT   = dT_dz * dz ;
   
   T_i_edge(1,i_time)       = T_su(1,i_time) ;
   T_i_edge(Nl+1,i_time) = T_b ;
   
   for k = 2:Nl
      T_i_edge(k,i_time) = T_i_edge(k-1,i_time) + zdT;
   end
   
   for k = 1:Nl
      T_i(k,i_time) = ( T_i_edge(k,i_time) + T_i_edge(k+1,i_time) ) / 2.;
   end
   
   zT  = T_i(:,i_time) - Kelvin;
   S_br(:,i_time) =  -18.7 * zT - 0.519 * zT.^2 - 0.00535 * zT.^3;
   
   zT  = T_i_edge(:,i_time) - Kelvin;
   S_br_edge(:,i_time) =  -18.7 * zT - 0.519 * zT.^2 - 0.00535 * zT.^3;
   
   
end

%--------------------------------------------------------------------------
%
% Salt dynamics
%
%--------------------------------------------------------------------------

S_i(1:Nl,1)  = phi_new * Sw        ;
S_bu(1,1) = mean( S_i(:,1) ) ;

%---------------------------------------------------------------------
%
if ( i_scheme == 1 ) % SI3 Empirical scheme --- based on exponentials
%
%---------------------------------------------------------------------

    for i_time = 2:N_ts

        swi_sm         = heaviside( T_su(1,i_time) * 2. ); % sado-maso switch for surface melt (no = 0, yes = 1)
        
        swi_bg         = heaviside( dhbo_dt(1,i_time) * 2. );
        
        % entrapment at sea ice base
        % S_new = ( S_old * h_old + phi_new * S_w * dh) / h_new
        if ( swi_bg > 0 )
            S_bu(1,i_time) = ( S_bu(1,i_time-1) .* h_i(1,i_time-1) + ...
                             swi_bg * phi_new * Sw * dhbo_dt(1,i_time) * dt ) ./ h_i(1,i_time); %% ;
        else
            S_bu(1,i_time) = S_bu(1,i_time-1);
        end
        
        % gravity drainage
        dS             = -max( S_bu(1,i_time) - 5., 0. ) * dt / 1.73e6 * ( 1. - swi_sm );
        S_bu(1,i_time) = S_bu(1,i_time) + dS;
        S_i(:,i_time)  = S_bu(1,i_time);
        
        % flushing
        dS             = -max( S_bu(1,i_time) - 2., 0. ) * dt / 8.64e5  * swi_sm;
        S_bu(1,i_time) = S_bu(1,i_time) + dS;
        S_i(:,i_time)  = S_bu(1,i_time);
        
    end

end

%---------------------------------------------------------------------
%
if ( i_scheme == 2 ) % Depth-dependent scheme
%
%---------------------------------------------------------------------
    
    %--------------------------
    %
    % preparing required fields
    %
    %--------------------------
    for i_time = 2:N_ts
        
        S_i(1:Nl,i_time) = S_i(1:Nl,i_time-1)        ;

        % gravity drainage
    
        % flushing
    
        %----------------------------------------
        %
        % entrapment during growth and remapping
        %
        %----------------------------------------
        zSh0 = zeros(Nl+1,1);
        dz0 = h_i(1,i_time-1) / Nl;
        zSh0(1:Nl,1) = S_i(:,i_time) * dz0;

        % basal salt entrapment
        zdhb = dhbo_dt(1,i_time) * dt;
        zSh0(Nl+1,1) = phi_new * Sw * zdhb;

        % old layers
        zl0 = zeros(Nl+2,1);
        dh0 = zeros(Nl+1,1);    
        zl0(1,1) = 0.;
        for k = 2:Nl+1
           zl0(k,1) = dz0 * (k-1);
        end
        zl0(Nl+2,1) = zl0(Nl+1,1) + zdhb;
        for k = 1:Nl+1
            dh0(k,1) = zl0(k+1,1) - zl0(k,1);
        end

        % new layers
        zl1 = zeros(Nl+1,1);
        dz1 = h_i(1,i_time) / Nl;
        dh1 = zeros(Nl,1);
        zl1(1,1) = 0.;
        for k = 2:Nl+1
           zl1(k,1) = dz1 * (k-1);
        end
        for k = 1:Nl
            dh1(k,1) = zl1(k+1,1) - zl1(k,1);
        end

        weights = zeros(Nl,Nl+1);

        for k1 = 1:Nl
            for k0=1:Nl+1
                zdenom = max ( [ dh0(k0,1) 1.0e-10 ] ); % ok
                zfac1  = min ( [ zl0(k0+1,1) zl1(k1+1,1) ] );
                zfac2  = max ( [ zl0(k0,1)   zl1(k1,1)   ] );
                znum   = zfac1 - zfac2;
                znum   = max ( [ 0. znum ] );
                weights(k1,k0) = znum / zdenom;
            end
        end

        zSh1 = zeros(Nl,1);
        for k1 = 1:Nl
           zSh1(k1,1) = sum(weights(k1,:).*zSh0(:,1)');
        end

        zarr = zSh1(:,1) ./ dh1(:,1)
        S_i(:,i_time) = zarr(:,1)';
        
        % diagnose bulk salinity
        S_bu(1,i_time) = mean( S_i(:,i_time) ) ;
    
    end
    
end
    
%
%---------------------------------------------------------------------


%==========================================================================
%
% Figures
%
%==========================================================================

           %-------------------
figure(1); % Model diagnostics
           %-------------------
        
subplot(2,3,1); hold on;
plot(day(2:N_ts),T_su(2:N_ts),'LineWidth', 2)
xlabel('day'); ylabel('T_{su} (°C)')
set(gca,'fontsize', 14); set(gca, 'FontName', 'Myriad Pro')

subplot(2,3,2); hold on
plot(day(2:N_ts),h_i(2:N_ts), 'LineWidth', 2)
xlabel('day'); ylabel('h_i (m)')
set(gca,'fontsize', 14); set(gca, 'FontName', 'Myriad Pro')

subplot(2,3,3); hold on;
plot(day(2:N_ts),Q_net(2:N_ts), 'k', 'LineWidth', 2)
plot(day(2:N_ts),solar_flux(doy(2:N_ts)), 'LineWidth', 2)
plot(day(2:N_ts),non_solar_flux(doy(2:N_ts)), 'LineWidth', 2)
plot(day(2:N_ts),-eps_sigma*(T_su(2:N_ts) - Kelvin + 273.15).^4,'r', 'LineWidth', 2)

legend('Q_{net}', 'Q_s', 'Q_{ns}', '-\epsilon\sigmaT_{su}^4' )
set(gca,'fontsize', 14); set(gca, 'FontName', 'Myriad Pro')
xlabel('day'); ylabel('Q (W/m^2)')

subplot(2,3,4); hold on;
plot(day(2:N_ts), dh_dt(2:N_ts)*86400. *100., 'LineWidth', 2 )
plot(day(2:N_ts), dhbo_dt(2:N_ts)*86400. *100., 'LineWidth', 2 )
plot(day(2:N_ts), dhsu_dt(2:N_ts)*86400. *100., 'LineWidth', 2 )
legend('tot', 'bot', 'surf')
set(gca,'fontsize', 14); set(gca, 'FontName', 'Myriad Pro')
xlabel('day'); ylabel('dh/dt (cm/day)')

subplot(2,3,5); hold on
plot(day(2:N_ts),S_bu(2:N_ts), 'LineWidth', 2)
xlabel('day'); ylabel('<S> (g/kg)')
set(gca,'fontsize', 14); set(gca, 'FontName', 'Myriad Pro')

%==========================================================================
%
% External functions
%
%==========================================================================

function [Fsw] = solar_flux(day)
  Fsw = 314. * exp(-(day-164).^2/4608.);
end

function [Fns] = non_solar_flux(day)
  Fns = 118. * exp(-0.5 * (day-206.).^2 / (53^2)) + 179.;
end

