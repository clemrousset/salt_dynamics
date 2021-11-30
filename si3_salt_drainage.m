%==================================================================================================
%                                        SI3_SALT_DRAINAGE
%                                        =================
%
% salinity drainage parameterization including simple ice growth model (with or without T_su change)
%
% we start with rn_himin ice thickness (typically 5cm) and make the ice grow such as:
% ==========================================================================
%   dh = - Fc * dt / ( rho_i * Lfus );
%
%      with   Fc  = ki/h * (T_su - T_b) conduction flux
%             ki  = heat conductivity in the ice (W/m/K)
%              h  = ice thickness
%          rho_i  = ice density
%           Lfus  = latent heat of fusion
%            Tsu  = surface temp
%            Tbo  = bottom temperature
%
% we do a vertical linear remapping on S (such as in SI3 for enthalpy)
% ======================================
%
% we then solve this equation:
% ===========================
%    dS/dt = -w dSbr/dz
%
%    with S   = sea ice salinity
%         Sbr = brine salinity
%         w   = upwelling Darcy velocity of the return flow (i.e. vertical velocity of the brines, positive downward => >0)
%
%    discrete form is solved using upward scheme (such as in CICE):
%    (S(t+dt)-S(t))/dt = -w(k) * (Sbr(k+1)-Sbr(k))/dz
%
% 2 schemes are proposed based on the paper from Thomas et al. (2020):
% ======================
%   0 |  ----------------------------------- surface
%     |                           
%     |  ----------------------------------- zc
%   z |
%     |      Ra > Rac => brine convection
%     |                          
%   h |  ------------------------------------ bottom
%     v
%
%     Ra = cl * g * beta_s * (Sbr(z) - Sw) * perm * (h-z) / (kl*visc)   [RWJ2014 formulation]
%
%        with Ra     : Rayleigh number
%             cl     : brine heat capacity         (J/m3/K)
%             g      : gravity                     (m/s2)
%             beta_s : saline density coefficient  (g/kg)-1
%             kl     : brine thermal conductivity  (W/m/K)
%             visc   : brine kinematic viscosity   (m2/s)
%             Sw     : ocean salinity              (g/kg)
%             zc     : critical depth below which convection occurs (m)
%             h      : total ice thickness         (m)
%             perm   : effective permeability      (m2)
%                      = 3.e-8    * (S/Sbr)^3      [SI3]     ( i_perm_for == 0 )
%                      = 1.995e-8 * (S/Sbr)^3.1    [Freitag] ( i_perm_for == 1 )
%                      = 1.0e-8   * (S/Sbr)^3.     [RJW2014] ( i_perm_for == 2 )
%
%    1) === Reese Jones & Worster 2014 (refer to as RJW2014) ===
%
%        w(z) = - alpha_rjw * kl / cl * max(Ra(z)-Rac) * (z-zc)/(h-zc)^2
%           with alpha_rjw ; intensity parameter
%
%    2) === Griewank & Notz 2013 (refer to as GN2013) ===
%
%        w(k) = - alpha_gn/rho * sum( (Ra(kk)-Rac) * dz(kk), [from kk=1 to k] )
%           with rho      : brine density       (kg/m3)
%                alpha_gn ; intensity parameter (kg/m3/s)
%
%    Other options:
%       liquidus formulation
%       S_br = - T / mu                                [linear liquidus]   ( i_liq == 1 )
%       S_br = -18.7 * T - 0.519 * T2 - 0.00535 * T3   [VC2019]            ( i_liq == 2 )
%       S_br = -17.6 * T - 0.389 * T2 - 0.00362 * T3   [Weast]             ( i_liq == 3 )
%
%==================================================================================================
% The 2 schemes can be compared to what is in the code right now     ( i_scheme = 3)
%                               to a scheme with prescribed velocity ( i_scheme = 0 )
%==================================================================================================
%
iprint=0
%
% --- calendar and layers --- %
Ndays = 100                 % Number of days
dt    = 10800               % Time step (3h)
Nts  = Ndays * 86400. / dt; % Total number of time steps (the first does not count, so 2 means 1 time step)
Nlays = 10;                 % Number of vertical layers (20 layers get unstable at 28800 time steps?)
Ni_S  = 20;                 % Number of iterations for salinity equation
                            %    needs 10  for GN13 with 3h time step and 5cm ice thickness (6 if 3 cm)
                            %    needs 11 for RJW14 with 3h time step and 5cm ice thickness (18 if 3cm)

% --- convection schemes --- %
i_scheme   = 1    % 0 = Prescribed velocity
                  % 1 = RJW2014
                  % 2 = GN2013
                  % 3 = SI3
i_explicit = 1    % 1 = explicit (upward)
                  % 0 = implicit (does not work for now, Mar 2019)

% --- liquidus relation --- %
i_liq  = 2;       % 1 = linear liquidus
                  % 2 = VC19 formulation
                  % 3 = Weast (used in RJW)

% --- rayleigh scheme --- %
i_Ra   = 1;       % 1 = RJW rayleigh formulation
                  % 2 = Vancoppenolle et al TCD2013
                  % for i_Ra = 2:
kappa = 1.2e-7;     % Thermal diffusivity of brine
mudyn = 1.9e-3;     % dynamical viscosity of brine
B_S   = 0.81;       % Sensitivity of density to salinity (kg/m3/(g/kg))              

% --- permeability --- %
i_perm_eff = 2;   % 1 = vertical minimum
                  % 2 = harmonic mean
i_perm_for = 2;   % 0 = SI3
                  % 1 = Freitag
                  % 2 = RJW

% --- temperature equation --- %
i_Tsu  = 0;       % 1 = computed temperature
                  % 0 = prescribed
                  % if i_Tsu == 1:
Flw        = 170.;
zeps_sigma = 0.95 * 5.67e-8; % Constant in Boltzman-law
Ni_T       = 5;              % Number of iterations for surface temperature

% --- Boundary conditions --- %
T_b      = -1.8 + Kelvin;  % bottom temperature
rn_himin = 0.05;           % mini ice thickness 

% --- Physical constants --- % 
Lfus        = 334000;             % Latent heat of freezing [J/kg]
rho_i       = 917;                % density of ice [kg/m^3]
k_i         = 2.14;               % heat conductivity of ice [W/(m K)] => clem: in SI3 it is around 2.0
Kelvin      = 273.15;             % 0C in Kelvins
mu          = 0.054;

% --- Prescribed velocity value (i_scheme = 0) --- %
w_pr = -2.0e-8;   % Prescribed brine velocity -2.0e-8 = guessed value

% --- RJW model physical constants (i_scheme = 1) --- %
cl     = 4.e6;    % Liquid volumetric heat capacity (J/m3)
g      = 9.81;    % Gravity (m/s^2)
beta_s = 7.5e-4;  % Beta_S (g/kg)-1
kl     = 0.523;   % W/m/K
visc   = 1.8e-6;  % Kinematic viscosity
                  % tuning constants (from Thomas et al. 2020)
Rc_RJW    = 2.9     % critical Rayleigh number
alpha_RJW = 0.13
                  % tuning constants (from Martin)
% $$$ Rc_RJW    = 2.7
% $$$ alpha_RJW = 0.037

% --- GN model constants (i_scheme = 2) --- %
rho_br_GN = 1020.    % Brine density (kg/m3)
                     % tuning constants (from Thomas et al. 2020)
Rc_GN     = 2.4        % Critical Rayleigh number
alpha_GN  = 6.7e-3     % Brine flow (kg/m3/s)
                     % tuning constants (from Martin)
% $$$ Rc_GN     = 1.01
% $$$ alpha_GN  = 1.56e-3

%--------------------------------------------------------------------------------------------------
% Vetal TCD 2013 Rayleigh constants

%==================================================================================================



if    ( i_scheme == 0 ); LeScheme = 'Scheme:CST';
elseif( i_scheme == 1 ); LeScheme = 'Scheme:RJW2014';
elseif( i_scheme == 2 ); LeScheme = 'Scheme:GN2013';
elseif( i_scheme == 3 ); LeScheme = 'Scheme:SI3'; end

if    ( i_liq == 1 ); LeLiq = 'Liq:linear';
elseif( i_liq == 2 ); LeLiq = 'Liq:VC19';
elseif( i_liq == 3 ); LeLiq = 'Liq:Weast'; end

if    ( i_perm_eff == 1 ); LePerm = 'Perm:mini';
elseif( i_perm_eff == 2 ); LePerm = 'Perm:harmo'; end

if    ( i_perm_for == 0 ); LePermF = 'PermF:SI3';
elseif( i_perm_for == 1 ); LePermF = 'PermF:Freitag';
elseif( i_perm_for == 2 ); LePermF = 'PermF:RJW2014'; end

LeTitre=[' ',LeScheme,' ',LeLiq,' ',LePerm,' ',LePermF,' '];

% Seawater salinity
if ( i_liq == 1 );    % linear liquidus
   Sw = -(T_b-Kelvin) / 0.054
elseif ( i_liq == 2); % Notz liquidus
   Sw = 32.0096
elseif ( i_liq == 3); % Weast liquidus
   Sw = 30.4408
end

% Arrays
h_i       = zeros(1,Nts);         % Ice thickness in m
T_su      = zeros(1,Nts);         % Surface temperature
w_br      = zeros(Nlays,Nts);        % Darcy velocity (positive downwards) (mid-points)
f_br      = zeros(Nlays,Nts);        % Brine fraction (mid-points)
% $$$ phi_edge  = zeros(Nlays+1,Nts);      % Brine fraction (edge of layers)
S_br      = zeros(Nlays,Nts);        % Brine salinity (mid-points)
% $$$ S_br_edge = zeros(Nlays+1,Nts);      % Brine salinity (edge-points)
S_i       = zeros(Nlays,Nts);        % Ice bulk salinity (mid-points)
T_i       = zeros(Nlays,Nts);        % Ice temperature (mid-moints)
T_i_edge  = zeros(Nlays,Nts);        % Ice temperature (edge-points)
z_mid     = zeros(Nlays,1);          % Edge-point normalized depth
z_edge    = zeros(Nlays+1,1);        % Edge-point temperature
days      = zeros(1,Nts);         % Days
Ra        = zeros(Nlays,Nts);        % Mid-point Rayleigh number
z_c       = zeros(1,Nts);         % Critical convective depth

delta_Sc  = zeros(1,Nts);         % Change in salt content (g/kg . m/s)
Fin       = zeros(1,Nts);         % Salt influx (g/kg.m/s)
Fout      = zeros(1,Nts);         % Salt outlufx (g/kg.m/s)
Rae       = zeros(1,Nts);         % Effective Ra
Ramax     = zeros(1,Nts);         % Mean Ra over convective layer

CFL       = zeros(1,Nts);         % CFL criterion
S_i_cvg1  = zeros(Nlays,(Nts-1)*Ni_S);
S_i_cvg2  = zeros(1,Nts);

%--------------------------------------------------------------------------------------------------

% Initial conditions
h_i (1,:) = rn_himin;
T_su(1,:) = -20. + Kelvin;
days(1,:) = [1:Nts]*dt/86400.;

% Initial temperatures
T_i_edge(1,1)       = T_su(1,1);
T_i_edge(Nlays+1,1) = T_b;
dT_dz = ( T_b - T_su(1,1) ) / h_i(1,1);
dz    = h_i(1,1) / Nlays;
zdT   = dT_dz * dz;
for k = 2:Nlays
   T_i_edge(k,1) = T_i_edge(k-1,1) + zdT;
end
for k = 1:Nlays
   T_i(k,1) = ( T_i_edge(k,1) + T_i_edge(k+1,1) ) / 2.;
end

% Initial salinities
zT  = T_i(:,1) - Kelvin
zT2 = zT .* zT 
zT3 = zT2 .* zT
if ( i_liq == 1 ); S_br(:,1) = - zT / mu; end
if ( i_liq == 2 ); S_br(:,1) = -18.7 * zT - 0.519 * zT2 - 0.00535 * zT3; end
if ( i_liq == 3 ); S_br(:,1) = -17.6 * zT - 0.389 * zT2 - 0.00362 * zT3; end

w_br(:,1) = 0.;
S_i (:,1) = Sw;
f_br(:,1) = S_i(:,1) ./ S_br(:,1);

%---------------------
% Vertical layers
%---------------------
z_edge(1,1) = 0;
for k = 2:Nlays+1
   z_edge(k,1) = (k-1)/Nlays;
end
for k = 1:Nlays
   z_mid(k,1) = (k-0.5)/Nlays;
end

%==================================================================================================
% Model loop
%==================================================================================================

for i_time = 2:Nts
   
   i_time
   
   Keff = k_i / h_i(1,i_time-1);
   %---------------------
   % Surface temperature
   %---------------------   
   if ( i_Tsu == 1 );
      zTsu = T_su(1,i_time-1);
      for i_iter = 1:Ni_T
         Fc = Keff * ( zTsu - T_b );
         Qnet       = Flw - zeps_sigma * zTsu^4- Fc;       % Net flux
         dQnet_dTsu = -4.*zeps_sigma * zTsu^3 - Keff;     % Net flux derivative
         delta_Tsu  = - Qnet / dQnet_dTsu;               % Tsu change
         zTsu       = zTsu + delta_Tsu;                  % Dummy Tsu
      end
      T_su(1,i_time) = min( [ zTsu, Kelvin ] );      
   end
   
   %---------------------
   % Ice Temperature
   %---------------------
   dT_dz = ( T_b - T_su(1,i_time) ) / h_i(1,i_time-1);
   dz    = h_i(1,i_time-1) / Nlays;
   zdT   = dT_dz * dz;
   
   T_i_edge(1,i_time)       = T_su(1,i_time);
   T_i_edge(Nlays+1,i_time) = T_b;
   
   for k = 2:Nlays
      T_i_edge(k,i_time) = T_i_edge(k-1,i_time) + zdT;
   end
   for k = 1:Nlays
      T_i(k,i_time) = ( T_i_edge(k,i_time) + T_i_edge(k+1,i_time) ) / 2.;
   end
   
   %---------------------
   % Brine salinity
   %---------------------
   
% $$$    % Edge of layers values
% $$$    zT  = T_i_edge(:,i_time) - Kelvin;    
% $$$    if ( i_liq == 1 ); % --- Linear liquidus
% $$$       S_br_edge(:,i_time) = - zT / mu; 
% $$$    elseif ( i_liq == 2 ); % --- VC19 liquidus
% $$$       S_br_edge(:,i_time) = -18.7 * zT - 0.519 * zT.^2 - 0.00535 * zT.^3; % ! gives high S value if freezing point at -1.8°C
% $$$    elseif ( i_liq == 3 ); % --- Weast 71 liquidus in RJW14
% $$$       S_br_edge(:,i_time) = -17.6 * zT - 0.389 * zT.^2 - 0.00362 * zT.^3;
% $$$    end
   
   % Mid-layer values
   zT  = T_i(:,i_time) - Kelvin;
   if ( i_liq == 1 ); % --- Linear liquidus
      S_br(:,i_time) = - zT / mu; 
   elseif ( i_liq == 2 ); % --- 3rd order liquidus, VC19
      S_br(:,i_time) =  -18.7 * zT - 0.519 * zT.^2 - 0.00535 * zT.^3; % ! gives high S value if freezing point at -1.8°C
   elseif ( i_liq == 3 ); % --- Weast 71 liquidus in RJW14
      S_br(:,i_time) =  -17.6 * zT - 0.389 * zT.^2 - 0.00362 * zT.^3;
   end
   
   %=======================
   % Salinity iteration
   %=======================
   
   S_i(:,i_time) = S_i(:,i_time-1); % iterated salinity
   for i_iter = 1:Ni_S
      
      %-------------------------------------------
      % 1) Brine volume fraction
      %-------------------------------------------
      
      % Brine volume fraction - mid points
      f_br(:,i_time) = S_i(:,i_time) ./ S_br(:,i_time);
      
      %-------------------------------------------
      % 2) Effective permeability
      %-------------------------------------------
      
      % Permeability (minimum)
      perm_eff      = zeros(Nlays,1);
      
      if ( i_perm_eff == 1 ); % Minimum
         
         for k1 = 1:Nlays; % --- mid points ---
            if ( i_perm_for == 0 )     % SI3
               perm_eff(k1,1) = 3.e-8 * (min(f_br(k1:Nlays,i_time)))^3;
            elseif ( i_perm_for == 1 )     % Freitag
               perm_eff(k1,1) = 1.995e-8 * (min(f_br(k1:Nlays,i_time)))^3.1;
            elseif ( i_perm_for == 2 ) % Rees Jones and Worster -> this case leads to bizzare results
               perm_eff(k1,1) = 1.0e-8 * (min(f_br(k1:Nlays,i_time)))^3.;
            end
         end
         
      elseif ( i_perm_eff == 2 ); % Harmonic Mean
         
         for k1 = 1:Nlays % --- mid points
            zsum = 0.;
            for k2 = k1:Nlays
               if ( i_perm_for == 0 );     % SI3
                  zperm = 3.e-8 * f_br(k2,i_time)^3;
               elseif ( i_perm_for == 1 );     % Freitag
                  zperm = 1.995e-8 * f_br(k2,i_time)^3.1;
               elseif ( i_perm_for == 2 ); % Rees Jones and Worster
                  zperm = 1.0e-8 * f_br(k2,i_time)^3.;
               end
               zsum = zsum + 1. / zperm;
            end
            zsum = 1. / zsum;
            perm_eff(k1,1) = ( Nlays - k1 + 1 ) * zsum;
         end
                  
      end
      
      %-------------------------------------------
      % 3) Rayleigh number
      %-------------------------------------------
      % Rayleigh number
      
      %--- Mid-points
      if ( i_Ra == 1 ); % formulation of Rees-Jones and Worster (JGR2014)
         
         z1 = cl * g * beta_s / ( kl * visc );
         z2 = S_br(:,i_time) - Sw;
         z3 = h_i(1,i_time-1) * ( 1. - z_mid(:,1) );
         z4 = perm_eff(:,1);
         
      elseif( i_Ra == 2 ); % formulation of Vancoppenolle et al (TCD2013)
         
         z1 = g / ( kappa * mudyn );
         z2 = B_S * ( S_br(:,i_time) - Sw );
         z3 = h_i(1,i_time-1) * ( 1. - z_mid(:,1) );
         z4 = perm_eff(:,1);
         
      end  
      Ra(:,i_time) = z1 .* z2 .* z3 .*z4;
      
      %---------------------
      % 4) Critical depth
      %---------------------
      % RJW
      % 1) if Ra is everywhere < Rc => no convection Rae =0
      % 2) else: convection until z_c
      % ---> if Ra(0) > Rc => full-depth convection (z_c = 0)
      % ---> if Ra(0) < Rc => find point where Ra >= Rc
      if (i_scheme == 1)        
         zRc_crit = Rc_RJW;
      elseif (i_scheme == 2)
         zRc_crit = Rc_GN;
      else
         zRc_crit = Rc_RJW;
      end

      %clem: rewrite because I do not understand how it is written
      i_c = Nlays+1;
      Rae(1,i_time) = 0;
      Ramax(1,i_time)=0;
      for k=Nlays:-1:1
         if ( Ra(k,i_time) >= zRc_crit );
            i_c = k;
            Rae  (1,i_time) = max( Rae  (1,i_time), Ra(k,i_time) - zRc_crit );
            Ramax(1,i_time) = max( Ramax(1,i_time), Ra(k,i_time) );

            z_c(1,i_time) = z_edge(k,1)*h_i(1,i_time-1);
         end
      end

      %-------------------------------------------
      % 5) Vertical velocity
      %-------------------------------------------
            
      if ( i_scheme == 0 );       %-------------------------
         w_br(:,i_time) = w_pr;   % 0 - Prescribed velocity
                                  %-------------------------
      elseif ( i_scheme == 1 );   % 1 - RJW 2014
                                  %-------------
         for k=1:Nlays
            if (k >= i_c && i_c<Nlays+1) %clem if ( k > i_c )
               zz = z_mid(k,1) * h_i(1,i_time-1);
               zc = z_c(1,i_time);
               zh = h_i(1,i_time-1);
               w_br(k,i_time) = - alpha_RJW * Rae(1,i_time) * ( kl / cl ) * ( zz - zc ) / ( zh - zc ).^2;
            else
               w_br(k,i_time) = 0.;
            end
         end
         %                         %------------
      elseif ( i_scheme == 2 );    % 2 - GN 2013
                                   %------------
         for k = 1:Nlays
            w_br(k,i_time) = - alpha_GN / rho_br_GN * sum ( ( Ra(1:k,i_time) - Rc_GN ) * dz );
         end
         
      end

      %---------------------
      % 7) Salinity
      %---------------------
      if ( i_explicit == 1 ); % Explicit Scheme

         if( i_scheme == 3 ) % SI3 scheme
            if( i_iter == 1 )
               for k = 1:Nlays
                  dS = -max( S_i(k,i_time) - 5., 0. ) * dt / 1.73e6;
                  S_i(k,i_time) = S_i(k,i_time) + dS;
               end
            end
            
         else
            Ak = zeros(Nlays,1);
            for k = 1:Nlays
               Ak(k) = w_br(k,i_time) * (dt/Ni_S) / dz;
            end
            
            % upstream scheme as in CICE: dS/dt = -w(k)/dz * ( Sbr(k+1) - Sbr(k) )
            for k = 1:Nlays-1
               dS = -Ak(k) * ( S_br(k+1,i_time) - S_br(k,i_time) );
               S_i(k,i_time) = S_i(k,i_time) + dS;
            end
            S_i(Nlays,i_time) = S_i(Nlays,i_time) - Ak(Nlays) * ( Sw - S_br(Nlays,i_time) );
            
            if( min(S_i(:,i_time)) <= 0. )
               'stop at i_explicit=1'
               return
            end
            
            CFL(1,i_time) = max( CFL(1,i_time), abs(max(Ak)) );
         end
% $$$             i_old = 0;
% $$$ 
% $$$             if ( i_old == 1 ) ; % Very first scheme I used
% $$$                for k = 1:Nlays
% $$$                   Xk = ( S_br_edge(k+1,i_time) - S_br_edge(k,i_time) ) / dz;
% $$$                   dS = - w_br(k,i_time) * Xk * dt;
% $$$                   %dS = max( [ -S_star(k,1) dS ] );
% $$$                   S_star(k,1) = S_i(k,i_time-1) + dS;
% $$$                end
% $$$ 
% $$$                zCFL_max = abs(max( w_br(:,i_time) * dt / dz ));
% $$$                CFL(1,i_time) = zCFL_max;
% $$$ 
% $$$             else % Scheme to be used in LIM1D with direct use of center-point values
% $$$                  % I havent figured why it is different from previous one.
% $$$                  % I have now: this is because of free or fixed surface brine salinity
% $$$                
% $$$                Ak = zeros(Nlays,1);
% $$$ 
% $$$                for k = 1:Nlays
% $$$                   Ak(k,1) = w_br(k,i_time) * dt / dz;
% $$$                   if ( k == 1 );
% $$$                      %... % S_br_1/2 is free in this version
% $$$                      %S_star(k,1) = S_br(k,i_time) * ( f_br(k,1) + Ak(k,1) ) +
% $$$                      %              S_br(k+1,i_time) * ( - Ak(k,1) ); %ok
% $$$                      %... % S_br_1/2 is fixed in this version
% $$$                      S_star(k,1) = S_br(k,i_time) * ( f_br(k,1) - Ak(k,1) / 2. ) + ...
% $$$                          S_br(k+1,i_time) * ( - Ak(k,1) / 2. ) + ...
% $$$                          S_br_edge(1,i_time) * Ak(k,1);
% $$$                   elseif ( k == Nlays );
% $$$                      S_star(k,1) = S_br(k-1,i_time) * ( Ak(k,1) / 2. ) + ...
% $$$                          S_br(k,i_time) * ( f_br(k,1) + Ak(k,1) / 2. ) + ...
% $$$                          Sw * ( - Ak(k,1) );
% $$$                   else;
% $$$                      S_star(k,1) = S_br(k-1,i_time) * Ak(k,1) / 2. + ...
% $$$                          S_br(k,i_time) * f_br(k,1) + ...
% $$$                          S_br(k+1,i_time) * ( - Ak(k,1) / 2. );
% $$$                   end;
% $$$                   % clem rewrite since I do not understand => no because we need w*gradS et pas div(w*S)
% $$$                   %Ak(k,1) = w_br(k,i_time) * dt / (Ni_S*dz);
% $$$                   %S_star(k,1) = S_star(k,1) + (dt/(Ni_S*dz)) * ( w_br(k+1,i_time) * S_br_edge(k+1,i_time) - w_br(k,i_time) * S_br_edge(k,i_time) ); 
% $$$ 
% $$$                end
% $$$                CFL(1,i_time) = max(CFL(1,i_time),max(abs(Ak)));
% $$$ 
% $$$             end % i_old
         
      elseif ( i_explicit == 0 ); % Implicit scheme --- This version does not work - dunno why ---
         
         Ak = zeros(Nlays,1);

         for k = 1:Nlays
            Ak(k,1) = w_br(k,i_time) * dt / dz;
         end
         
         ztrid = zeros (Nlays,Nlays);
         zind  = zeros (Nlays,1);
         
         ztrid(1,1) = f_br(1,i_time) - Ak(1,1);
         ztrid(1,2) = Ak(1,1);
         zind(1,1)  = S_i(1,i_time);
         
         for k = 2:Nlays-1
            ztrid(k,k-1) = - Ak(k,1) / 2.;
            ztrid(k,k) = f_br(k,i_time);
            ztrid(k,k+1) = Ak(k,1)/2.;
            zind(k,1) = S_i(k,i_time);
         end
         
         ztrid(Nlays,Nlays-1) = - Ak(Nlays,1)/2.;
         ztrid(Nlays,Nlays)   = f_br(Nlays,i_time) - Ak(Nlays,1)/2.;
         zind(Nlays,1)     = S_i(Nlays,i_time) - Ak(Nlays,1)*Sw;
         
         zSbrnew = zind' / ztrid;
         
         S_i(:,i_time) = zSbrnew' .* f_br(:,i_time);
         
         CFL(1,i_time) = abs(max(Ak));
         
      end

      % record convergence
      S_i_cvg1(:,(i_time-2)*Ni_S+i_iter) = S_i(:,i_time);
      
   end % end of iterative procedure

   % record convergence
   S_i_cvg2(:,i_time) = max(S_i_cvg1(:,(i_time-1)*Ni_S-1) - S_i_cvg1(:,(i_time-1)*Ni_S));
            
   %-------------------------------
   % Diagnose conservation of salt
   %-------------------------------
   %delta_Sc  = zeros(1,Nts);         % Change in salt content (g/kg . m/s)
   %Fin       = zeros(1,Nts);         % Salt influx (g/kg.m/s)
   %Fout      = zeros(1,Nts);         % Salt outlufx (g/kg.m/s)
   
   delta_Sc(1,i_time)   = ( sum(S_i(:,i_time))  - ...
                            sum(S_i(:,i_time-1)) ) ...
                            * dz / dt;
   
   Fin(1,i_time)  = -w_br(Nlays,i_time) * Sw;
   Fout(1,i_time) = delta_Sc(1,i_time) - Fin(1,i_time);
   
   %---------------------
   % Thickness change
   %---------------------
   Fc   = Keff * ( T_su(1,i_time) - T_b );
   zdhb = - Fc * dt / ( rho_i * Lfus );
   h_i(1,i_time) = h_i(1,i_time-1) + zdhb;
   
   dz0 = h_i(1,i_time-1) / Nlays;
   dz1 = h_i(1,i_time) / Nlays;

   if( i_scheme == 3)
      % New bottom ice salinity (Cox & Weeks, JGR88 )
      zgrr = min( 1.e-3, max ( zdhb / dt , 1.e-10 ) );
      if( zgrr < 2.e-8 )
         zfracs = 0.12;
      elseif( zgrr > 3.6e-7 )
         zfracs = min( 0.5, 0.26 / ( 0.26 + 0.74 * exp( - 724300.0 * zgrr ) ) );
      else
         zfracs = 0.8925 + 0.0568 * log( 100.0 * zgrr );
      end
   else
      zfracs = 0.75; % just a guess of the new salinity associated with basal growth   
   end
   S_i_new = zfracs * Sw;
   
   %----------------
   % Remapping
   %----------------
   zh_cum0 = zeros(Nlays+2,1);
   zh_cum1 = zeros(Nlays+1,1);
   zsh_cum0 = zeros(Nlays+2,1);
   zsh_cum1 = zeros(Nlays+1,1);
   % first level = surface interface = 0
   zsh_cum0(1) = 0; zh_cum0 (1) = 0;
   % Nlays levels cumulative sum 
   zsh_cum0(2:Nlays+1) = cumsum( S_i(:,i_time) * dz0 );
   zh_cum0 (2:Nlays+1) = dz0 * [1:Nlays];
   % Nlays+2 level = ice-ocean interface = new thickness and salinity
   zsh_cum0(Nlays+2) = zsh_cum0(Nlays+1) + S_i_new * zdhb;
   zh_cum0 (Nlays+2) = zh_cum0 (Nlays+1) + zdhb;
   
   % cumulative sum of the new salt content
   zsh_cum1(1) = 0; zh_cum1 (1) = 0;
   zh_cum1 (2:Nlays+1) = dz1 * [1:Nlays];
   for k0 = 2:Nlays+2
      for k1 = 2:Nlays
         if( zh_cum1(k1) <= zh_cum0(k0) && zh_cum1(k1) > zh_cum0(k0-1) )
            zsh_cum1(k1) = ( zsh_cum0(k0-1) * ( zh_cum0(k0) - zh_cum1(k1  ) ) +  ...
                             zsh_cum0(k0  ) * ( zh_cum1(k1) - zh_cum0(k0-1) ) )  ...
                             / ( zh_cum0(k0) - zh_cum0(k0-1) );
         end
      end
   end
   zsh_cum1(Nlays+1) = zsh_cum0(Nlays+2); % ensure strict conservation

   % remapped salinity
   for k1 = 1:Nlays
      S_i(k1,i_time) = ( zsh_cum1(k1+1) -  zsh_cum1(k1) ) / dz1;
   end

   % SI3 scheme
   if( i_scheme == 3 )
      zs = mean( S_i(:,i_time) );
      S_i(:,i_time) = zs;
   end
end

%==================================================================================================
% PLOTS
%==================================================================================================
if( i_scheme == 3 )
   zdh_dt = [1.:0.5:1000]*1.e-8;
   % New bottom ice salinity (Cox & Weeks, JGR88 )
   zgrr = zdh_dt;
   if( zgrr < 2.e-8 )
      zfracs = 0.12; 
   elseif( zgrr > 3.6e-7 )
      zfracs = min( 0.5, 0.26 / ( 0.26 + 0.74 * exp( - 724300.0 * zgrr ) ) );
   else
      zfracs = 0.8925 + 0.0568 * log( 100.0 * zgrr );
   end
   initfig('land')
   plot(zdh_dt, zfracs*Sw); grid on; xlabel('dh/dt'); ylabel('Snew');
   if iprint==1
      printfig('/usr/home/crlod/WORK/NEMO/PROG/fig_local/SALT/',...
               ['salt_newice_scheme',num2str(i_scheme)])
   end
end

%--------------
% convergence
%--------------
initfig('land')

subplot(1,2,1)
plot(days,CFL); xlabel('days'); ylabel('CFL'); grid on
subplot(1,2,2)
plot(days,S_i_cvg2); xlabel('days'); ylabel('salinity convergence (g/kg)');grid on

% title
ax = [0.3 .92 0.45 0.02];
axes('position',ax);
title(['Iterations=',num2str(Ni_S),' (',LeTitre,')'])
axis off;

% prints
if iprint==1
   printfig('/usr/home/crlod/WORK/NEMO/PROG/fig_local/SALT/',...
            ['salt_cvg_scheme',num2str(i_scheme),'_liq',num2str(i_liq),'_permeff',num2str(i_perm_eff),'_permfor',num2str(i_perm_for),'_lay',num2str(Nlays)])
end


%--------------
% Time series
%--------------
initfig('land')

if ( i_Tsu == 0 );
   h_stefan = sqrt (- 2*k_i/(rho_i * Lfus) * (T_su-T_b) .* days * 86400.);
end
subplot(2,3,1); grid on
plot(days,T_su-Kelvin); ylabel('T_{su}(C)')

subplot(2,3,2); hold on; grid on
plot(days,h_i); xlabel('days'); ylabel('h_i (m)')
if ( i_Tsu == 0 );
   plot(days,h_stefan, 'k.')
end
plot(days, z_c, 'r--')

subplot(2,3,3); hold on; grid on
plot(days,T_i_edge-Kelvin); xlabel('days'); ylabel('T_i(C)');
plot(days,T_i-Kelvin, '--');

subplot(2,3,4); hold on; grid on
plot(days,S_br); xlabel('days'); ylabel('S_{br}(g/kg)');

subplot(2,3,5); hold on; grid on
plot(days,S_i); xlabel('days'); ylabel('S_{i}(g/kg)');
plot(days,mean(S_i,1), 'k', 'LineWidth', 1);

subplot(2,3,6); hold on; grid on
plot(days,w_br); xlabel('days'); ylabel('w_{br}(m/s)');

% title
ax = [0.3 .92 0.45 0.02];
axes('position',ax);
title([LeTitre])
axis off;

% prints
if iprint==1
   printfig('/usr/home/crlod/WORK/NEMO/PROG/fig_local/SALT/',...
            ['salt_ts_scheme',num2str(i_scheme),'_liq',num2str(i_liq),'_permeff',num2str(i_perm_eff),'_permfor',num2str(i_perm_for),'_lay',num2str(Nlays)])
end


%----------
% Profiles
%----------
initfig('land') 
subplot(2,3,1) %--- S ---
hold on; grid on
xlabel('S_i (g/kg)'); ylabel('z (m)');

for i_time = 1:Nts
   if ( days(i_time) == 0.5 | ...
        days(i_time) == 1.  | ...
        days(i_time) == 2.  | ...
        days(i_time) == 5.  | ...
        days(i_time) == 10. | ...
        days(i_time) == 20. )
      plot(S_i(:,i_time),-z_mid(:,1)*h_i(1,i_time), 'LineWidth', 1); 
   end
end
xlim([0. 31.]); ylim([-0.7 0]);


subplot(2,3,2) %--- T ---
hold on; grid on
xlabel('T_i'); ylabel('z (m)');
xlim([-20 0]); ylim([-0.7 0]);
for i_time = 1:Nts
   if ( days(i_time) == 0.5 | ...
        days(i_time) == 1.  | ...
        days(i_time) == 2.  | ...
        days(i_time) == 5.  | ...
        days(i_time) == 10. | ...
        days(i_time) == 20. )
      plot(T_i(:,i_time)-Kelvin,-z_mid(:,1)*h_i(1,i_time), 'LineWidth', 1); 
   end
end

subplot(2,3,3) %--- Ra ---
hold on; grid on
xlabel('Ra'); ylabel('z (m)')
for i_time = 1:Nts
   if ( days(i_time) == 0.5 | ...
        days(i_time) == 1.  | ...
        days(i_time) == 2.  | ...
        days(i_time) == 5.  | ...d
      days(i_time) == 10. | ...
          days(i_time) == 20. )
      plot(Ra(:,i_time),-z_mid(:,1)*h_i(1,i_time), 'LineWidth', 1); 
   end
end

subplot(2,3,4) %--- w ---
hold on; grid on
xlabel('w (m/s)'); ylabel('z (m)')
for i_time = 1:Nts
   if ( days(i_time) == 0.5 | ...
        days(i_time) == 1.  | ...
        days(i_time) == 2.  | ...
        days(i_time) == 5.  | ...
        days(i_time) == 10. | ...
        days(i_time) == 20. )
      plot(w_br(:,i_time),-z_mid(:,1)*h_i(1,i_time), 'LineWidth', 1); 
   end
end

subplot(2,3,5) %--- D ---
hold on; grid on
xlabel('\phi'); ylabel('z (m)')
for i_time = 1:Nts
   if ( days(i_time) == 0.5 | ...
        days(i_time) == 1.  | ...
        days(i_time) == 2.  | ...
        days(i_time) == 5.  | ...
        days(i_time) == 10. | ...
        days(i_time) == 20. )
      plot(f_br(:,i_time),-z_mid(:,1)*h_i(1,i_time), 'LineWidth', 1); 
   end
end
legend('0.5d', '1d', '2d', '5d', '10d', '20d','location','best' )
legend('boxoff')


% title
ax = [0.3 .92 0.45 0.02];
axes('position',ax);
title([LeTitre])
axis off;

% prints
if iprint==1
   printfig('/usr/home/crlod/WORK/NEMO/PROG/fig_local/SALT/',...
            ['salt_profiles_scheme',num2str(i_scheme),'_liq',num2str(i_liq),'_permeff',num2str(i_perm_eff),'_permfor',num2str(i_perm_for),'_lay',num2str(Nlays)])
end


%
% $$$ initfig('portrait')
% $$$ 
% $$$ subplot(2,4,1)
% $$$ plot(Fin, 'linewidth', 1); ylabel('F_{in} (g/kg.m/s)')
% $$$ set(gca,'fontsize', 12, 'FontName', 'Myriad Pro')
% $$$ grid on
% $$$ 
% $$$ subplot(2,4,2)
% $$$ plot(Fout, 'linewidth', 1); ylabel('F_{out} (g/kg.m/s)');
% $$$ set(gca,'fontsize', 12, 'FontName', 'Myriad Pro')
% $$$ grid on
% $$$ 
% $$$ subplot(2,4,3)
% $$$ plot(delta_Sc, 'linewidth', 1); ylabel('\Delta_{Sc} (g/kg.m/s)')
% $$$ set(gca,'fontsize', 12, 'FontName', 'Myriad Pro')
% $$$ grid on
% $$$ 
% $$$ subplot(2,4,4)
% $$$ plot(Fout+Fin-delta_Sc, 'linewidth', 1), ylabel('Conservation test (g/kg.m/s)')
% $$$ set(gca,'fontsize', 12, 'FontName', 'Myriad Pro')
% $$$ grid on
%%% Fc vs Ramax
% $$$ subplot(2,4,5); box on; hold on
% $$$ zaddr = find(z_c == 0); %full-depth mode
% $$$ zccc = [ 1 0 0];
% $$$ plot(Ramax(1,zaddr),-(Fin(1,zaddr)+Fout(1,zaddr)), 'ko', 'MarkerFaceColor', zccc); 
% $$$ zaddr = find(z_c > 0); %
% $$$ zccc = [ 0 0 1];
% $$$ plot(Ramax(1,zaddr),-(Fin(1,zaddr)+Fout(1,zaddr)), 'ko', 'MarkerFaceColor', zccc); 
% $$$ legend('full-depth', 'bottom')
% $$$ 
% $$$ xlabel('Ra_{max}'); ylabel('-(F_{in}+F_{out}) (g/kg.m/s)');
% $$$ set(gca,'fontsize', 12, 'FontName', 'Myriad Pro')
% $$$ xlim([0. 15.]) ;
% $$$ ylim([0. 1.0e-4]);

%%% Fc vs 1/hi
% $$$ subplot(2,4,6); hold on; box on
% $$$ zaddr = find(z_c == 0);
% $$$ zccc = [ 1 0 0];
% $$$ plot(1./h_i(1,zaddr),-(Fin(1,zaddr)+Fout(1,zaddr)), 'ko', 'MarkerFaceColor', zccc); 
% $$$ zaddr = find(z_c > 0);
% $$$ zccc = [ 0 0 1];
% $$$ plot(1./h_i(1,zaddr),-(Fin(1,zaddr)+Fout(1,zaddr)), 'ko', 'MarkerFaceColor', zccc); 
% $$$ 
% $$$ legend('full-depth', 'bottom')
% $$$ xlabel('h_i^{-1} (m^{-1})'); ylabel('-(F_{in}+F_{out}) (g/kg.m/s)');
% $$$ set(gca,'fontsize', 12, 'FontName', 'Myriad Pro')
% $$$ xlim([0. 15.]) ;
% $$$ ylim([0. 1.0e-4]);

%plot(Rae); plot(Ramax); ylabel('Ra_{e} / Ra_{max}'); xlabel('time step')
%set(gca,'fontsize', 12, 'FontName', 'Myriad Pro')

% $$$ subplot(2,4,7)
% $$$ zzz = diff(h_i);
% $$$ plot(zzz,delta_Sc(1:Nts-1),'ko'); xlabel('dh/dt'); ylabel('\Delta Sc (g/kg.m/s')
% $$$ set(gca,'fontsize', 12, 'FontName', 'Myriad Pro')
% $$$ grid on
% $$$ 
% $$$ subplot(2,4,8)
% $$$ plot(f_br,-w_br,'o')
% $$$ xlabel('\phi'), ylabel('w_{br}')
% $$$ xlim([0. 1.]); ylim([0. 1.0e-6]);
% $$$ if ( i_scheme == 1 ); title('RJW14'); end
% $$$ if ( i_scheme == 2 ); title('GN13'); end
% $$$ 
% $$$ set(gca,'fontsize', 12, 'FontName', 'Myriad Pro')
% $$$ grid on
% $$$ 
% $$$ if ( i_scheme == 1 ); % RJW SPECIFIC PLOT TO UNDERSTAND DYNAMICS
% $$$ 
% $$$    initfig('portrait'); title('Two modes of drainage for RJW')
% $$$    subplot(1,3,1); hold on; box on
% $$$    plot(h_i-z_c, 'LineWidth', 1); ylabel('h_c (m) & h (m)')
% $$$    plot(h_i, 'LineWidth', 2); 
% $$$    set(gca,'fontsize', 12, 'FontName', 'Myriad Pro')
% $$$ 
% $$$    subplot(1,3,2); hold on; box on
% $$$    for i = 1350:1360
% $$$       plot(Ra(:,i),-z_mid*h_i(1,i), 'LineWidth', 1)
% $$$       plot( zRc_crit , -z_c(1,i) , 'ko', 'MarkerFaceColor', 'k' )
% $$$    end
% $$$    plot( [zRc_crit zRc_crit ], [-1. 0], ':' )
% $$$    xlabel('Ra'), ylabel('z (m)')
% $$$ 
% $$$    subplot(1,3,3); hold on; box on
% $$$    for i = 1350:1360
% $$$       plot(w_br(:,i),-z_mid*h_i(1,i), 'LineWidth', 1)
% $$$    end
% $$$    xlabel('w (m/s)'), ylabel('z (m)')
% $$$ 
% $$$ end
