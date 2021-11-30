%==================================================================================================
%                                        SI3_SALT_FLUSHING
%                                        =================
%
% salinity flushing parameterization including simple ice melt model (with or without T_su change)
%    based on M. Vancoppenolle et al. 2007
% 
% we start with rn_himax ice thickness (typically 1m) and make the ice melt such as:
% ========
%   dh = - (1-Ftr) * Qsolar * dt / ( rho_i * Lfus );
%
%      with  Qsolar  = surface heat flux
%               Ftr  = part of solar flux that penetrates in the ice
%             rho_i  = ice density
%              Lfus  = latent heat of fusion
%
% we compute new ice temperature giving a fraction of solar flux penetrating thru the ice
% ==============================
%   Qabs = Ftr * ( - exp( -kappa * zh(k+1) ) + exp( -kappa * zh(k) ) )
%
%      with kappa = extinction coefficient
%              zh = depth 
%
%   dT is calcuated from Qabs = rho * cp * h * dT/dt
%
% we compute the new brine salinity
% =================================
%
% we then solve this equation:
% ===========================
%    dS/dt = -w dSbr/dz
%
%    with Sbr = brine salinity
%         w   = upwelling velocity (i.e. vertical velocity of the brines, negative upward => < 0)
%
%         w   = Fmass / rhob            if fr_b > fr_bc (= 5%)
%             = 0                       otherwise
%
%           with Fmass = -Flush * rhoi * dh / dt : mass flux (kg/m2/s, >0 since dh<0)
%                rhob = rhow * ( 1 + c*Sbr )     : brine density
%                fr_b = S / Sbr                  : brine volume fraction
%                rhoi                            : ice density
%                rhow                            : fresh water density (kg/m3)
%                c                               : empirical coef (0.8e-3 ‰-1)
%           tuning parameters:
%                Flush                           : fraction of melt water allowed to percolate thru the ice (30%)
%                fr_bc                           : critical brine volume above which there is flushing (5%)
%
%    discrete form is solved using upward scheme (such as in CICE):
%    (S(t+dt)-S(t))/dt = -w(k) * (Sbr(k-1)-Sbr(k))/dz
%
%    Other options:
%       liquidus formulation
%       S_br = - T / mu                                [linear liquidus]   ( i_liq == 1 )
%       S_br = -18.7 * T - 0.519 * T2 - 0.00535 * T3   [VC2019]            ( i_liq == 2 )
%       S_br = -17.6 * T - 0.389 * T2 - 0.00362 * T3   [Weast]             ( i_liq == 3 )
%
%==================================================================================================
%
iprint=0
%
% Calendar and layers
%
Ndays = 100                  % Number of days
dt    = 10800                % Time step (3h)
Nts   = Ndays * 86400. / dt; % Total number of time steps (the first does not count, so 2 means 1 time step)
Nlays = 10;                  % Number of vertical layers (20 layers get unstable at 28800 time steps?)
Ni_S  = 1;                   % Number of iterations for salinity equation
                             %    needs just 1 with 3h time step and 5cm minimum ice thickness
%==================================================================================================

i_melt = 1 % 0=fake melting, 0=real melting
i_Sini = 1 % 0=homogeneous, 1=C-shape

% --- flushing schemes --- %
i_scheme   = 1             % 1 = VC2007
                           % 2 = SI3

% --- liquidus relation --- %
i_liq  = 2;                 % 1 = linear liquidus
                            % 2 = VC19 formulation
                            % 3 = Weast (used in RJW)

% --- Boundary conditions --- %
Qsolar   = 50.;               % total solar flux reaching the ice (W/m2)
Ftr      = 0.3;               % part of solar flux transmitted thru surface scattering layer
Kelvin   = 273.15;            % 0C in Kelvins
T_b      = -1.8 + Kelvin;
rn_himin = 0.05;              % mini ice thickness 

% --- Physical constants --- % 
hi_ssl  = 0.0;               % surface scattering layer (10cm for sea ice) => 0 otherwise inversion temp.
kappa_i = 1.;                % exctinction coefficient in sea ice (m-1)
rcpi    = 2067;              % specific heat of fresh ice (J/kg/K)
Lfus    = 334000;            % Latent heat of freezing [J/kg]
rho_i   = 917;               % density of ice [kg/m^3]
rho_w   = 1000;              % density of water [kg/m^3]
rho_b   = 1000;              % density of the brines
                             %    in theory = rhow*(1+c*Sbr), with c = 0.8e-3
mu      = 0.054;

% --- tuning parameters --- % 
f_flush    = 0.3;
e_treshold = 0.05;
%==================================================================================================

if    ( i_scheme == 1 ); LeScheme = 'Scheme:VC2007';
elseif( i_scheme == 2 ); LeScheme = 'Scheme:SI3'; end

if    ( i_liq == 1 ); LeLiq = 'Liq:linear';
elseif( i_liq == 2 ); LeLiq = 'Liq:VC19';
elseif( i_liq == 3 ); LeLiq = 'Liq:Weast'; end

LeTitre=[' ',LeScheme,' ',LeLiq,' '];

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
S_br      = zeros(Nlays,Nts);        % Brine salinity (mid-points)
S_i       = zeros(Nlays,Nts);        % Ice bulk salinity (mid-points)
T_i       = zeros(Nlays,Nts);        % Ice temperature (mid-moints)
T_i_edge  = zeros(Nlays,Nts);        % Ice temperature (edge-points)
z_mid     = zeros(Nlays,1);          % Edge-point normalized depth
z_edge    = zeros(Nlays+1,1);        % Edge-point temperature
days      = zeros(1,Nts);         % Days


CFL       = zeros(1,Nts);         % CFL criterion
S_i_cvg1  = zeros(Nlays,(Nts-1)*Ni_S);
S_i_cvg2  = zeros(1,Nts);

%--------------------------------------------------------------------------------------------------

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

% Initial conditions
h_i (:)   = 1.;
T_su(:)   = 0. + Kelvin;
days(:)   = [1:Nts]*dt/86400.;
if( i_Sini == 0 )
   S_i (:,1) = 7;
elseif( i_Sini == 1 )
   S_i (:,1) = 10*(1.4-sqrt(1-(2*z_mid(:)-1).^2)); % C shape init
end
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

% Initial brine salinities
zT  = T_i(:,1) - Kelvin
zT2 = zT .* zT 
zT3 = zT2 .* zT
if ( i_liq == 1 ); S_br(:,1) = - zT / mu; end
if ( i_liq == 2 ); S_br(:,1) = -18.7 * zT - 0.519 * zT2 - 0.00535 * zT3; end
if ( i_liq == 3 ); S_br(:,1) = -17.6 * zT - 0.389 * zT2 - 0.00362 * zT3; end

f_br(:,1) = S_i(:,1) ./ S_br(:,1);

%==================================================================================================
% Model loop
%==================================================================================================

for i_time = 2:Nts
   
   i_time

   %---------------------
   % Thickness change
   %---------------------
   zdhtot = - (1.-Ftr) * Qsolar * dt / ( rho_i * Lfus );

   if( i_melt == 0 )
      h_i(1,i_time) = h_i(1,i_time-1);
      zdh = 0;
   else
      h_i(1,i_time) = max( rn_himin, h_i(1,i_time-1) + zdhtot );
      zdh = h_i(1,i_time) - h_i(1,i_time-1);
   end
   
   dz0 = h_i(1,i_time-1) / Nlays;
   dz1 = h_i(1,i_time)   / Nlays;
   
   %----------------
   % Remapping
   %----------------
   zh_cum0  = zeros(Nlays+1,1); zh_cum1  = zeros(Nlays+1,1);
   zsh_cum0 = zeros(Nlays+1,1); zsh_cum1 = zeros(Nlays+1,1);
   % Nlays levels cumulative sum
   zhmelt = zdh;
   for k0 = 2:Nlays+1
      zh_cum0 (k0) = zh_cum0 (k0-1) + max( 0., dz0 + zhmelt);
      zsh_cum0(k0) = zsh_cum0(k0-1) + S_i(k0-1,i_time-1) * max( 0., dz0 + zhmelt);
      zhmelt = min( 0., zhmelt + dz0 );
   end
   zh_cum1(2:Nlays+1) = [1:Nlays] * dz1;
   % cumulative sum of the new salt content
   zsh_cum1(1) = 0; zh_cum1 (1) = 0;
   for k0 = 2:Nlays+1
      for k1 = 2:Nlays
         if( zh_cum1(k1) <= zh_cum0(k0) && zh_cum1(k1) > zh_cum0(k0-1) )
            zsh_cum1(k1) = ( zsh_cum0(k0-1) * ( zh_cum0(k0) - zh_cum1(k1  ) ) +  ...
                             zsh_cum0(k0  ) * ( zh_cum1(k1) - zh_cum0(k0-1) ) )  ...
                             / ( zh_cum0(k0) - zh_cum0(k0-1) );
         end
      end
   end
   zsh_cum1(Nlays+1) = zsh_cum0(Nlays+1); % ensure strict conservation

   % remapped salinity
   for k1 = 1:Nlays
      Snew = ( zsh_cum1(k1+1) -  zsh_cum1(k1) ) / dz1;
      %if( Snew > S_i(k1,i_time) )
      S_i(k1,i_time) = Snew;
   end

   % SI3 scheme
   if( i_scheme == 2 )
      zs = mean( S_i(:,i_time) );
      S_i(:,i_time) = zs;
   end

   if( min(S_i(:,i_time)) <= 0. )
      'stop s_i<0'
      return
   end

   %---------------------
   % Ice Temperature
   %---------------------
   for k = 1:Nlays
      Qabs = Ftr * Qsolar * ( - exp( -kappa_i * max( 0., zh_cum1(k+1) - hi_ssl ) ) ...
                              + exp( -kappa_i * max( 0., zh_cum1(k)   - hi_ssl ) ) );
      dT = Qabs * dt / (rho_i*rcpi*dz1);
      T_i(k,i_time) = min( -0.01 + Kelvin, T_i(k,i_time-1) + dT );
   end
   
   %---------------------
   % Brine salinity
   %---------------------
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
   for i_iter = 1:Ni_S
      
      %-------------------------------------------
      % 1) Brine volume fraction
      %-------------------------------------------
      f_br(:,i_time) = S_i(:,i_time) ./ S_br(:,i_time);
      

      %---------------------
      % 7) Salinity
      %---------------------
      if( i_scheme == 2 ) % SI3 scheme
         if( i_iter == 1 )
            for k = 1:Nlays
               dS = -max( S_i(k,i_time) - 2., 0. ) * dt / 8.64e5;
               S_i(k,i_time) = S_i(k,i_time) + dS;
            end
         end
         
      else

         if( min( f_br(:,i_time) ) >= e_treshold )

            %-------------------------------------------
            % 5) Vertical velocity
            %-------------------------------------------
            w_br(:,i_time) = -f_flush * (zdhtot/dt) * (rho_i/rho_b); %>0         

            Ak = zeros(Nlays,1);
            for k = 1:Nlays
               Ak(k) = w_br(k,i_time) * (dt/Ni_S) / dz1;
            end
            
            % upstream scheme as in CICE: dS/dt = -w(k)/dz * ( Sbr(k) - Sbr(k-1) )
            S_i(1,i_time) = S_i(1,i_time) - Ak(1) * ( S_br(1,i_time) - 0. );            
            for k = 2:Nlays
               dS = -Ak(k) * ( S_br(k,i_time) - S_br(k-1,i_time) );
               S_i(k,i_time) = S_i(k,i_time) + dS;
               if( dS > 0 )
                  'stop d(s_i)>0'
                  return
               end
            end
            
            if( min(S_i(:,i_time)) <= 0. )
               'stop s_i<0'
               return
            end
            
            CFL(1,i_time) = max( CFL(1,i_time), abs(max(Ak)) );
         end

      end

      % record convergence
      S_i_cvg1(:,(i_time-2)*Ni_S+i_iter) = S_i(:,i_time);
      
   end % end of iterative procedure

   % record convergence
   if( Ni_S > 1 )
      S_i_cvg2(:,i_time) = max(S_i_cvg1(:,(i_time-1)*Ni_S-1) - S_i_cvg1(:,(i_time-1)*Ni_S));
   end
end

%==================================================================================================
% PLOTS
%==================================================================================================
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
            ['flushing_cvg_scheme',num2str(i_scheme),'_liq',num2str(i_liq),'_lay',num2str(Nlays)])
end


%--------------
% Time series
%--------------
initfig('land')

subplot(2,3,1); grid on
plot(days,T_su-Kelvin); ylabel('T_{su}(C)')

subplot(2,3,2); hold on; grid on
plot(days,h_i); xlabel('days'); ylabel('h_i (m)')

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
            ['flushing_ts_scheme',num2str(i_scheme),'_liq',num2str(i_liq),'_lay',num2str(Nlays)])
end


%----------
% Profiles
%----------
initfig('land') 
subplot(2,2,1) %--- S ---
hold on; grid on
xlabel('S_i (g/kg)'); ylabel('z frac');

for i_time = 1:Nts
   if ( days(i_time) == 1.   | ...
        days(i_time) == 20.  | ...
        days(i_time) == 40.  | ...
        days(i_time) == 60.  | ...
        days(i_time) == 80.  | ...
        days(i_time) == 100. )
      plot(S_i(:,i_time),-z_mid(:,1),'-o', 'LineWidth', 1); 
   end
end
xlim([0. 10.]); %ylim([-0.7 0]);


subplot(2,2,2) %--- T ---
hold on; grid on
xlabel('T_i'); ylabel('z frac');
xlim([-2 0]); %ylim([-0.7 0]);
for i_time = 1:Nts
   if ( days(i_time) == 1.   | ...
        days(i_time) == 20.  | ...
        days(i_time) == 40.  | ...
        days(i_time) == 60.  | ...
        days(i_time) == 80.  | ...
        days(i_time) == 100. )
      plot(T_i(:,i_time)-Kelvin,-z_mid(:,1), 'LineWidth', 1); 
   end
end


subplot(2,2,3) %--- w ---
hold on; grid on
xlabel('w (m/s)'); ylabel('z frac')
for i_time = 1:Nts
   if ( days(i_time) == 1.   | ...
        days(i_time) == 20.  | ...
        days(i_time) == 40.  | ...
        days(i_time) == 60.  | ...
        days(i_time) == 80.  | ...
        days(i_time) == 100. )
      plot(w_br(:,i_time),-z_mid(:,1), 'LineWidth', 1); 
   end
end

subplot(2,2,4) %--- D ---
hold on; grid on
xlabel('\phi'); ylabel('z frac')
for i_time = 1:Nts
   if ( days(i_time) == 1.   | ...
        days(i_time) == 20.  | ...
        days(i_time) == 40.  | ...
        days(i_time) == 60.  | ...
        days(i_time) == 80.  | ...
        days(i_time) == 100. )
      plot(f_br(:,i_time),-z_mid(:,1), 'LineWidth', 1); 
   end
end
legend('1d', '20d', '40d', '60d', '80d', '100d','location','best' )
legend('boxoff')


% title
ax = [0.3 .92 0.45 0.02];
axes('position',ax);
title([LeTitre])
axis off;

% prints
if iprint==1
   printfig('/usr/home/crlod/WORK/NEMO/PROG/fig_local/SALT/',...
            ['flushing_profiles_scheme',num2str(i_scheme),'_liq',num2str(i_liq),'_lay',num2str(Nlays)])
end

