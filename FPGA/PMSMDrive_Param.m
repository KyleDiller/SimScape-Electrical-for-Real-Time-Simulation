
% Control Parameters

%Ts_sim        =100e-9*1e+0/2       ; % Calculation Time [sec]
%Ts_sim=0.5e-6;
Tsc           =50e-6;      ; 
exp_temp     =   20.0                 ; % Temparature [deg]
Inertia_rate =  0.157               ; % Inertia rate [times]


% Motor Parameters

Power       = 7500                      ; % Rated Power [W]
Pole        = 6                         ; % Poles
Wrpm_100    = 1800                      ; % Based speed [r/min]
Wrpm_max    = 2700                      ; % MAX speed [r/min]
Jf          = 124.1060+20.5             ; % Inertia [10-4*kgm2]
R1_20       = 0.120                     ; % Resistance (Phase) [ohm]
temp        = 20                        ; % temp. [deg] 
Irms        = 26.00                     ; % Rated Current [Arms]
Turn        = 26                        ; % Coil Turns [Turn]
Ld          = 2.984                      ; % d-Inductance (Phase) [mH]
Lq          = 4.576                      ; % q-Inductance (Phase) [mH]
KE          = 97.5988        ; % Inducted Voltage constant [mV(rms)/rpm]


% Vector Control

Torque_Lim  = 50           ; % Rate of the Torque Current Limit [%]
Weak_Lim    = 50          ; % Rate of the Weak Current Limit [%]
wsc         =  200          ; % Frequency Response of the Speed Controller
wcc         = 2000          ; % Response Frequency of the Current Controller
carrier     = 9000          ; % Carrier Frequency
wspi        = wsc / 10      ; % Break Frequency of the Speed Controller
wcpi        = wcc / 10      ; % Break Frequency of the Current Controller


% Fundamenal Parameters of the Circuit (Indispensable)           
Ron_T       = 0.010         ;
Vf_T        = 1.000         ;
VDC         = 200 * sqrt(2)     ;
V_upper     = 400               ;
AC_freq     = 60            ;
PN_Cap      = 3000e-6       ;
Td          = 5.0e-6        ;


% Useful Constant Values (Auto)
sqrt2       = sqrt( 2 ) ;
sqrt3       = sqrt( 3 ) ;
sqrt23      = sqrt( 2 / 3 ) ;
tr32        = sqrt23 * [ 1 -1/2 -1/2 ; 0 sqrt(3)/2 -sqrt(3)/2 ];
tr23        = sqrt23 * [ 1 0 ; -1/2 sqrt(3)/2 ; -1/2 -sqrt(3)/2 ];
rad_to_rpm  = 60 / 2 / pi ;
rpm_to_rad  = 1 / 60 * 2 * pi ;


% Automatic calculation (Auto)

Pm          = Pole / 2 ;
Wr_100      = Wrpm_100 * rpm_to_rad ;
Wr_max      = Wrpm_max * rpm_to_rad ;
R1          = ( 234.5 + exp_temp ) / ( 234.5 + temp ) * R1_20 ;
Ld_100      = Ld / 1000 ;
inv_Ld_100  = 1 / Ld_100 ;
Lq_100      = Lq / 1000 ;
inv_Lq_100  = 1 / Lq_100 ;
Ke          = KE / 1000 * 60 / 2 / pi ; 
Fai         = Ke / Pm ;
Tm_100      = Power / Wr_100 ;
Tm_max      = Tm_100 * Torque_Lim / 100 ;
iq_100      = Irms * sqrt3 ;
inv_iq_100  = 1 / iq_100 ;
iq_max      = iq_100 * Torque_Lim / 100 ;
id_100      = iq_100 ;
inv_id_100  = 1 / id_100 ;
id_max      = id_100 * Weak_Lim / 100 ;
Irms        = sqrt( id_100^2 + iq_100^2 ) / sqrt(3) ;
Irms_max    = sqrt( id_max^2 + iq_max^2 ) / sqrt(3) ;
Irms_peak   = Irms_max * sqrt(2) ;
Kt          = Ke ;
inv_Kt      = 1 / Kt ;
J1          = Jf / 10000 ;
J2          = J1 * Inertia_rate ;
Jall        = J1 + J2 ;
inv_Jall    = 1 / Jall ;
VDC_2       = VDC / 2 ;


% TMC Control

alpha = 0.866 ;                       ; % Limitter - alpha
beta = 84.00 ;                        ; % Limitter - beta
Kd_TMC = 200 ;           ; % Gain - Kd
Kq_TMC = 40 ;            ; % Gain - Kq
