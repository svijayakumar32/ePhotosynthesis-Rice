This folder contains functions for calculating Rubisco-related parameters, which are then integrated within the e-Photosynthesis model by modifying `CM_Rate.m` as follows:

## Calculation of temperature-dependent Rubisco kinetics using Arrhenius equations
c_c = 38.9;% Kc c rice, umol mol-1

dHa_c = 83.1;% Kc dHa rice, umol mol-1

c_air = 30.5;% Kcair c rice,umol mol-1

dHa_air = 60.5;% Kcair dHa rice, umol mol-1

O2_uM_25  =  210 * 0.0013 * 1000;

PsKM11_0  =  0.0080;                           % initial rice Kc	(liquid phase)

PsKM12_0    =  (O2_uM_25/((17.3/8.0)-1))/1000; % initial rice Ko (liquid phase)	 

Kc_25 = exp(c_c-dHa_c * 1000/(R * (25+273.15)));%Kc at 25C in umol mol -1

Kcair_25 = exp(c_air-dHa_air * 1000/(R * (25+273.15))); %Kcair at 25C in umol mol -1

Ko_25 = ((21000/((Kcair_25/Kc_25) - 1))/1000)*10; %Ko at 25C in mbar or mmol mol -1 

Kc_Tp = exp(c_c-dHa_c * 1000/(R*(Tp+273.15))); %Kc at Tp in umol mol-1

Kcair_Tp = exp(c_air-dHa_air * 1000/(R * (Tp+273.15))); %Kcair at Tp in umol mol -1

Ko_Tp = ((21000/((Kcair_Tp/Kc_Tp) - 1))/1000) * 10; %Ko at Tp in mbar or mmol mol -1 

PsKM11 = PsKM11_0 * (Kc_Tp)/(Kc_25); % temperature-corrected rice Kc	(liquid phase)

PsKM12 = PsKM12_0 * (Ko_Tp)/(Ko_25); % temperature-corrected rice Kc	(liquid phase)

## Temperature response of Rubisco activation state using quadratic polynomial fit (Fig 4A, Makino and Sage, 2007) 

- see `Temp-resp-Rubisco.R` for details)

Ru_Act = - 0.0007851 * Tp ^2 + 0.0319 * Tp + 0.5233;

## Temperature response of ratio between Rubisco oxygenation and carboxylation (PrV111/PsV1 or Vomax/Vcmax) 

PsV1 * Vo/Vc ratio for rice at 25C (Makino et al 1988) * Arrhenius function substituting c = 0.1217 and dHa = -5.22 from a normalised non-linear fit of Bernacchi Vo/Vc against temperature 

- see `Vmax_temp_adj.R` for details

PrV111 = PsV1 * (0.58/1.77)*(0.1217 * exp(-(-5.22) / (0.008314 * (Tp + 273.15))));
