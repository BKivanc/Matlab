clc;
clear;
% Reference Aircraft - Cessna CitationJet-M2 (CJ1)
% Name of the Designed Aircraft: Apex-HBH0705
format long
Sref=22.3; %do not change while iteration
bref=14.4; %do not change while iteration
W0_guess = 31629; %change
W_crew=2*900;
W_PL=(5*900)+2000;
alt_cruise = 38000; %change
alt_loiter = 7000;
b = 18.72; %change
S = 18.003; %change
g=9.80665;
T0 = 8000; %change
tmax = 45/60;
RhVmax = 2100000;
CD0 = 0.016;
AR = (b^2)/S;
e = 0.8;
K = (1/(pi*e*AR));
[temp, sos, press, rho, rhoh] = atmos_1976_SI(alt_cruise);
CL_star = sqrt(CD0/K);
LD = 1/(2*sqrt(CD0*K));
W1 = 0.995*W0_guess;
W2 = (1-0.000001*alt_cruise)*W1;
VmaxR = sqrt((2*W2)/(rho*S*sqrt(CD0/(3*K))));
[T_cruise,ct_cruise] = turbofan_SI(T0,VmaxR,alt_cruise);
CL2 = (2*W2)/(rho*(VmaxR^2)*S);
A=((RhVmax*ct_cruise*(1/3600))/(2*VmaxR*LD));
B= (atan(CL2/CL_star))-A;
CL3=CL_star*tan(B);
W3 = (CL3*rho*(VmaxR^2)*S)/2;
W4 = 0.993*W3;
[temp_loiter, sos_loiter, press_loiter, rho_loiter, rhoh_loiter] = atmos_1976_SI(alt_loiter);
VLDmax=sqrt((2*W4)/(rho_loiter*S*CL_star));
[T_loiter,ct_loiter] = turbofan_SI(T0,VLDmax,alt_loiter);
W5= W4*exp((-tmax*ct_loiter)/LD);
W6 = 0.995*W5;
W_misfuel=W0_guess-W6;
W_total_fuel=1.05*W_misfuel;
W_empty=0.7446*W0_guess^(0.9713);
W0_calculated=W_crew+W_PL+W_total_fuel+W_empty
diff = abs(W0_guess-W0_calculated);

M_cruise = VmaxR/sos;
C = 0.488*M_cruise^(0.728)
T_W_ratio = T0/W0_guess


alt_sealevel = 0;
[temp_sealevel, sos_sealevel, press_sealevel, rho_sealevel, rhoh_sealevel] = atmos_1976_SI(alt_sealevel);
CLmax_landing = 2;


mur_takeoff = 0.03;
CL=0.1;
hb = 0.08;
hOB = 50*0.3048;
CD_gear = 0.02;
CD_flaps = 0.015;
CLmax_takeoff = 2;
G = ((16*hb)^2)/(1+(16*hb)^2);
Vstall_takeoff = sqrt((2*W0_guess) / (rho_sealevel * S * CLmax_takeoff));
VLO = 1.2*Vstall_takeoff;
R = (Vstall_takeoff^2)*4.8/g;
gammaOB = acos((R-hOB)/R);
sa_takeoff = R*sin(gammaOB)
CD = CD0 + CD_gear + CD_flaps + (G*K*CL^2);

coeff_a = (1/W0_guess)*(((1/2)*rho_sealevel*S)*(CD-mur_takeoff*CL));
coeff_t = ((T0/W0_guess)-mur_takeoff);
sg_takeoff = (-1/(2*g*coeff_a))*log((1 - (coeff_a/coeff_t)*VLO^2))
Takeoff_total = sg_takeoff + sa_takeoff

Vstall_landing = sqrt((2*W6)/(rho_sealevel*S*CLmax_landing));
VTD = 1.15*Vstall_landing;
mur_landing = 0.2;
sg_landing = (VTD^2)/(2*g*(mur_landing))


if T_W_ratio<=C
    disp("Thrust-to-Weight condition satisfied")
else 
    disp("Thrust-to-Weightcondition not satisfied")
end

if diff<=0.1
    disp("Weight condition satisfied")
else 
    disp("Weight condition not satisfied")
end

if Takeoff_total<=700
    disp("Takeoff distance satisfied")
else 
    disp("Takeoff distance not satisfied")
end

if sg_landing<=400
    disp("Ground roll distance satisfied")
else 
    disp("Ground roll distance not satisfied")
end

if Sref*0.7<=S & S<=Sref*1.3
    disp("S referance satisfied")
else 
    disp("S referance not satisfied")
end

if bref*0.7<=b & b<=bref*1.3
    disp("b referance satisfied")
else 
    disp("b referance not satisfied")
end

if 20000<=alt_cruise & alt_cruise<=38000
    disp("altitude satisfied")
else 
    disp("altitude not satisfied")
end
