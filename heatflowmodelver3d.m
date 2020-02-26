clear all
clc

% Definitions
% Node g (Regolith)
% Node 0 (Regolith boundary layer beneath greenhouse)
% Node 1 (Isolation in greenhouse floor)
% Node 2 (Regolith in greenhouse floor)
% Node 3 (Air inside greenhouse)
% Node 4 (Inner roof layer)
% Node 5 (Isolating layer)
% Node 6 (Outer roof layer)
% Node 7 (Boundary layer of atmosphere above roof)
% Node 8 (Atmosphere)

% User defined constants
radius=3.0;     % m    Radius of semi-cylindrical roof
Length=1;%57.0;    % m    Length of greenhouse
Ngr=1;%6;          % pcs  Number og greenhouses
pin=100000;     % Pa   Inside greenhouse pressure
Nlayers=2;      % pcs  Number of polyimide layers in roof
pisolate=50000; % Pa   Pressure in isolating layer
%d8=0.01;        % m  Assumed thermally active thickness of atmosphere
d7=0.5;        % m  Thickness of boundary layer between outer layer of roof and atmosphere
dRoof=0.025E-3; % m    Thickness of outer and inner layer of plastic
disolate=0.01;  % m    Thickness of isolating CO2-layer in roof
dFloor=0.30;    % m    Thickness of regolith layer in floor
dIsolationFloor=0.10;  % m Thickness of insolating layer in floor
d0=0.05;        % m    Thickness of boundary layer between greenhouse bottom and regolith.

% Meteorological constants
Tmaxavg=273;  % K
Tminavg=193;  % k
Tground=233;  % K Average temperature of Martian surface. Assumed to be fairly constant.
pout=1000;    % Pa   Outside pressure

% Natural contants
sol=88775;            % seconds (24*3600 + 39*60 + 35)
sigma=5.670E-8;       % W/(m^2*K^4)     Stefan - boltzmann-constant
R=8.314;              % Pa*m^3/(mol*K)  Gas constant
CCO2=37.135;          % J/(K*mol)       Specific heat capacity CO2
CAir=718;             % J/(kg*K)        Cv @300K
CLarc=1090;           % J/(kg*K)
Cregolith=1.2E6;      % J/(m^3*K)       @250 K, Phoenix
CEpolystyrene=1300;   % J/(kg*K)
MCO2=44.01;           % g/mol           Molar mass CO2
densityCO2RTp=1.98;   % kg/m^3          @standard T and p
densityAirRTp=1.225;  % Kg/m^3          @15C and 101300 Pa
densityLarc=1540;     % kg/m^3
densitySoil=1500;     % kg/m^3
densityEPolystyrene=50; % kg/m^3
JsunE=1361;           % W/m^2
Re=1;                 % AU
Rm=1.524;             % AU
Jsun0=JsunE*(Re^2)/(Rm^2);    % W/m^2

% Thermal constants
lambdaCO2200K=9.6E-3;     % W/m*k Heat conductivity of C02 @200K
lambdaCO2300K=16.8E-3;    % W/m*k Heat conductivity of C02 @300K
lambdaAir233K=0.0209;     % W/m*k Heat conductivity of Air @233K
lambdaAir300K=0.0260;     % W/m*k Heat conductivity of Air @300K
lambdaLarc=0.12;          % W/m*k Heat conductivity of NeXolve LaRC-CP1 Polyimide
lambdaRegolith=0.12;      % W/m*k Heat conductivity of Martian regolith @250K, results from Phoenix
lambdaIsolation=0.033;   % W/m*k Heat conductivity of expanded polystyrene

% Optical constants
avisLarc=0.08; % Visible, 0.025 mm Absorbtance
tvisLarc=0.83; % Visible, 0.025 mm Transmittance
rvisLarc=0.09; % Visible, 0.025 mm Reflectivity
ELaRC=0.45;    % IR, 0.025 mm, hemispherical Emissivity
AirLarcAg=0.1; % Ag-coated IR, GUESS Absorbtance
TirLarcAg=0.2; % Ag-coated IR, GUESS Transmittance
RirLarcAg=0.7; % Ag-coated IR, GUESS Reflectivity
avisFloor=0.8; % Absorbtance
tvisFloor=0.0; % Transmittance
rvisFloor=0.2; % Reflectivity
EFloor=0.95;   % Emissivity

% Derived constants
Ar=Ngr*(Length*pi*radius);                %m^2 Area of roof of greenhouses
Af=Ngr*(Length*2*radius);                 %m^2 Area og floor of greenhouse
densityCO2rel=densityCO2RTp*300/101300;   % Relative density of CO2

m7=Ngr*d7*Length*densityCO2rel*pout/Tminavg; % Mass of layer 7
n7=m7*1000/MCO2;                          % Number of moles of layer 7
m6=Ngr*Ar*dRoof*densityLarc;              %kg Mass of outer polyimide layer
m5=Ngr*disolate*Length*densityCO2rel*pisolate/Tminavg; % Mass of layer 5
n5=m5*1000/MCO2;                          % Number of moles of layer 5
m4=m6;                                    % kg Mass layer 4
m3=Ngr*(1/2)*pi*radius^2*Length*densityAirRTp; % Mass of layer 3
V2=Ngr*2*radius*Length*dFloor;            % m^3 Volume layer 2
m1=Ngr*2*radius*Length*dIsolationFloor*densityEPolystyrene;
V0=Ngr*2*radius*Length*d0;

C0=V0*Cregolith;
C1=m1*CEpolystyrene;    %Heat Capacities (J/K)
C2=V2*Cregolith;
C3=m3*CAir;
C4=m4*CLarc;
C5=n5*CCO2;
C6=m6*CLarc;
C7=n7*CCO2;

h0=Af*lambdaRegolith/d0;                % Thermal conductance of layer 0
h1=Af*lambdaIsolation/dIsolationFloor;  % Thermal conductance of layer 1
h2=Af*lambdaRegolith/dFloor;            % Thermal conductance of soil
h3=Ar*lambdaAir300K/radius;             % Thermal conductance of layer 3
h4=Ar*lambdaLarc/dRoof;                  % Thermal conductance of layer 4
h6=Ar*lambdaLarc/dRoof;

r24=(ELaRC*EFloor/(ELaRC+EFloor-ELaRC*EFloor))*Af*sigma;
r42=(ELaRC*EFloor/(ELaRC+EFloor-ELaRC*EFloor))*Ar*sigma;
r46=(ELaRC*ELaRC/(ELaRC+ELaRC-ELaRC*ELaRC))*Ar*sigma;
r64=(ELaRC*ELaRC/(ELaRC+ELaRC-ELaRC*ELaRC))*Ar*sigma;



%%%%  Time loop
dt=0.1;             % Step size in seconds (sunrise after 22196 seconds = 0.25*sol)
N=ceil(100/dt);
t=zeros(1,N);

T0=zeros(1,N);
T1=zeros(1,N);
T2=zeros(1,N);
T3=zeros(1,N);
T4=zeros(1,N);
T5=zeros(1,N);
T6=zeros(1,N);
T7=zeros(1,N);
%T8=zeros(1,N);

Tg(1)=Tground;
T0(1)=233;                         %Initial temperatures
T1(1)=2.317296852e+02; 
T2(1)=2.296915976e+02; 
T3(1)=2.293836553e+02; 
T4(1)=2.290858479e+02; 
T5(1)=2.288983311e+02; 
T6(1)=2.287108143e+02; 
T7(1)=2.110571376e+02; 
T8(1)=Tminavg; 

for i=1:N
% Update time
t(i+1)=t(i)+dt;

% Calculate solar input
theta(i)=2*pi*t(i)/sol-pi/2;                            % rad Solar elevation angle
Qsun(i)=0;%Jsun0*(1+sin(theta(i)))*Length*radius*Ngr;      % w/m^2 The amount of sunlight that reach the roof. (dA = phistep*Length*radius.)
if sin(theta(i)) < 0 
 Qsun(i)=0;
endif

% Calculate various time dependent parameters
T8(i)=Tminavg;%(Tmaxavg-Tminavg)*0.5*sin(2*pi*t(i)/sol-pi/2) - (Tmaxavg-Tminavg)*0.5+273.15; %Approximate diournal temperature oscillation

lambdaCO25(i)=lambdaCO2200K;%((lambdaCO2300K-lambdaCO2200K)/(300-200)*(T5(i)-300) + lambdaCO2300K); %Linear interpolation of lambda CO2
lambdaCO27(i)=lambdaCO2200K;%((lambdaCO2300K-lambdaCO2200K)/(300-200)*(T7(i)-300) + lambdaCO2300K); %Linear interpolation of lambda CO2
%lambdaCO28(i)=lambdaCO2200K;%((lambdaCO2300K-lambdaCO2200K)/(300-200)*(T8(i)-300) + lambdaCO2300K); %Linear interpolation of lambda CO2

h5(i)=Ar*lambdaCO25(i)/disolate;    % Thermal conductance of layer 5
h7(i)=Ar*lambdaCO27(i)/d7;          % Thermal conductance of layer 7
%h8(i)=Ar*lambdaCO28(i)/d8;          % Thermal conductance of layer 8

Qa2(i)=Qsun(i)*tvisLarc^Nlayers*avisFloor;      % Absorptance of visible light in regolith

% Runge Kutta solution (Non-Linearized version):
% Define funktions [dT/dt = f(t, T0, T1, ...) = Q1 + Q2 + ...]
fg=@(t,Tg      )        0;
f0=@(t,   T0,T2)        (            - h1*(T0-T2)                              )/C0;
f1=@(t,Tg,T1,T3)        (-h0*(T1-Tg) - h2*(T1-T3)                              )/C1;
f2=@(t,T0,T2,T4)        (-h1*(T2-T0) - h3*(T2-T4) - r24*(T2^4 - T4^4) + Qa2(i) )/C2;
f3=@(t,T1,T3,T5)        (-h2*(T3-T1) - h4*(T3-T5)                              )/C3;
f4=@(t,T2,T4,T6)        (-h3*(T4-T2) - h5(i)*(T4-T6) - r42*(T4^4 - T2^4) - r46*(T4^4 - T6^4))/C4;
f5=@(t,T3,T5,T7)        (-h4*(T5-T3) - h6*(T5-T7)                              )/C5;
f6=@(t,T4,T6,T8)        (-h5(i)*(T6-T4) - h7(i)*(T6-T8) - r64*(T6^4 - T4^4)    )/C6;
f7=@(t,T5,T7   )        (-h6*(T7-T5)                                           )/C7;
f8=@(t,      T8)        0;

% RK step 1
k1g=dt*fg(t(i),Tg(i));
k10=dt*f0(t(i),T0(i),T2(i));
k11=dt*f1(t(i),Tg(i),T1(i),T3(i));
k12=dt*f2(t(i),T0(i),T2(i),T4(i));
k13=dt*f3(t(i),T1(i),T3(i),T5(i));
k14=dt*f4(t(i),T2(i),T4(i),T6(i));
k15=dt*f5(t(i),T3(i),T5(i),T7(i));
k16=dt*f6(t(i),T4(i),T6(i),T8(i));
k17=dt*f7(t(i),T5(i),T7(i));
k18=dt*f8(t(i),T8(i));

% RK step 2
k2g=dt*fg(t(i)+dt/2,Tg(i)+k1g/2);
k20=dt*f0(t(i)+dt/2,T0(i)+k10/2,T2(i)+k12/2);
k21=dt*f1(t(i)+dt/2,Tg(i)+k1g/2,T1(i)+k11/2,T3(i)+k13/2);
k22=dt*f2(t(i)+dt/2,T0(i)+k10/2,T2(i)+k12/2,T4(i)+k14/2);
k23=dt*f3(t(i)+dt/2,T1(i)+k11/2,T3(i)+k13/2,T5(i)+k15/2);
k24=dt*f4(t(i)+dt/2,T2(i)+k12/2,T4(i)+k14/2,T6(i)+k16/2);
k25=dt*f5(t(i)+dt/2,T3(i)+k13/2,T5(i)+k15/2,T7(i)+k17/2);
k26=dt*f6(t(i)+dt/2,T4(i)+k14/2,T6(i)+k16/2,T8(i)+k18/2);
k27=dt*f7(t(i)+dt/2,T5(i)+k15/2,T7(i)+k17/2);
k28=dt*f8(t(i)+dt/2,T8(i)+k18/2);

% RK step 3
k3g=dt*fg(t(i)+dt/2,Tg(i)+k2g/2);
k30=dt*f0(t(i)+dt/2,T0(i)+k20/2,T2(i)+k22/2);
k31=dt*f1(t(i)+dt/2,Tg(i)+k2g/2,T1(i)+k21/2,T3(i)+k23/2);
k32=dt*f2(t(i)+dt/2,T0(i)+k20/2,T2(i)+k22/2,T4(i)+k24/2);
k33=dt*f3(t(i)+dt/2,T1(i)+k21/2,T3(i)+k23/2,T5(i)+k25/2);
k34=dt*f4(t(i)+dt/2,T2(i)+k22/2,T4(i)+k24/2,T6(i)+k26/2);
k35=dt*f5(t(i)+dt/2,T3(i)+k23/2,T5(i)+k25/2,T7(i)+k27/2);
k36=dt*f6(t(i)+dt/2,T4(i)+k24/2,T6(i)+k26/2,T8(i)+k28/2);
k37=dt*f7(t(i)+dt/2,T5(i)+k25/2,T7(i)+k27/2);
k38=dt*f8(t(i)+dt/2,T8(i)+k28/2);

% RK step 4
k4g=dt*fg(t(i)+dt,Tg(i)+k3g);
k40=dt*f0(t(i)+dt,T0(i)+k30,T2(i)+k32);
k41=dt*f1(t(i)+dt,Tg(i)+k3g,T1(i)+k31,T3(i)+k33);
k42=dt*f2(t(i)+dt,T0(i)+k30,T2(i)+k32,T4(i)+k34);
k43=dt*f3(t(i)+dt,T1(i)+k31,T3(i)+k33,T5(i)+k35);
k44=dt*f4(t(i)+dt,T2(i)+k32,T4(i)+k34,T6(i)+k36);
k45=dt*f5(t(i)+dt,T3(i)+k33,T5(i)+k35,T7(i)+k37);
k46=dt*f6(t(i)+dt,T4(i)+k34,T6(i)+k36,T8(i)+k38);
k47=dt*f7(t(i)+dt,T5(i)+k35,T7(i)+k37);
k48=dt*f8(t(i)+dt,T8(i)+k38);

% RK Iteration
Tg(i+1)=Tg(i);
T0(i+1)=T0(i) + (1/6)*(k10 + 2*k20 + 2*k30 + k40);
T1(i+1)=T1(i) + (1/6)*(k11 + 2*k21 + 2*k31 + k41);
T2(i+1)=T2(i) + (1/6)*(k12 + 2*k22 + 2*k32 + k42);
T3(i+1)=T3(i) + (1/6)*(k13 + 2*k23 + 2*k33 + k43);
T4(i+1)=T4(i) + (1/6)*(k14 + 2*k24 + 2*k34 + k44);
T5(i+1)=T5(i) + (1/6)*(k15 + 2*k25 + 2*k35 + k45);
T6(i+1)=T6(i) + (1/6)*(k16 + 2*k26 + 2*k36 + k46);
T7(i+1)=T7(i) + (1/6)*(k17 + 2*k27 + 2*k37 + k47);
T8(i+1)=T8(i);

endfor

%%%%%%% Plots
%T6(i+1)=T6(i);
%T5(i+1)=T5(i);
%T4(i+1)=T4(i);
%T3(i+1)=T3(i);
%T2(i+1)=T2(i);
%T1(i+1)=T1(i);

%t(i+1)=t(i)+dt;
%T8(i+1)=T8(i);
%Qsun(i+1)=Qsun(i);
%normdiff(i+1)=normdiff(i);
%theta(i+1)=theta(i);
%SSiterations(i+1)=SSiterations(i);
%Qr24(i+1)=Qr24(i);
%Qr42(i+1)=Qr42(i);
%Qc21(i+1)=Qc21(i);
%Qc10(i+1)=Qc10(i);

%plot(t,T6)
%grid on

