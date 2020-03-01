clear all
clc

% User defined constants
dt=1;           % seconds
radius=3.0;     % m    Radius of semi-cylindrical roof
Length=1;%57.0;    % m    Length of greenhouse
Ngr=1;%6;          % pcs  Number og greenhouses
pin=100000;     % Pa   Inside greenhouse pressure
Nlayers=2;      % pcs  Number of polyimide layers in roof
pisolate=1000; % Pa   Pressure in isolating layer
d8=0.02;        % m  Assumed thermally active thickness of atmosphere
d7=0.02;        % m  Thickness of boundary layer between outer layer of roof and atmosphere
dRoof=0.025E-3; % m    Thickness of outer and inner layer of plastic
disolate=0.02;  % m    Thickness of isolating CO2-layer in roof
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


%%%% Time loop Sunrise: after 22196 seconds
t=[0*dt:dt:6000*dt];%[0:dt:0.01*sol]; %seconds
T0=Tground; %Initial temperature layer 0
T1(1)=224.2105189; %Initial temperature layer 1
T2(1)=210.1087141; %Initial temperature layer 2
T3(1)=207.2326940; %Initial temperature layer 3
T4(1)=204.4513295; %Initial temperature layer 4
T5(1)=202.7927152; %Initial temperature layer 5
T6(1)=201.1341010; %Initial temperature layer 6
T7(1)=198.4608916; %Initial temperature layer 7
T8=Tminavg; %Initial temperature layer 8

h01=1/((dIsolationFloor/2)/(Af*lambdaIsolation)+(d0/2)/(Af*lambdaRegolith));        % Thermal conductance of layer 0 to 1

h10=1/((dIsolationFloor/2)/(Af*lambdaIsolation)+(d0/2)/(Af*lambdaRegolith));        % Thermal conductance of layer 1 to 0
h12=1/((dIsolationFloor/2)/(Af*lambdaIsolation)+(dFloor/2)/(Af*lambdaRegolith)); % Thermal conductance of layer 1 to 2

h21=1/((dFloor/2)/(Af*lambdaRegolith)+(dIsolationFloor/2)/(Af*lambdaIsolation)); % Thermal conductance of layer 2 to 1
h23=1/((dFloor/2)/(Af*lambdaRegolith)+(radius/2)/(Ar*lambdaAir300K));    % Thermal conductance of layer 2 to 3
k24=(ELaRC*EFloor/(ELaRC+EFloor-ELaRC*EFloor))*Af*sigma;

h32=1/((radius/2)/(Ar*lambdaAir300K)+(dFloor/2)/(Af*lambdaRegolith));  % Thermal conductance of soil
h34=1/((radius/2)/(Ar*lambdaAir300K)+(dRoof/2)/(Ar*lambdaLarc));    % Thermal conductance of layer 3 to 4

k42=(ELaRC*EFloor/(ELaRC+EFloor-ELaRC*EFloor))*Ar*sigma;
h43=1/((dRoof/2)/(Ar*lambdaLarc)+(radius/2)/(Ar*lambdaAir300K));   % Thermal conductance of layer 4 to 3
k46=(ELaRC*ELaRC/(ELaRC+ELaRC-ELaRC*ELaRC))*Ar*sigma;

k64=(ELaRC*ELaRC/(ELaRC+ELaRC-ELaRC*ELaRC))*Ar*sigma;

lambdaCO28=((lambdaCO2300K-lambdaCO2200K)/(300-200)*(T8-300) + lambdaCO2300K); %Linear interpolation of lambda CO2

dt0=0;
dt1=dt/(m1*CEpolystyrene);
dt2=dt/(V2*Cregolith);
dt3=dt/(m3*CAir);
dt4=dt/(m4*CLarc);
dt5=dt/(n5*CCO2);
dt6=dt/(m6*CLarc);
dt7=dt/(n7*CCO2);
dt8=0;

for i=1:length(t)
% Calculate the amount of sunlight that hit the greenhouse
theta(i)=2*pi*t(i)/sol-pi/2; % rad Solar elevation angle
Qsun(i)=0;%Jsun0*(1+sin(theta(i)))*Length*radius*Ngr;%*tvisLarc^Nlayers; % w/m^2 dA = phistep*Length*radius. The amount of sunlight that is transmitted through the roof.
if sin(theta(i)) < 0 
 Qsun(i)=0;
endif

lambdaCO27(i)=lambdaCO2200K;%((lambdaCO2300K-lambdaCO2200K)/(300-200)*(T7(i)-300) + lambdaCO2300K); %Linear interpolation of lambda CO2
h78(i)=1/((d7/2)/(Ar*lambdaCO27(i))+(d8/2)/(Ar*lambdaCO28));   % Thermal conductance from layer 7 to 8
h76(i)=1/((d7/2)/(Ar*lambdaCO27(i))+(dRoof/2)/(Ar*lambdaLarc));   % Thermal conductance of layer 7 to 6

h67(i)=h76(i);

lambdaCO25(i)=lambdaCO2200K;%((lambdaCO2300K-lambdaCO2200K)/(300-200)*(T5(i)-300) + lambdaCO2300K); %Linear interpolation of lambda CO2
h65(i)=1/((dRoof/2)/(Ar*lambdaLarc)+(disolate/2)/(Ar*lambdaCO25(i)));   % Thermal conductance of layer 6 to 5

h56(i)=h65(i);

h54(i)=1/((disolate/2)/(Ar*lambdaCO25(i))+(dRoof/2)/(Ar*lambdaLarc)); % Thermal conductance of layer 5 to 4

h45(i)=h54(i);

Qa2(i)=Qsun(i)*tvisLarc^Nlayers*avisFloor;                                  % Absorptance of visible light in regolith

% Node 0 (Regolith beneath greenhouse)
% Node 1 (Isolation in greenhouse floor)
% Node 2 (Regolith in greenhouse floor)
% Node 3 (Air inside greenhouse)
% Node 4 (Inner roof layer)
% Node 5 (Isolating layer)
% Node 6 (Outer roof layer)
% Node 7 (Boundary layer of atmosphere above roof)
% Node 8 (Atmosphere)

Ttemp1=T1(i);
Ttemp2=T2(i);
Ttemp3=T3(i);
Ttemp4=T4(i);
Ttemp5=T5(i);
Ttemp6=T6(i);
Ttemp7=T7(i);

%%%%% RK - step 1


H0=[1,0,0,0,0,0,0,0,0];
H1=[-h10,h10+h12-1/dt1,-h12,0,0,0,0,0,0];
H2=[0,-h21,h23+h21+k24*4*T2(i)^3-1/dt2,-h23,-k24*4*T4(i)^3,0,0,0,0];
H3=[0,0,-h32,h32+h34-1/dt3,-h34,0,0,0,0];
H4=[0,0,-k42*4*T2(i)^3,-h43,h43+h45(i)+k42*4*T4(i)^3+k46*4*T4(i)^3-1/dt4,-h45(i),-k46*4*T6(i)^3,0,0];
H5=[0,0,0,0,-h54(i),h54(i)+h56(i)-1/dt5,-h56(i),0,0];
H6=[0,0,0,0,-k64*4*T4(i)^3,-h65(i),h65(i)+h67(i)+k64*4*T6(i)^3-1/dt6,-h67(i),0];
H7=[0,0,0,0,0,0,-h76(i),h76(i)+h78(i)-1/dt7,-h78(i)];
H8=[0,0,0,0,0,0,0,0,1];

K=[T0;0-T1(i)/dt1;k24*(3*T2(i)^4-3*T4(i)^4)-T2(i)/dt2+Qa2(i);0-T3(i)/dt3;k42*(3*T4(i)^4-3*T2(i)^4)+k46*(3*T4(i)^4-3*T6(i)^4)-T4(i)/dt4;0-T5(i)/dt5;k64*(3*T6(i)^4-3*T4(i)^4)-T6(i)/dt6;0-T7(i)/dt7;T8];

H=[H0;H1;H2;H3;H4;H5;H6;H7;H8];

Tk1=H^-1*K - [T0;T1(i);T2(i);T3(i);T4(i);T5(i);T6(i);T7(i);T8];

%%%%% RK - step 2

T1(i)=Ttemp1+Tk1(2)/2;
T2(i)=Ttemp2+Tk1(3)/2;
T3(i)=Ttemp3+Tk1(4)/2;
T4(i)=Ttemp4+Tk1(5)/2;
T5(i)=Ttemp5+Tk1(6)/2;
T6(i)=Ttemp6+Tk1(7)/2;
T7(i)=Ttemp7+Tk1(8)/2;

dt0=0.5*0;
dt1=0.5*dt/(m1*CEpolystyrene);
dt2=0.5*dt/(V2*Cregolith);
dt3=0.5*dt/(m3*CAir);
dt4=0.5*dt/(m4*CLarc);
dt5=0.5*dt/(n5*CCO2);
dt6=0.5*dt/(m6*CLarc);
dt7=0.5*dt/(n7*CCO2);
dt8=0.5*0;

H0=[1,0,0,0,0,0,0,0,0];
H1=[-h10,h10+h12-1/dt1,-h12,0,0,0,0,0,0];
H2=[0,-h21,h23+h21+k24*4*T2(i)^3-1/dt2,-h23,-k24*4*T4(i)^3,0,0,0,0];
H3=[0,0,-h32,h32+h34-1/dt3,-h34,0,0,0,0];
H4=[0,0,-k42*4*T2(i)^3,-h43,h43+h45(i)+k42*4*T4(i)^3+k46*4*T4(i)^3-1/dt4,-h45(i),-k46*4*T6(i)^3,0,0];
H5=[0,0,0,0,-h54(i),h54(i)+h56(i)-1/dt5,-h56(i),0,0];
H6=[0,0,0,0,-k64*4*T4(i)^3,-h65(i),h65(i)+h67(i)+k64*4*T6(i)^3-1/dt6,-h67(i),0];
H7=[0,0,0,0,0,0,-h76(i),h76(i)+h78(i)-1/dt7,-h78(i)];
H8=[0,0,0,0,0,0,0,0,1];

K=[T0;0-T1(i)/dt1;k24*(3*T2(i)^4-3*T4(i)^4)-T2(i)/dt2+Qa2(i);0-T3(i)/dt3;k42*(3*T4(i)^4-3*T2(i)^4)+k46*(3*T4(i)^4-3*T6(i)^4)-T4(i)/dt4;0-T5(i)/dt5;k64*(3*T6(i)^4-3*T4(i)^4)-T6(i)/dt6;0-T7(i)/dt7;T8];

H=[H0;H1;H2;H3;H4;H5;H6;H7;H8];

Tk2=H^-1*K - [T0;T1(i);T2(i);T3(i);T4(i);T5(i);T6(i);T7(i);T8];

%%%%% RK - step 3

T1(i)=Ttemp1+Tk2(2)/2;
T2(i)=Ttemp2+Tk2(3)/2;
T3(i)=Ttemp3+Tk2(4)/2;
T4(i)=Ttemp4+Tk2(5)/2;
T5(i)=Ttemp5+Tk2(6)/2;
T6(i)=Ttemp6+Tk2(7)/2;
T7(i)=Ttemp7+Tk2(8)/2;


H0=[1,0,0,0,0,0,0,0,0];
H1=[-h10,h10+h12-1/dt1,-h12,0,0,0,0,0,0];
H2=[0,-h21,h23+h21+k24*4*T2(i)^3-1/dt2,-h23,-k24*4*T4(i)^3,0,0,0,0];
H3=[0,0,-h32,h32+h34-1/dt3,-h34,0,0,0,0];
H4=[0,0,-k42*4*T2(i)^3,-h43,h43+h45(i)+k42*4*T4(i)^3+k46*4*T4(i)^3-1/dt4,-h45(i),-k46*4*T6(i)^3,0,0];
H5=[0,0,0,0,-h54(i),h54(i)+h56(i)-1/dt5,-h56(i),0,0];
H6=[0,0,0,0,-k64*4*T4(i)^3,-h65(i),h65(i)+h67(i)+k64*4*T6(i)^3-1/dt6,-h67(i),0];
H7=[0,0,0,0,0,0,-h76(i),h76(i)+h78(i)-1/dt7,-h78(i)];
H8=[0,0,0,0,0,0,0,0,1];

K=[T0;0-T1(i)/dt1;k24*(3*T2(i)^4-3*T4(i)^4)-T2(i)/dt2+Qa2(i);0-T3(i)/dt3;k42*(3*T4(i)^4-3*T2(i)^4)+k46*(3*T4(i)^4-3*T6(i)^4)-T4(i)/dt4;0-T5(i)/dt5;k64*(3*T6(i)^4-3*T4(i)^4)-T6(i)/dt6;0-T7(i)/dt7;T8];

H=[H0;H1;H2;H3;H4;H5;H6;H7;H8];

Tk3=H^-1*K - [T0;T1(i);T2(i);T3(i);T4(i);T5(i);T6(i);T7(i);T8];

%%%%% RK - step 4

T1(i)=Ttemp1+Tk3(2);
T2(i)=Ttemp2+Tk3(3);
T3(i)=Ttemp3+Tk3(4);
T4(i)=Ttemp4+Tk3(5);
T5(i)=Ttemp5+Tk3(6);
T6(i)=Ttemp6+Tk3(7);
T7(i)=Ttemp7+Tk3(8);

H0=[1,0,0,0,0,0,0,0,0];
H1=[-h10,h10+h12-1/dt1,-h12,0,0,0,0,0,0];
H2=[0,-h21,h23+h21+k24*4*T2(i)^3-1/dt2,-h23,-k24*4*T4(i)^3,0,0,0,0];
H3=[0,0,-h32,h32+h34-1/dt3,-h34,0,0,0,0];
H4=[0,0,-k42*4*T2(i)^3,-h43,h43+h45(i)+k42*4*T4(i)^3+k46*4*T4(i)^3-1/dt4,-h45(i),-k46*4*T6(i)^3,0,0];
H5=[0,0,0,0,-h54(i),h54(i)+h56(i)-1/dt5,-h56(i),0,0];
H6=[0,0,0,0,-k64*4*T4(i)^3,-h65(i),h65(i)+h67(i)+k64*4*T6(i)^3-1/dt6,-h67(i),0];
H7=[0,0,0,0,0,0,-h76(i),h76(i)+h78(i)-1/dt7,-h78(i)];
H8=[0,0,0,0,0,0,0,0,1];

K=[T0;0-T1(i)/dt1;k24*(3*T2(i)^4-3*T4(i)^4)-T2(i)/dt2+Qa2(i);0-T3(i)/dt3;k42*(3*T4(i)^4-3*T2(i)^4)+k46*(3*T4(i)^4-3*T6(i)^4)-T4(i)/dt4;0-T5(i)/dt5;k64*(3*T6(i)^4-3*T4(i)^4)-T6(i)/dt6;0-T7(i)/dt7;T8];

H=[H0;H1;H2;H3;H4;H5;H6;H7;H8];

Tk4=H^-1*K - [T0;T1(i);T2(i);T3(i);T4(i);T5(i);T6(i);T7(i);T8];

Tinf=Tk1/6 + Tk2/3 + Tk3/3 + Tk4/6 + [T0;T1(i);T2(i);T3(i);T4(i);T5(i);T6(i);T7(i);T8];


%normdiff(i)=abs(norm(Tinf)-sqrt(T0^2+T1(i)^2+T2(i)^2+T3(i)^2+T4(i)^2+T5(i)^2+T6(i)^2+T7(i)^2+T8^2));

%{
dummy=0;
while normdiff(i) > 1
dummy=dummy+1;

%m*C*dT/dt = Qa
%Qa*dt/m*C = dT


H0=[1,0,0,0,0,0,0,0,0];
H1=[-h10,h10+h12,-h12,0,0,0,0,0,0];
H2=[0,-h21,h23+h21+k24*4*T2(i)^3,-h23,-k24*4*T4(i)^3,0,0,0,0];
H3=[0,0,-h32,h32+h34,-h34,0,0,0,0];
H4=[0,0,-k42*4*T2(i)^3,-h43,h43+h45(i)+k42*4*T4(i)^3+k46*4*T4(i)^3,-h45(i),-k46*4*T6(i)^3,0,0];
H5=[0,0,0,0,-h54(i),h54(i)+h56(i),-h56(i),0,0];
H6=[0,0,0,0,-k64*4*T4(i)^3,-h65(i),h65(i)+h67(i)+k64*4*T6(i)^3,-h67(i),0];
H7=[0,0,0,0,0,0,-h76(i),h76(i)+h78(i),-h78(i)];
H8=[0,0,0,0,0,0,0,0,1];

K=[T0;0;k24*(3*T2(i)^4-3*T4(i)^4)+Qa2(i);0;k42*(3*T4(i)^4-3*T2(i)^4)+k46*(3*T4(i)^4-3*T6(i)^4);0;k64*(3*T6(i)^4-3*T4(i)^4);0;T8];

H=[H0;H1;H2;H3;H4;H5;H6;H7;H8];

Tinf=H^-1*K;

normdiff(i)=abs(norm(Tinf)-sqrt(T0^2+T1(i)^2+T2(i)^2+T3(i)^2+T4(i)^2+T5(i)^2+T6(i)^2+T7(i)^2+T8^2));

T1(i)=Tinf(2); %New temperature layer 1
T2(i)=Tinf(3); %New temperature layer 2
T3(i)=Tinf(4); %New temperature layer 3
T4(i)=Tinf(5); %New temperature layer 4
T5(i)=Tinf(6); %New temperature layer 5
T6(i)=Tinf(7); %New temperature layer 6
T7(i)=Tinf(8); %New temperature layer 7

endwhile
SSiterations(i)=dummy;
%}

T1(i)=Ttemp1;
T2(i)=Ttemp2;
T3(i)=Ttemp3;
T4(i)=Ttemp4;
T5(i)=Ttemp5;
T6(i)=Ttemp6;
T7(i)=Ttemp7;



%T1(i+1)=Tinf(2); %New temperature layer 1
%T2(i+1)=Tinf(3); %New temperature layer 2
%T3(i+1)=Tinf(4); %New temperature layer 3
%T4(i+1)=Tinf(5); %New temperature layer 4
%T5(i+1)=Tinf(6); %New temperature layer 5
%T6(i+1)=Tinf(7); %New temperature layer 6
%T7(i+1)=Tinf(8); %New temperature layer 7


%{
konst=1.2E-3;

if abs(Tinf(2)-T1(i)) < konst
  T1(i+1)=T1(i);
else
  T1(i+1)=(Tinf(2)+T1(i))/2;
endif

if abs(Tinf(3)-T2(i)) < konst
  T2(i+1)=T2(i);
else
  T2(i+1)=(Tinf(3)+T2(i))/2;
endif

if abs(Tinf(4)-T3(i)) < konst
  T3(i+1)=T3(i);
else
  T3(i+1)=(Tinf(4)+T3(i))/2;
endif

if abs(Tinf(5)-T4(i)) < konst
  T4(i+1)=T4(i);
else
  T4(i+1)=(Tinf(5)+T4(i))/2;
endif

if abs(Tinf(6)-T5(i)) < konst
  T5(i+1)=T5(i);
else
  T5(i+1)=(Tinf(6)+T5(i))/2;
endif

if abs(Tinf(7)-T6(i)) < konst
  T6(i+1)=T6(i);
else
  T6(i+1)=(Tinf(7)+T6(i))/2;
endif

if abs(Tinf(8)-T7(i)) < konst
  T7(i+1)=T7(i);
else
  T7(i+1)=(Tinf(8)+T7(i))/2;
endif
%}

T1(i+1)=(Tinf(2)+T1(i))/2; %New temperature layer 1
T2(i+1)=(Tinf(3)+T2(i))/2; %New temperature layer 2
T3(i+1)=(Tinf(4)+T3(i))/2; %New temperature layer 3
T4(i+1)=(Tinf(5)+T4(i))/2; %New temperature layer 4
T5(i+1)=(Tinf(6)+T5(i))/2; %New temperature layer 5
T6(i+1)=(Tinf(7)+T6(i))/2; %New temperature layer 6
T7(i+1)=(Tinf(8)+T7(i))/2; %New temperature layer 7

%Still need radiative contributions
% Modellen i bogen for to lag, der udveksler strålevarme med hinanden tager ikke højde for, at der 
% også er en del af strålingen, der bliver transmitteret. Nok fordi det sædvanligvis ikke erf
% relevant inden for rumfart. LaRC Polyimide er derimod delvist transparent i IR, men det kan 
% vi tilsyneladende stort set undgå ved at bruge en metalcoating. Så det vil bare sige at emissiviteten
% = absorbitviteten for taget skal vi bruge den for coated materiale.

endfor

%%%%%%% Plots
%T6(i+1)=T6(i);
%T5(i+1)=T5(i);
%T4(i+1)=T4(i);
%T3(i+1)=T3(i);
%T2(i+1)=T2(i);
%T1(i+1)=T1(i);

t(i+1)=t(i)+dt;
%T8(i+1)=T8(i);
Qsun(i+1)=Qsun(i);
%normdiff(i+1)=normdiff(i);
theta(i+1)=theta(i);
%SSiterations(i+1)=SSiterations(i);
%Qr24(i+1)=Qr24(i);
%Qr42(i+1)=Qr42(i);
%Qc21(i+1)=Qc21(i);
%Qc10(i+1)=Qc10(i);

plot(t,T3)
%grid on

