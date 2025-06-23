% Matlab program for fiber optic SPR sensor
% Gupta, B. D., and Anuj K. Sharma. "Sensitivity evaluation of a 
% multi-layered surface plasmon resonance-based fiber optic sensor: 
% a theoretical study." Sensors and Actuators B: Chemical 107, no. 1 
% (2005): 40-46.

% Parameters to be varied 
% 1. Incident angle (theta)
% 2. Wavelength (lmda)
% 3. Core diameter (D)
% 4. Sensing length (L)
% 5. Material refractive index (n)
% 6. Matherial thickness (d)

clc
clear all
close all

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 10)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 10)



c=2.998*1e08;
Dcore=600e-6; % Core diameter of fiber

d0=0;       % Thickness of substrate
d1=7.5*3*1e-9; % Thickness of silver, d1+d2=30;    % Thickness of silver (d1) + gold (d2)
d2=d1/3;    % Thickness of gold
d_sens=0;   % Thickness of sensing medium

NumA=0.3; % Numerical aperture of fiber
lmda_max=600*1e-9;   % Absoption maximum wavelength
w0=(2*pi*c)/lmda_max; % Absorption frequency
D_lmda_max=75*1e-9; % Full width at half-maximum of absorption spectrum of the sensing layer sample
gama=(2*pi*c)*(D_lmda_max/lmda_max^2);
NA=6.023*1e23; % Avogadro number
e=1.602*1e-19; %Charge of electron
me=9.109*1e-31; % Mass of electron
e_sens_inf= (1.341)^2; % Background dielectric constant of the sensing layer
e_0=8.8541878128*1e-12; % Permittivity of vaccum
f=1.0;  % Oscillator strength (function of molar extinction coefficient)


index2 = 0;
for x=-0.75:0.25:0.5  
index2 = index2 + 1;
w=w0-x*gama;
lmda=(2 * pi*c)/w;
lmda1=lmda*1e+6;
k0 = w/c;

% Refractive index of core for lmda=0.5 to 1.6 micrometer
l=0.035;
C0=1.4508554;
C1=-0.0031268;
C2=-0.0000381;
C3=0.0030270;
C4=-0.0000779;
C5=0.0000018;
n_core=C0+(C1*lmda1^2)+(C2*lmda1^4)+(C3/((lmda1^2)-l))+(C4/((lmda1^2)-l)^2)+(C5/((lmda1^2)-l)^3);

theta_cr= acos(NumA/n_core); % Criticl angle

% Refractive index of plasmonic metal 
lmdap1=1.4541e-7*1e+6; % for silver
lmdac1=1.7614e-5*1e+6; % for silver
n_Ag=sqrt(1-(((lmda1^2)*lmdac1)/((lmdap1^2)*(lmdac1+i*lmda1))));

lmdap2=1.6826e-7*1e+6; % for gold
lmdac2=8.9342e-6*1e+6; % for gold
n_Au=sqrt(1-(((lmda1^2)*lmdac2)/((lmdap2^2)*(lmdac2+i*lmda1))));
        
index1 = 0;
for Con=0.005:0.005:0.055    % C=[0.005 0.055] molar concentration of absorption oscillators

index1 = index1 + 1;

% Refractive index of sensing (absorbing) layer
N=1e3*NA*Con; 
e_sens=e_sens_inf+((N*e^2)/(me*e_0 ))*(f/(w0^2-w^2-i*w*gama));
n_sens=sqrt(e_sens);

d=[d0 d1 d2 d_sens];
n=[n_core n_Ag n_Au n_sens];
 
dtheta1=1e-5;
index0 = 0;
for theta1= theta_cr:dtheta1:pi/2;   % theta1_max=acosd(NumA/n_core)
index0 = index0 + 1;

L=30*Dcore;   % Length of the exposed sensing region
Nref=L/(Dcore*tan(theta1)); %Total number of reflections performed by a ray making angle ? with the normal to the core-metal layer interface in the sensing region

Pnum=(n_core^2)*(sin(theta1))*(cos(theta1));
Pden=(1-(n_core^2)*(cos(theta1))^2)^2;
P1=Pnum/Pden;

M=length(n); % Total Number of Layers in the Sensors

% Calculate propagation angles within each layer
master_matrix=eye(2);
for q=1:M
if q>=2
costheta(q) = sqrt(1 - (n(1)/n(q))^2 * (sin(theta1))^2);
end
% Calulate q terms for reflection at boundaries
if q==1
react(q) = (n(q) * k0 * cos(theta1));
else
react(q) = (n(q) * k0 * costheta(q));
end
% Calculation of layer terms
Z(q) = (react(q)) / (w * n(q)^2);
if q>=2||q<M
qprod(q) = d(q) * react(q);
A(q) = cos(qprod(q));
B(q) = (sin(qprod(q)) * (-j * Z(q)));
C(q) = ((-j * sin(qprod(q))) / Z(q));
D(q) = cos(qprod(q));
% Solve for master matrix values
matrix = [A(q) , B(q); C(q) , D(q)];
master_matrix = master_matrix*matrix; 
end
end
A = master_matrix(1,1);
B = master_matrix(1,2);
C = master_matrix(2,1);
D = master_matrix(2,2);
% Calculate the reflection value  
Rnum = A + (B/Z(M)) - Z(1)*(C+(D/Z(M)));
Rden = A + (B/Z(M)) + Z(1)*(C+(D/Z(M)));
R = (Rnum/Rden);
R=abs(R);
R = R^2;

P2=(R^Nref)*P1;

P11(index0)=P1;
P21(index0)=P2;
Theta(index0)=theta1;
end  
Tnum=trapz(Theta,P21);
Tden=trapz(Theta,P11);

% Tnum=dtheta1*trapz(P21);
% Tden=dtheta1*trapz(P11);

T=Tnum/Tden;

if index1==1
    Ti=T;
    Ci=Con;
end
S=abs((T-Ti)/(Con-Ci));
Sc(index1) = S;
Cc(index1)=Con;
Pt(index1)=T;
end
figure(1)
plot(Cc,Sc)
xlabel('Concentration (M)')
ylabel('Sensitivity (a.u.)')
hold on
figure(2)
plot(Cc,Pt)
xlabel('Concentration (M)')
ylabel('Normalized transmitted power(a.u.)')
hold on
end
legend('\omega = \omega_0- 0.75\gamma','\omega = \omega_0- 0.50\gamma',...
    '\omega=\omega_0- 0.25\gamma','\omega = \omega_0',...
    '\omega = \omega_0+ 0.25\gamma','\omega = \omega_0+ 0.50\gamma')
  
