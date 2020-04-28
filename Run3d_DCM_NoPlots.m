function [ChiHa,PsiHa,KI_A,KIII_A,Kb0] = Run3d_DCM (data,H,x0,y0)
L = 2.5e3;

x = data(:,1);
y = data(:,2);
z = data(:,3);
u = data(:,4);
v = data(:,5);
w = data(:,6);

upper_A = find( (y > y0) & (x<x0) );
lower_A = find( (y < y0) & (x<x0) );
upper_B = find( (y > y0) & (x>x0) );
lower_B = find( (y < y0) & (x>x0) );

u1_A = u(upper_A);
u2_A = u(lower_A);
v1_A = v(upper_A);
v2_A = v(lower_A);
w1_A = w(upper_A);
w2_A = w(lower_A);

u1_B = u(upper_B);
u2_B = u(lower_B);
v1_B = v(upper_B);
v2_B = v(lower_B);
w1_B = w(upper_B);
w2_B = w(lower_B);

z1_A = z(upper_A);
z2_A = z(lower_A);
z1_B = z(upper_B);
z2_B = z(lower_B);

[z1_A,i1_A] = sort(z1_A);
[z2_A,i2_A] = sort(z2_A);
[z1_B,i1_B] = sort(z1_B);
[z2_B,i2_B] = sort(z2_B);

du_A = u1_A(i1_A) - u2_A(i2_A);
dv_A = v1_A(i1_A) - v2_A(i2_A);
dw_A = w1_A(i1_A) - w2_A(i2_A);

du_B = u1_B(i1_B) - u2_B(i2_B);
dv_B = v1_B(i1_B) - v2_B(i2_B);
dw_B = w1_B(i1_B) - w2_B(i2_B);

Wrift=10;
a = 4*Wrift;
b = 8*Wrift;

rho = 916; g_const = 9.8;
rhow=1024;
G = 3.6e9;
nu = 0.33;
E = 2*G*(1+nu);
kappa = 3-4*nu;

ra = sqrt( (Wrift/2)^2 + a^2  );
rb = sqrt( (Wrift/2)^2 + b^2  );

KI_A = sqrt(2*pi/ra) * G/(kappa+1) * dv_A;
KII_A = sqrt(2*pi/ra) * G/(kappa+1) * du_A;
KIII_A = sqrt(2*pi/ra) * G * dw_A;

KI_B = sqrt(2*pi/rb) * G/(kappa+1) * dv_B;
KII_B = sqrt(2*pi/rb) * G/(kappa+1) * du_B;
KIII_B = sqrt(2*pi/rb) * G * dw_B;

z=z1_A;
KI_B = interp1(z1_B,KI_B,z);
KII_B = interp1(z1_B,KII_B,z);
KIII_B = interp1(z1_B,KIII_B,z);

KI   = KI_A   + ra/(rb-ra) * (KI_B   - KI_A);
KII  = KII_A  + ra/(rb-ra) * (KII_B  - KII_A);
KIII = KIII_A + ra/(rb-ra) * (KIII_B - KIII_A);

r = rho/rhow; 
phi = 3*r - 2*r^2 - 1;
D0 = E * H^3 / 12 / (1-nu^2);
lambda = ( D0 /  ( rhow * g_const))^(1/4);
sigb = phi * rho * g_const * H/2;
sig0 = (1-r) * rho * g_const * H / 2;

Kb = sigb * sqrt(lambda) * sqrt ( (1-nu)*(3+nu) ) / sqrt(2*pi);
Kb0 = sigb * sqrt(lambda);

Kb_3d = -Kb*((z-(H/2))/(H/2));

Kh =0;
K0 = sig0 * sqrt(pi * L);
Kh_FEM = mean(KI);
ChiH = Kh_FEM / K0;
PsiH = mean(KII) / K0;

PsiHa = mean(KII_A)/K0;
PsiHb = mean(KII_B)/K0;

ChiHa = mean(KI_A)/K0;
ChiHb = mean(KI_B)/K0;

disp(' ');
disp(['                    Psi:  ' num2str(PsiH)]);
disp(['                    Chi:  ' num2str(ChiH)]);
disp(['                     KI:  ' num2str(mean(KI),4)]);
disp(['                    KII:  ' num2str(mean(KII),4)]);
disp(['                   KIII:  ' num2str(mean(KIII),4)]);
disp(['                 lambda:  ' num2str(lambda,4)]);
disp(['                     Kb:  ' num2str(Kb,4)]);
%disp(['                     r2:  ' num2str(list(i)/1e4,4)]);

%figure
%plot(KI,z)
%hold on
%plot(Kb_3d - 1.3e6,z);
%figure(2); plot(KIII,z)

