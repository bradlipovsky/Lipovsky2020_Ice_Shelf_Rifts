clear;

load SimuK-IceTongue

L = 2500;
G = 3.6e9;
nu = 0.33;
E = 2*G*(1+nu);
E1 = E/(1-nu^2);
E2 = E/(1+nu);
r = 916/1024;
phi = 3*r - 2*r^2 - 1;
rho = 916; g_const = 9.8;
Kc = 100e3;

A_List = unique(alpha);
if ~exist('delta')
    delta = ones(size(alpha));
end
D_List = unique(delta);
H=200;

D0 = E * H^3 / 12 / (1-nu^2);
lambda = ( D0 / (4 * rho * g_const))^(1/4);
sigb = phi * rho * g_const * H/2;
Kb = sigb * sqrt(lambda);
sig0 = (1-r) * rho * g_const * H / 2;

ThisD = D_List(1);
SubSubset = find( (delta==ThisD)  );
[~,srt]=sort(alpha(SubSubset));
TheseAlphas = alpha(SubSubset(srt));
TheseChis=chi(SubSubset(srt));
ThesePsis=psi(SubSubset(srt));
	
KI = TheseChis * sig0 * sqrt(pi*L) - Kb;
KII= ThesePsis * sig0 * sqrt(pi*L);

th = -2*atan2(8*KII, (-2*KI + 2*sqrt(KI.^2 + 8*KII.^2))  );
MM = cos(th/2).*(KI/Kc .* cos(th/2).^2 - 1.5*KII/Kc.*sin(th));

subplot(221); plot(TheseAlphas, KI/Kc,'linewidth',2); hold on;
xlabel('Rift position (distance from ice front), W/Ly');
title('A. Mode I Stress Intensity Factor, K_I / K_{Ic}');
axis tight;

subplot(222); plot(TheseAlphas, KII/Kc,'linewidth',2); hold on;
xlabel('Rift position (distance from ice front), W/Ly');
title('B. Mode II Stress Intensity Factor, K_{II} / K_{Ic}');

th = th*180/pi;
th(th>120) = th(th>120) - 360;
subplot(223); plot(TheseAlphas, th,'linewidth',2); hold on;
xlabel('Rift position (distance from ice front), W/Ly');
ylabel('Degrees');
title('C. Propagation Angle');
axis tight;

subplot(224); plot(TheseAlphas, MM,'linewidth',2); hold on;
xlabel('Rift position (distance from ice front), W/Ly');
title('D. Optimally Oriented Mode I SIF, K_{I}^{Op} / K_{Ic}');
axis tight;
