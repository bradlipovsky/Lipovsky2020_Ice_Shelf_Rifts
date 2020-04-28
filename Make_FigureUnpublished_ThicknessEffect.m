clear;

H_List = [25,50,100,200,400,800,1600];
KI = zeros(size(H_List));
KIII = zeros(size(H_List));
Kb0 = zeros(size(H_List));

for i = 1:numel(H_List)
	H = H_List(i);
	commnd = ['dcm_StickOut_3d_H' num2str(H) 'm'];
	eval(commnd)
	[~,~,K1,K3,Kb0(i)] = Run3d_DCM_NoPlots (data,H,97535,40000);
	pI = polyfit(linspace(0,H,numel(K1))',K1,1);
	pIII = polyfit(linspace(0,H,numel(K3))',K3,1);
	%KI(i) = -pI(1) * H;
	%KIII(i) = -pIII(1) * H;
	KI(i)   = K1(1) - mean(K1);
	KIII(i) = K3(1) - mean(K3);
end

p1 = polyfit( log10(H_List) , log10(KI) , 1);
p3 = polyfit( log10(H_List) , log10(KIII) , 1);

figure(1);
subplot(1,2,1); loglog(H_List,KI,'o'); hold on
subplot(1,2,2); loglog(H_List,KIII,'o'); hold on;

subplot(121); loglog(H_List, 10^(p1(2)).*H_List.^p1(1) );
title(['KIb ~ H \^ ' num2str(p1(1))]);
subplot(122); loglog(H_List, 10^(p3(2)).*H_List.^p3(1) );
title(['KIIIb ~ H \^ ' num2str(p3(1))]);
