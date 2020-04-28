% Make figure comparing marginal rifts.
figure(1); clf;
clear; dat_StickOut_Margin; Plot_MixedMode; ylim(100*[-1 1]);
	for i=[1,2,4]; subplot(2,2,i); ylim([-75 75]); end
clear; dat_FixedSides_Margin; Plot_MixedMode; ylim(100*[-1 1]);
	for i=[1,2,4]; subplot(2,2,i); ylim([-75 75]); end
clear; dat_WeakMargins; Plot_MixedMode; ylim(100*[-1 1]);
	for i=[1,2,4]; subplot(2,2,i); ylim([-75 75]); end

