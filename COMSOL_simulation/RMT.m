%% Eigenvalues Distribution
load('eigenvalues.mat')
[modes,~] = size(eig8);
tran_uniform = linspace(0,1,modes);
bins = 60;
set(groot, 'DefaultAxesFontSize',20);
figure(1)
plot(1:modes,eig8,'r-',...
    1:modes,eig6,'g-',...
    1:modes,eig4,'b-')
legend('8 $\mu$m thick slab','6 $\mu$m thick slab','4 $\mu$m thick slab')
title('Transmission Eigenvalues')
xlabel('Eigenmode Index, $i$')
ylabel('$T_i$')
ylim([0,1.1])
xlim([0,243])
set(groot, 'DefaultAxesFontSize',24);

figure(2)
plot(tran_uniform,distribution8,'r-',...
    tran_uniform,distribution6,'g-',...
    tran_uniform,distribution4,'b-')
% set(gca, 'Yscale',  'log')
legend('8 $\mu$m thick slab','6 $\mu$m thick slab','4 $\mu$m thick slab')
ylim([0 distribution8(end-1)])
set(gca,'YTickLabel',[]);
title('Theoretical Bimodal Distribution')

figure(3)
histogram(eig8,bins,'FaceColor','r')
hold on
plot(tran_uniform,distribution8,'r-')
hold off
xlim([0,1])
ylim([0,50])
title('8 $\mu$m thick slab')
legend('Simulated Eigenvalues','Theoretical Distribution')

figure(4)
histogram(eig6,bins,'FaceColor','g')
hold on
plot(tran_uniform,distribution6,'g-')
hold off
xlim([0,1])
ylim([0,50])
title('6 $\mu$m thick slab')
legend('Simulated Eigenvalues','Theoretical Distribution')

figure(5)
histogram(eig4,bins,'FaceColor','b')
hold on
plot(tran_uniform,distribution4,'b-')
hold off
xlim([0,1])
ylim([0,50])
title('4 $\mu$m thick slab')
legend('Simulated Eigenvalues','Theoretical Distribution')

scaling = 1;
figure(6)
histogram(eig8,bins*scaling,'FaceColor','r')
hold on
plot(tran_uniform,distribution8/scaling,'r-')
histogram(eig6,bins*scaling,'FaceColor','b')
plot(tran_uniform,distribution6/scaling,'b-')
histogram(eig4,bins*scaling,'FaceColor','g')
plot(tran_uniform,distribution4/scaling,'g-')
hold off
legend('8 $\mu$m Eigenvalues', 'Dist for 8 $\mu$m','6 $\mu$m Eigenvalues','Dist for 6 $\mu$m','4 $\mu$m Eigenvalues','Dist for 4 $\mu$m')
xlim([0,1])
ylim([0,50])
% set(gca, 'Yscale',  'log')

maxim8= mean(eig8);
set(groot, 'DefaultAxesFontSize',20);
figure(7)
histogram(eig8/mean(eig8),bins*scaling,'FaceColor','r')
hold on
histogram(eig6/mean(eig6),bins*scaling,'FaceColor','b')
histogram(eig4/mean(eig4),bins*scaling,'FaceColor','g')
hold off
leg1 = legend(['8 $\mu$m thick, var$_{\tilde{T}_{8}}$ = ' num2str(var(eig8/mean(eig8)))],...
    [ '6 $\mu$m thick, var$_{\tilde{T}_{6}}$ = ' num2str(var(eig6/mean(eig6)))],...
    ['4 $\mu$m thick, var$_{\tilde{T}_{4}}$ = ' num2str(var(eig4/mean(eig4)))]);
xlim([0,5])
set(gca, 'XTick', 0:1:5)
ylim([0,50])
ylabel('no. of modes')
xlabel('$\tilde{T}$, $T$ normalized by its mean')
set(groot, 'DefaultAxesFontSize',24);

% [N,edges] = histcounts(eig8)
% summed = 0;
% combos = combntns(1:modes,2);
% [num_combo,~] =  size(combos);
% eigchosen = eig4;
% for i = 1:num_combo
%     temp = eigchosen(combos(i,1))-eigchosen(combos(i,2));
%     summed = summed + temp;
% end
% summed/num_combo/sum(eigchosen)

figure(8)
set(groot, 'DefaultAxesFontSize',20);
plot(1:modes,eig8/mean(eig8),'r-')
hold on
plot(1:modes,eig6/mean(eig6),'b-')
plot(1:modes,eig4/mean(eig4),'g-')
plot(1:modes,ones(243,1),'k-')
hold off
legend(['8 $\mu$m thick slab, var$_{\tilde{T}_{8}}$ = ' num2str(var(eig8/mean(eig8)))],...
    ['6 $\mu$m thick slab, var$_{\tilde{T}_{6}}$ = ' num2str(var(eig6/mean(eig6)))],...
    ['4 $\mu$m thick slab, var$_{\tilde{T}_{4}}$ = ' num2str(var(eig4/mean(eig4)))],...
    ['free space, var$_{\tilde{T}}$ = 0'])
xlim([0,243])
ylim([0,5])
title('Normalized Transmission Eigenvalues')
xlabel('Eigenmode Index, $i$')
ylabel('$\tilde{T}_i$')
set(groot, 'DefaultAxesFontSize',24);
% set(gca, 'Yscale',  'log')

