% function plot_performance()
% plot the calculation burden, residual error, computation time and 1-D plot
%
function plot_performance(logData2)

Nsim=length(logData2);
lw=1.5;
fz=11;
figure;
set(gcf,'outerposition',[10 10 1000 800]);

%% TE and TR
subplot(2,2,1)
plot(logData2{1}.sequence.echo,'linewidth',lw); hold on; 
plot(logData2{1}.sequence.time,'linewidth',lw); hold on; ylim([0,12]);
xlim([1,length(logData2{1}.sequence.angle)]);
ylabel('t [ms]');
yyaxis right
plot(logData2{1}.sequence.angle./180,'linewidth',lw); hold off; 
ylabel('\alpha/\pi');

legend('TE','TR','FA','Location','southeast');xlabel('t / TR');
xlabel('Num of RF pulses');
title('Protocol parameters');set(gca,'FontSize',fz);
%% NKvec vs. Necho
subplot(2,2,2);
a=1/2*(3.^(1:logData2{1}.sequence.npulse)-1);
for i=1:Nsim
semilogy(1.5*logData2{i}.NKperEcho,'linewidth',lw);hold on; %adding kz for each echo, which is skipped for the output file.
legendInfo{i} = ['Kg ',mat2str(logData2{i}.sequence.ktolerance_phys*1e6),' rad/m']; 
end
semilogy(a(1:20),'k--','linewidth',lw);hold off;    
legend(legendInfo,'Location','southeast');
xlabel('Num of RF pulses');
ylabel('Num of K vectors');
title('Burden');set(gca,'FontSize',fz);
%% Computation Time and Erroneous vs. K-grid radius
subplot(2,2,3)
for i=Nsim:-1:1
    xlab(i)=logData2{i}.sequence.ktolerance_phys*1e6;
    YY(i)= logData2{i}.t; 
end
plot(YY,'-*','linewidth',lw); hold on;
ylabel('Computation Time [s]');

error=zeros(1,Nsim);error(1)=nan;
for i=Nsim:-1:2
    error(i) = abs(sum(logData2{i}.sig3d(:,end))-sum(logData2{i-1}.sig3d(:,end)))./sum(logData2{i-1}.sig3d(:,end))*100;
end
yyaxis right
plot(error,'-*','linewidth',lw); hold off;
ylabel('Successive iterate changes [%]');

xticks(1:length(YY));xticklabels(xlab);
xlabel('K-grid radius (rad/m)');
title('Efficiency');set(gca,'FontSize',fz);

%% pSSFP  Spectrum
center=size(logData2{1}.sig3d,1)/2;
subplot(2,2,4)
for i=1:Nsim
    plot(logData2{i}.sig3d(:,end),'linewidth',lw);hold on;
end
xticks([center-250,center-125,center,center+125,center+250])
set(gca,'xticklabels',{'-500','-250','0','250','500'});ylabel('|M_\perp| / M_0')
xlabel('\omega [rad/s]');xlim([1,628*2]);%ylim([0,0.3])
title(['pSSFP  Spectrum']);
set(gca,'FontSize',fz);
end
