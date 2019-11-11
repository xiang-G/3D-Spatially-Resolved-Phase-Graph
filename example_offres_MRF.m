% example_offres_MRF.m
% fast offresonance simulation for MRF, bSSFP vs. original MRF vs. pSSFP 
% author: xiang gao
% xiang.gao@uniklinik-freiburg.de

clear;clc;clear global;
global FOV;
global Mat;
global Kmax2Pi;     %Skope of 3D-K space
global TDImage;     %Output 3D image or 1D signal intensity
TDImage=false;
Rep = 4;            %the number of tissues in simulation (def=4) 
                    %csf,white matter,gray matter,fat                    
%==========================================================================
%% simulation
%==========================================================================
NPulse=51;
gamma = 42.5756; %MHz/T
FOV = 128; %mm
Mat = floor(200*pi);
Kmax2Pi=0;
offres=200;%Hz
dG=offres./gamma./FOV;%mT/m
MRF=cell(3,1);
for Simu_MRF_type=1:3   %1:bSSFP; 2:original_MRF; 3:pSSFP    
    MRF{Simu_MRF_type} = Simu_MRF(Simu_MRF_type,dG,NPulse,Rep);
end
%==========================================================================
%% output
%==========================================================================
MRF_SSFP=MRF{1};
MRF_DMa=MRF{2};
MRF_Jacob=MRF{3};
center=size(MRF_SSFP{1}.sig3d,1)/2;
piece=':';

figure;
subplot(5,1,1)
plot(MRF_SSFP{1}.sig3d(piece,end),'c');hold on;plot(MRF_SSFP{2}.sig3d(piece,end),'b');hold on;
plot(MRF_SSFP{3}.sig3d(piece,end),'r');hold on;plot(MRF_SSFP{4}.sig3d(piece,end),'k');
xticks([center-250,center-125,center,center+125,center+250]),ylim([0,0.5])
set(gca,'xticklabels','');legend('CSF','WM','GM','Fat'),ylabel('|M_\perp| / M_0')
subplot(5,1,2)
plot(MRF_DMa{1}.sig3d(piece,end-1),'c');hold on;plot(MRF_DMa{2}.sig3d(piece,end-1),'b');hold on;
plot(MRF_DMa{3}.sig3d(piece,end-1),'r');hold on;plot(MRF_DMa{4}.sig3d(piece,end-1),'k');
xticks([center-250,center-125,center,center+125,center+250]),ylim([0,0.5])
set(gca,'xticklabels','');ylabel('|M_\perp| / M_0')
subplot(5,1,3)
plot(MRF_DMa{1}.sig3d(piece,end),'c');hold on;plot(MRF_DMa{2}.sig3d(piece,end),'b');hold on;
plot(MRF_DMa{3}.sig3d(piece,end),'r');hold on;plot(MRF_DMa{4}.sig3d(piece,end),'k');
xticks([center-250,center-125,center,center+125,center+250]),ylim([0,0.5])
set(gca,'xticklabels','');ylabel('|M_\perp| / M_0')
subplot(5,1,4)
plot(MRF_Jacob{1}.sig3d(piece,end-1),'c');hold on;plot(MRF_Jacob{2}.sig3d(piece,end-1),'b');hold on;
plot(MRF_Jacob{3}.sig3d(piece,end-1),'r');hold on;plot(MRF_Jacob{4}.sig3d(piece,end-1),'k');
xticks([center-250,center-125,center,center+125,center+250]),ylim([0,0.5])
set(gca,'xticklabels','');ylabel('|M_\perp| / M_0')
subplot(5,1,5)
plot(MRF_Jacob{1}.sig3d(piece,end),'c');hold on;plot(MRF_Jacob{2}.sig3d(piece,end),'b');hold on;
plot(MRF_Jacob{3}.sig3d(piece,end),'r');hold on;plot(MRF_Jacob{4}.sig3d(piece,end),'k');
xticks([center-250,center-125,center,center+125,center+250]),ylim([0,0.5])
set(gca,'xticklabels',{'-500','-250','0','250','500'});ylabel('|M_\perp| / M_0')
xlabel('\omega [rad/s]')

%enlargement
piece=[center-23:center+23];
figure;set(gcf,'outerposition',[500 500 120 150]);
plot(MRF_DMa{1}.sig3d(piece,end-1),'c','linewidth',2);hold on;
plot(MRF_DMa{2}.sig3d(piece,end-1),'b','linewidth',2);hold on;
plot(MRF_DMa{3}.sig3d(piece,end-1),'r','linewidth',2);hold on;
plot(MRF_DMa{4}.sig3d(piece,end-1),'k','linewidth',2);
xticks([center-250,center-125,center,center+125,center+250]),ylim([0,0.3])
set(gca,'xticklabels','','yticklabels','');
figure;set(gcf,'outerposition',[500 500 120 150]);
plot(MRF_DMa{1}.sig3d(piece,end),'c','linewidth',2);hold on;
plot(MRF_DMa{2}.sig3d(piece,end),'b','linewidth',2);hold on;
plot(MRF_DMa{3}.sig3d(piece,end),'r','linewidth',2);hold on;
plot(MRF_DMa{4}.sig3d(piece,end),'k','linewidth',2);
xticks([center-250,center-125,center,center+125,center+250]),ylim([0,0.3])
set(gca,'xticklabels','','yticklabels','');
figure;set(gcf,'outerposition',[500 500 120 150]);
plot(MRF_Jacob{1}.sig3d(piece,end-1),'c','linewidth',2);hold on;
plot(MRF_Jacob{2}.sig3d(piece,end-1),'b','linewidth',2);hold on;
plot(MRF_Jacob{3}.sig3d(piece,end-1),'r','linewidth',2);hold on;
plot(MRF_Jacob{4}.sig3d(piece,end-1),'k','linewidth',2);
xticks([center-250,center-125,center,center+125,center+250]),ylim([0,0.3])
set(gca,'xticklabels','','yticklabels','');
figure;set(gcf,'outerposition',[500 500 120 150]);
plot(MRF_Jacob{1}.sig3d(piece,end),'c','linewidth',2);hold on;
plot(MRF_Jacob{2}.sig3d(piece,end),'b','linewidth',2);hold on;
plot(MRF_Jacob{3}.sig3d(piece,end),'r','linewidth',2);hold on;
plot(MRF_Jacob{4}.sig3d(piece,end),'k','linewidth',2);
xticks([center-250,center-125,center,center+125,center+250]),ylim([0,0.3])
set(gca,'xticklabels','','yticklabels','');

%
figure; %RF,TR&TE protocol settings
subplot(2,2,1)
plot(MRF_DMa{1}.sequence.angle./180,'b-'); hold on;ylim([0,0.4]);ylabel('\alpha / \pi')
title('MRF'),xlim([1,length(MRF_DMa{1}.sequence.angle)]);
subplot(2,2,2)
plot(MRF_Jacob{1}.sequence.angle./180,'b-');hold on;ylim([0,0.4]);title('pSSFP')
xlim([1,length(MRF_DMa{1}.sequence.angle)]);
subplot(2,2,3)
plot(MRF_DMa{1}.sequence.echo,'r-'); hold on;ylim([0,12]);
plot(MRF_DMa{1}.sequence.time,'b-'); hold on;ylim([0,12]);
xlim([1,length(MRF_DMa{1}.sequence.angle)]);ylabel('t [ms]');xlabel('t / TR');
subplot(2,2,4)
plot(MRF_Jacob{1}.sequence.echo,'r-'); hold on;ylim([0,12]);
plot(MRF_Jacob{1}.sequence.time,'b-'); hold off;ylim([0,12]);
xlim([1,length(MRF_DMa{1}.sequence.angle)]);legend('TR','TE');xlabel('t / TR');



