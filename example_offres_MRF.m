% example_offres_MRF.m
% fast offresonance simulation for MRF, with pSSFP scheme
% author: xiang gao
% xiang.gao@uniklinik-freiburg.de

clear;clc;clear global;
global FOV;
global Mat;
global TDImage;     %Output 3D spatial modulation function or 1D signal intensity
TDImage=false;

%==========================================================================
%% simulation
%==========================================================================
NPulse=100;
offres=200;%Hz
FOV = 128; %mm
Mat = floor(offres*pi);
TR=10;%ms

% % pSSFP % Jakob Asslaender, MRM, 2017
disp('Simu MRF type: MRF from Jakob Asslaender; pSSFP')
rng(42);
FA=zeros(4,250);
FA(1,:)=10+sin(2*pi/500*[1:250])*50+normrnd(0,5,[1,250]);FA(2,:)=FA(1,:)./2;
FA(3,:)=10+sin(2*pi/500*[1:250])*50+normrnd(0,5,[1,250]);FA(4,:)=FA(3,:)./2;
theta0=[0,col(FA')']./2;theta1=[col(FA')',0]./2;alpha=theta0+theta1;
TRssfp=TR;TEa=zeros(1,NPulse);TRa=TEa;
theta_ratio=sin(theta0./180*pi)./sin(theta1./180*pi);
TEa(1)=0;
for i=1:NPulse
    if theta_ratio(1+i)>1
        TRa(i)=TEa(i)+TRssfp/2/theta_ratio(1+i);
    else
        TRa(i)=TEa(i)+TRssfp/2;
    end
    TEa(i+1)=(TRa(i)-TEa(i)).*theta_ratio(i+1);
end
TRa(i+1)=TRssfp;
TEa=TEa(1:NPulse);
TRa=TRa(1:NPulse);
Rf=alpha(1:NPulse);
%==========================================================================
%% simulation preparation
%==========================================================================
Tissue{1}.T1=4500;Tissue{1}.T2=2200; %csf
Tissue{2}.T1=1084;Tissue{2}.T2=69;   %white matter
Tissue{3}.T1=1820;Tissue{3}.T2=99;   %gray matter
Tissue{4}.T1=371; Tissue{4}.T2=133;  %fat
type=2;
b_autoKG=false;

%% offset freq gradient
gamma = 42.5756; %MHz/T
dG_freqoff=offres./gamma./FOV;%mT/m
Gy = 0;  Gz = 0;
[Gx] = deal(dG_freqoff);

sequence.npulse=NPulse;
sequence.time=TRa(1:NPulse);        % ms
sequence.echo=TEa(1:NPulse); 
sequence.T1=Tissue{type}.T1;
sequence.T2=Tissue{type}.T2;
sequence.gradx=Gx.*ones(NPulse,1);
sequence.grady=Gy.*ones(NPulse,1);
sequence.gradz=Gz.*ones(NPulse,1);
sequence.velocityx=zeros(NPulse,1);     % transport velocity in mm/s
sequence.velocityy=zeros(NPulse,1);     % transport velocity in mm/s
sequence.velocityz=zeros(NPulse,1);     % transport velocity in mm/s
sequence.angle=Rf;
if rem(NPulse,2)
    sequence.axes=[180;repmat([0;180],floor(NPulse/2),1)];
else
    sequence.axes=repmat([0;180],NPulse/2,1);
end
% sequence.diff=0;                       % mcm^2/ms  %2 for water
sequence.diffTensor=[0;0;0;0;0;0;0;0;0];  % mcm^2/ms
sequence.Omega=0;                      % Hz, frequency offset
sequence.ktolerance=1.e-8;             % computational parameter:
sequence.ktolerance_phys=2.e-6;       % grid merging size 1e-3 = 1*(1/mm)
sequence.nOutput = 1;
matlabseq = ['matlabseq_',int2str(sequence.nOutput),'.txt'];
matlab2c(sequence, matlabseq);
string = ['3D_SR-PG_tensor.exe ',matlabseq];
%==========================================================================
%% simulation execution
%==========================================================================
if b_autoKG    % simulation execution (auto adjusting k-grid)
    error_tolerance = 0.01; % signal changes below tolerance =1%
    ktolerance_phys_step = 0.2; % every step, ktol set as 20% of previous one
    ktolerance_phys_hard_quit = 2e-8;
    
    error_calc = 1;sig_last = 0;cnt=1;
    while (error_calc>error_tolerance)
        disp(['searching ktolerance_phys= ',mat2str(sequence.ktolerance_phys)]);
        tic;system(string);t=toc;
        logData=KmapSignal2(['output',int2str(sequence.nOutput),'.txt']);%faster, but need preallocate memory(!).
        logData=IMG_DFT(logData);
        logData.t=t;
        sig_cur=sum(logData.sig3d(:,end));
        error_calc = abs(sig_cur-sig_last)./sig_last;
        disp(['ERROR_calc= ',mat2str(error_calc)]);
        sig_last = sig_cur;
        
        logData.sequence = sequence;
        logData2{cnt}=logData;cnt=cnt+1;
        
        sequence.ktolerance_phys=sequence.ktolerance_phys*ktolerance_phys_step;
        if sequence.ktolerance_phys < ktolerance_phys_hard_quit
            break;
        end
        matlabseq = ['matlabseq_',int2str(sequence.nOutput),'.txt'];
        matlab2c(sequence, matlabseq);
    end
    disp(['auto ktolerance_phys= ',mat2str(sequence.ktolerance_phys)]);
    logData=IMG_DFT(logData);
    plot_performance(logData2);
else
    tic;system(string);t=toc;disp(t);
    logData=KmapSignal2(['output',int2str(sequence.nOutput),'.txt']);%faster, but need preallocate memory(!).
    logData.sequence = sequence;
    logData=IMG_DFT(logData);
    center=size(logData.sig3d,1)/2;
    figure;
    plot(logData.sig3d(:,end));hold on
    xticks([center-250,center-125,center,center+125,center+250]);xlim([0,Mat])
    set(gca,'xticklabels',{'-500','-250','0','250','500'});ylabel('|M_\perp| / M_0')
end
return

%%
% load('logData2.mat'); 
% plot_performance(logData2);
% cla(subplot(2,2,4))
% % pSSFP  Spectrum
% center=size(logData2{1}.sig3d,1)/2;
% subplot(2,2,4)
% plot(logData2{end}.sig3d(:,end),'linewidth',1.5);hold on;
% load('logData2_diff.mat');
% plot(logData2{end}.sig3d(:,end),'linewidth',1.5);hold off;
% legend('D=0','D=0.2 cm^2/s')
% xticks([center-250,center-125,center,center+125,center+250])
% set(gca,'xticklabels',{'-500','-250','0','250','500'});ylabel('|M_\perp| / M_0')
% xlabel('\omega [rad/s]'),ylim([0,0.25])
% title(['pSSFP  Spectrum']);
% set(gca,'FontSize',11);




