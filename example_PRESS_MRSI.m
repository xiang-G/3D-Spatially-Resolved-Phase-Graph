% example_PRESS_MRSI.m
% signal evolution in 3d k-space domain & image domain
% author: xiang gao
% xiang.gao@uniklinik-freiburg.de

clear;clc;clear global;
global FOV;
global Mat;
global TDImage;     %Output 3D spatial modulation function or 1D signal intensity
TDImage=true;

%==========================================================================
%% sequence preparation
%==========================================================================
Nphaccyc=1; %num of phase cycle;
gamma = 42.5756; %MHz/T
FOV = 48; %mm
Mat = 16;
Position = [];%phantomP;
RF_BandWidth = 4000; %Hz
Gs = RF_BandWidth/gamma/FOV;%mT/m
Gs_Duration = 4; %ms
KGs = Gs*Gs_Duration*gamma;
TE1hf = 7; %ms
TE2hf = 8; %ms
KSpoil1 = 20000./FOV; %1/m
KSpoil2 = 20000./FOV; %1/m
Gx = [-KSpoil1/TE1hf;         -KSpoil1/TE1hf;         -KSpoil2/TE2hf;         -KSpoil2/TE2hf;0;0]./gamma;
Gy = [(-KSpoil1-KGs/2)/TE1hf; (-KSpoil1-KGs/2)/TE1hf; -KSpoil2/TE2hf;         -KSpoil2/TE2hf;0;0]./gamma;
Gz = [-KSpoil1/TE1hf;         -KSpoil1/TE1hf;         (-KSpoil2-KGs/2)/TE2hf; (-KSpoil2-KGs/2)/TE2hf;0;0]./gamma;

% freqoff=[0,0,0];
freqoff=0.1*[1,0.5,1.5]; %hz/cm
dG = freqoff./10./gamma;
Gx = Gx+dG(:,1);  Gy = Gy+dG(:,2);  Gz = Gz+dG(:,3);

%==========================================================================
%% simulation preparation
%==========================================================================
sequence=cell(Nphaccyc,1);
for i=1:Nphaccyc
    NPulse=6;
    sequence{i}.npulse=NPulse;
    sequence{i}.time=[TE1hf;TE1hf;TE2hf;TE2hf;250;250];% ms
    sequence{i}.echo=[TE1hf;TE1hf;TE2hf;TE2hf;250;250];% ms
    sequence{i}.gradx=Gx;
    sequence{i}.grady=Gy;
    sequence{i}.gradz=Gz;
    sequence{i}.velocityx=zeros(1,NPulse);     % transport velocity in mm/s
    sequence{i}.velocityy=zeros(1,NPulse);     % transport velocity in mm/s
    sequence{i}.velocityz=zeros(1,NPulse);     % transport velocity in mm/s
    sequence{i}.angle=[90/2;180/2;0;180/2;0;0];
    sequence{i}.axes=[0,0*ones(1,NPulse-1)];  % rotation axis.
    [phs1,phs2,phs3]=press_phacyc(i,Nphaccyc);
    sequence{i}.axes=[phs1,phs2,0,phs3,0,0];  % rotation axis.
    % sequence.diff=0;                        % mcm^2/ms  %2 for water
    sequence{i}.diffTensor=[0;0;0;0;0;0;0;0;0];  % mcm^2/ms  %2 for water
    sequence{i}.T1=0;                          % ms, 0 equivalent to infinity
    sequence{i}.T2=0;                          % ms, 0 equivalent to infinity
    sequence{i}.Omega=0;                       % Hz, frequency offset
    sequence{i}.ktolerance=1.e-8;              % computational parameter:
    sequence{i}.ktolerance_phys=1.e-6;         % 3d k-vectors gridding size 1e-3 = 1 rad/mm
    sequence{i}.nOutput=i;
    matlabseq = ['matlabseq_',int2str(sequence{i}.nOutput),'.txt'];
    matlab2c(sequence{i}, matlabseq);
    string{i} = ['3D_SR-PG_tensor.exe ',matlabseq];
end
%==========================================================================
%% simulation execution
%==========================================================================
sig_phacyc=0;
for i=1:Nphaccyc
    tic
    system(string{i});
    toc
    logData=KmapSignal2(['output',int2str(sequence{i}.nOutput),'.txt']);%faster, but need preallocate memory(!).
    logData.sequence = sequence{i};
    sig_phacyc=logData.signal+sig_phacyc;
end
logData.signal=sig_phacyc/Nphaccyc;
logData=clean_signal(logData,sequence{1}.ktolerance);%clean residual signal(noise) from sig_phacyc/Nphaccyc
logData=IMG_DFT(logData,Position);
%==========================================================================
%% simulation execution (auto adjusting k-grid)
%==========================================================================
% error_tolerance = 0.05; % MSE error changes below tolerance =5%
% ktolerance_phys_step = 0.2; % every step, ktol set as 20% of previous one
% ktolerance_phys_hard_quit = sequence{1}.ktolerance;
%
% error_calc = 1;
% sig_last = 0;
% while error_calc>error_tolerance
%     sig_phacyc=0;
%     for i=1:Nphaccyc
%         system(string{i});
%         logData=KmapSignal2(['output',int2str(sequence{i}.nOutput),'.txt']);%faster, but need preallocate memory(!).
%         sig_phacyc=logData.signal+sig_phacyc;
%     end
%     logData.signal=sig_phacyc/Nphaccyc;
%     logData=clean_signal(logData,sequence{1}.ktolerance);%suppose signal noise threshold = kcalc precision.
%     sig_cur = sum(abs(logData.signal(:,1)+1i*logData.signal(:,2)));
%     error_calc = abs(sig_cur-sig_last)./sig_last;
%     disp(['MSE_calc= ',mat2str(error_calc)]);
%     disp(['searching ktolerance_phys= ',mat2str(sequence{1}.ktolerance_phys)]);
%     if sequence{1}.ktolerance_phys < ktolerance_phys_hard_quit
%         break;
%     end
%     sig_last = sig_cur;
%
%     for i=1:Nphaccyc
%         sequence{i}.ktolerance_phys=sequence{i}.ktolerance_phys*ktolerance_phys_step;
%         matlabseq = ['matlabseq_',int2str(sequence{i}.nOutput),'.txt'];
%         matlab2c(sequence{i}, matlabseq);
%     end
% end
% disp(['auto ktolerance_phys= ',mat2str(sequence{1}.ktolerance_phys)]);
% logData.sequence = sequence{1};
% logData=IMG_DFT(logData,Position);
%==========================================================================
%% OUTPUT
%==========================================================================
%% plot 3D kspace graph
% we plot from second echo
Nimg=NPulse-2;
nth = [2,4,5,6];

figure; set(gcf,'outerposition',[0 0 1000 800]);
for i=1:Nimg
    subplot(2,Nimg/2,i);
    plot_kmap3d(logData,nth(i),2000,i);
end

%% plot projected image
figure; set(gcf,'outerposition',[0 0 500 400]);colormap gray;
for i=1:Nimg
    subplot(3,Nimg,i)
    imagesc(single(logData.sig3d(:,:,nth(i),1)));caxis([0,0.4]);
    title(['Time ',int2str(i)]);
    if i==1;ylabel('XY','FontWeight','bold');end
    if i==Nimg;colorbar('Manual', 'position', [0.93 0.1 0.02 0.81]);end
    xticks(1:(Mat-1)/2:Mat);xticklabels({int2str(-FOV/2),'0',int2str(FOV/2)});
    yticks(1:(Mat-1)/2:Mat);yticklabels({int2str(-FOV/2),'0',int2str(FOV/2)});
end
for i=1:Nimg
    subplot(3,Nimg,Nimg+i)
    imagesc(single(logData.sig3d(:,:,nth(i),2)));caxis([0,0.4]);
    if i==1;ylabel('YZ','FontWeight','bold');end
    xticks(1:(Mat-1)/2:Mat);xticklabels({int2str(-FOV/2),'0',int2str(FOV/2)});
    yticks(1:(Mat-1)/2:Mat);yticklabels({int2str(-FOV/2),'0',int2str(FOV/2)});
end
for i=1:Nimg
    subplot(3,Nimg,2*Nimg+i)
    imagesc(single(logData.sig3d(:,:,nth(i),3)));caxis([0,0.4]);
    if i==1;ylabel('XZ','FontWeight','bold');end
    xticks(1:(Mat-1)/2:Mat);xticklabels({int2str(-FOV/2),'0',int2str(FOV/2)});
    yticks(1:(Mat-1)/2:Mat);yticklabels({int2str(-FOV/2),'0',int2str(FOV/2)});
end

return
