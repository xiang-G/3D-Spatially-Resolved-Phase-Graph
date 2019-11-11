function logData=Simu_MRF(type,dG,NPulse,Rep)
% Simulate signal spectrum for different tissue types
% In:  type:    MRF_types
%      dG:      static background gradient
%      NPulse:  the number of pulses in simulation
%      Rep:     the number of tissues in simulation (4) 
% Out: logData
%==========================================================================
%% sequence preparation
%==========================================================================
sequence=cell(Rep,1); logData = cell(Rep,1); string=cell(Rep,1);
TR=10;%ms
Gx = dG;  Gy = 0;  Gz = 0; 
%% MRF_type
switch type
    case 1  % 1 constant RF, bSSFP
        disp('Simu MRF type: bSSFP')
        FA=26.3; Rf = [FA/2;FA.*ones(NPulse-1,1)];
        for i=1:Rep
            sequence{i}.time=[TR/2;TR.*ones(NPulse-1,1)];      % ms
            sequence{i}.echo=sequence{1}.time./2;
            sequence{i}.gradx=Gx.*ones(1,NPulse);
        end
    case 2  % 2 VFL %D Ma, Nature, 2013
        disp('Simu MRF type: original MRF from DMa with fixed TE&TR')
        FA=zeros(4,250);
        FA(1,:)=10+sin(2*pi/500.*[1:250])*50+normrnd(0,(5),[1,250]);FA(2,:)=FA(1,:)./2;
        FA(3,:)=10+sin(2*pi/500.*[1:250])*50+normrnd(0,(5),[1,250]);FA(4,:)=FA(3,:)./2;
        save('DMA_FA.mat','FA');
        Rf = [180,FA(1,:),zeros(1,50),FA(2,:),zeros(1,50),FA(3,:),zeros(1,50),FA(4,:)];
        Rf = Rf(1:NPulse);Rf(Rf<0)=0;
        Gxx = Gx.*[0,ones(1,250),zeros(1,50),ones(1,250),zeros(1,50),ones(1,250),zeros(1,50),ones(1,250)];
        for i=1:Rep
            sequence{i}.gradx=Gxx(1:NPulse);
            sequence{i}.time=[0,TR.*ones(1,NPulse-1)];      % ms
            sequence{i}.echo=sequence{1}.time./2;
        end
    case 3  % 3 VFL %Jakob Asslaender, MRM, 2017
        disp('Simu MRF type: MRF from Jakob Asslaender; pSSFP')
        if exist('DMA_FA.mat')
            load DMA_FA.mat
        else
            FA=zeros(4,250);
            FA(1,:)=10+sin(2*pi/500*[1:250])*50+normrnd(0,5,[1,250]);FA(2,:)=FA(1,:)./2;
            FA(3,:)=10+sin(2*pi/500*[1:250])*50+normrnd(0,5,[1,250]);FA(4,:)=FA(3,:)./2;
        end
        theta0=[0,col(FA')']./2;theta1=[col(FA')',0]./2;alpha=theta0+theta1;
        Rf = [180,alpha(1:NPulse-1)];Rf(Rf<0)=0;
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
        for i=1:Rep
            sequence{i}.time=[0,TRa(2:NPulse)];      % ms
            sequence{i}.echo=[0,TEa(2:NPulse)];
            sequence{i}.gradx=[0,Gx.*ones(1,NPulse-1)];
        end
end
%==========================================================================
%% simulation preparation
%==========================================================================
sequence{1}.T1=4500;sequence{1}.T2=2200; %csf
sequence{2}.T1=1084;sequence{2}.T2=69;   %white matter
sequence{3}.T1=1820;sequence{3}.T2=99;   %gray matter
sequence{4}.T1=371; sequence{4}.T2=133;  %fat
for irep=1:Rep
    sequence{irep}.npulse=NPulse;
    sequence{irep}.grady=Gy.*ones(NPulse,1);
    sequence{irep}.gradz=Gz.*ones(NPulse,1);
    sequence{irep}.velocityx=zeros(NPulse,1);     % transport velocity in mm/s
    sequence{irep}.velocityy=zeros(NPulse,1);     % transport velocity in mm/s
    sequence{irep}.velocityz=zeros(NPulse,1);     % transport velocity in mm/s
    sequence{irep}.angle=Rf;
    if rem(NPulse,2)
        sequence{irep}.axes=[180;repmat([0;180],floor(NPulse/2),1)];
    else
        sequence{irep}.axes=repmat([0;180],NPulse/2,1);
    end
    sequence{irep}.diff=0;                       % mcm^2/ms  %2 for water
    sequence{irep}.Omega=0;                      % Hz, frequency offset
    sequence{irep}.ktolerance=1.e-8;             % computational parameter:
    sequence{irep}.ktolerance_phys=1.e-6;        % merging size (cell)  1e-3 = 1*(1/mm)
    sequence{irep}.nOutput = 20+irep;
    matlabseq = ['matlabseq_',int2str(sequence{irep}.nOutput),'.txt'];
    matlab2c(sequence{irep}, matlabseq);

    if ismac
        string{irep} = ['./3D_SR-PG_ios.out ',matlabseq];
    elseif isunix
        string{irep} = ['./3D_SR-PG.out ',matlabseq];
    else
        string{irep} = ['3D_SR-PG.exe ',matlabseq];
    end
end
%% simulation execution && signal calculation
parfor irep=1:Rep
    tic
    system(string{irep});
    toc
    logData{irep}=KmapSignal2(['output',int2str(sequence{irep}.nOutput),'.txt']);%faster, but need preallocate memory(!).
end
for irep=1:Rep
    logData{irep}.sequence = sequence{irep};
    logData{irep}=IMG_DFT(logData{irep});
end

%==========================================================================
%% OUTPUT
%==========================================================================
% Signal Spectrum
% center=size(logData{1}.sig3d,1)/2;
% figure;
% plot(logData{1}.sig3d(:,end-1),'c');hold on;plot(logData{2}.sig3d(:,end-1),'b');hold on;
% plot(logData{3}.sig3d(:,end-1),'r');hold on;plot(logData{4}.sig3d(:,end-1),'k');
% legend('CSF','WM','GM','Fat')
% xticks([center-250,center-125,center,center+125,center+250])
% xticklabels({'-500','-250','0','250','500'}),ylim([0,0.5]);
% 
% figure;
% plot(logData{1}.sig3d(:,end),'c');hold on;plot(logData{2}.sig3d(:,end),'b');hold on;
% plot(logData{3}.sig3d(:,end),'r');hold on;plot(logData{4}.sig3d(:,end),'k');
% legend('CSF','WM','GM','Fat')
% xticks([center-250,center-125,center,center+125,center+250])
% xticklabels({'-500','-250','0','250','500'}),ylim([0,0.5]);
