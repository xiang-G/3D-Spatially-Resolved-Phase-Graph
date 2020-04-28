% example_PRESS_MRSI.m
% conventional spoiler scheme vs. DOCTOPS spoiler scheme
% author: xiang gao
% xiang.gao@uniklinik-freiburg.de

clear;clc;clear global;
global FOV;
global Mat;
global Kmax2Pi;     %Skope of 3D-K space
global TDImage;     %Output 3D image or 1D signal intensity
TDImage=true;
ConventionalSpoiler = true;
saveimage=false;

%==========================================================================
%% sequence preparation
%==========================================================================
gamma = 42.5756; %MHz/T
FOV = 48; %mm
Mat = 16;
Kmax2Pi= 1.5; %default Kmax2Pi=3;
Position = [];%phantomP;
RF_BandWidth = 1200; %Hz
RF_Duration = 5; %ms
Gs = RF_BandWidth/gamma/FOV;%mT/(m*ms)
KGs = Gs*RF_Duration*gamma;
TE1hf = 7; %ms
TE2hf = 8; %ms
KSpoil1 = 10000./FOV; %1/m
KSpoil2 = 10000./FOV; %1/m
if ConventionalSpoiler
    Gx = [-KSpoil1/TE1hf;         -KSpoil1/TE1hf;         -KSpoil2/TE2hf;         -KSpoil2/TE2hf;0;0]./gamma;
    Gy = [(-KSpoil1-KGs/2)/TE1hf; (-KSpoil1-KGs/2)/TE1hf; -KSpoil2/TE2hf;         -KSpoil2/TE2hf;0;0]./gamma;
    Gz = [-KSpoil1/TE1hf;         -KSpoil1/TE1hf;         (-KSpoil2-KGs/2)/TE2hf; (-KSpoil2-KGs/2)/TE2hf;0;0]./gamma;
else
    Gx = [-KSpoil1/TE1hf;         -KSpoil1/TE1hf;         -KSpoil2/TE2hf;         -KSpoil2/TE2hf;0;0]./gamma;
    Gy = [(-KSpoil1-KGs/2)/TE1hf; (-KSpoil1-KGs/2)/TE1hf; -KSpoil2/TE2hf;         -KSpoil2/TE2hf;0;0]./gamma;
    Gz = [-KSpoil1/TE1hf;         0/TE1hf;                (0-KGs/2)/TE2hf;        (KSpoil2-KGs/2)/TE2hf;0;0]./gamma;
end
dG = 0.00*Gs.*[1,2,3];
Gx = Gx+dG(1);  Gy = Gy+dG(2);  Gz = Gz+dG(3);

%==========================================================================
%% simulation preparation
%==========================================================================
NPulse=6;
sequence.npulse=NPulse;
sequence.time=[TE1hf;TE1hf;TE2hf;TE2hf;20;40];% ms
sequence.echo=[TE1hf;TE1hf;TE2hf;TE2hf;20;40];% ms
sequence.gradx=Gx;
sequence.grady=Gy;
sequence.gradz=Gz;
sequence.velocityx=zeros(1,NPulse);     % transport velocity in mm/s
sequence.velocityy=zeros(1,NPulse);     % transport velocity in mm/s
sequence.velocityz=zeros(1,NPulse);     % transport velocity in mm/s
sequence.angle=[90/2;180/2;0;180/2;0;0];
sequence.axes=[0,90*ones(1,NPulse-1)];  % rotation axis.
sequence.diff=0;                        % mcm^2/ms  %2 for water
sequence.T1=0;                          % ms, 0 equivalent to infinity
sequence.T2=0;                          % ms, 0 equivalent to infinity
sequence.Omega=0;                       % Hz, frequency offset
sequence.ktolerance=1.e-8;              % computational parameter:
sequence.ktolerance_phys=1.e-6;         % 3d k-vectors gridding size 1e-3 = 1*(1/mm)
sequence.nOutput=1;
matlabseq = 'matlabseq_1.txt';
matlab2c(sequence, matlabseq);
if ismac
    string = ['./3D_SR-PG_ios.out ',matlabseq];    
elseif isunix
    string = ['./3D_SR-PG.out ',matlabseq];
elseif ispc
    string = ['3D_SR-PG.exe ',matlabseq];
else
    disp('Platform not supported')    
end

%==========================================================================
%% simulation execution && signal calculation
%==========================================================================
tic
system(string);
toc
logData=KmapSignal2(['output',int2str(sequence.nOutput),'.txt']);%faster, but need preallocate memory(!).
logData.sequence = sequence;
logData=IMG_DFT(logData,Position);

% 3D-Kmap
% 3d k-vectors
dKx = 1/FOV*1e3;%1/m  (post-gradding)
dKy = 1/FOV*1e3;%1/m
dKz = 1/FOV*1e3;%1/m
NKperEcho=logData.NKperEcho;
kvector = logData.kvector./2/pi*1e6; %1/m
signal = abs(logData.signal(:,1)+1i*logData.signal(:,2));
if Kmax2Pi
    Kmax=[1,1,1]*Kmax2Pi*1e3; %1/m
else
    Kmax=max(abs(kvector(1:end,:))); %1/m
    disp(Kmax)
end
Kgrid3D=zeros(2*ceil(Kmax(1)/dKx)+1,2*ceil(Kmax(2)/dKy)+1,2*ceil(Kmax(3)/dKz)+1,numel(NKperEcho)-1);
for necho=1:numel(NKperEcho)-1
    Kgrid=zeros(2*ceil(Kmax(1)/dKx)+1,2*ceil(Kmax(2)/dKy)+1,2*ceil(Kmax(3)/dKz)+1);
    shift=[abs(ceil(Kmax(1)/dKx)+1),abs(ceil(Kmax(2)/dKy)+1),abs(ceil(Kmax(3)/dKy)+1)];
    for i=NKperEcho(necho)+1:NKperEcho(necho+1)
        bins = [myfloor(kvector(i,1)/dKx),myfloor(kvector(i,2)/dKy),myfloor(kvector(i,3)/dKz)]+shift;
        Kgrid(bins(1),bins(2),bins(3))=Kgrid(bins(1),bins(2),bins(3))+signal(i);
    end
    Kgrid3D(:,:,:,necho) = Kgrid;
end
% receive zone
background = zeros(size(Kgrid3D));
bg1_st = ceil((size(Kgrid3D,1)-Mat)/2);
bg2_st = ceil((size(Kgrid3D,2)-Mat)/2);
bg3_st = ceil((size(Kgrid3D,3)-Mat)/2);
if Mat<=1
    for necho=1:numel(NKperEcho)-1
        background(bg1_st+1,bg2_st+1,bg3_st+1,necho)=1;
    end
else
    for necho=1:numel(NKperEcho)-1
        background(bg1_st:Mat/2:bg1_st+Mat,bg2_st:Mat/2:bg2_st+Mat,bg3_st:Mat/2:bg3_st+Mat,necho)=1;
    end
end

%==========================================================================
%% OUTPUT
%==========================================================================
% 3D kspace graph
figure;hold on;
set(gcf,'outerposition',[0 0 1000 800]);
for i=1:numel(NKperEcho)-1
    subplot(2,ceil((numel(NKperEcho)-1)/2),i)
    plot_3dkmap(squeeze((Kgrid3D(:,:,:,i))),background(:,:,:,i)./10,i)%3D(x,y,z) on third echo
end
hold off
if saveimage
    if ConventionalSpoiler
        print(gcf, '-dpng', '-r150', ['conventional_mat',mat2str(Mat),'_K.png']);
    else
        print(gcf, '-dpng', '-r150', ['new_mat',mat2str(Mat),'_K.png']);
    end
end

% projected images
figure;
set(gcf,'outerposition',[0 0 500 400]);
for i=1:3*(NPulse-1)
    subplot(3,NPulse-1,i)
    imagesc(single(logData.sig3d(:,:,i)));
    title(['Time ',int2str(i)]);
    if i==1
        ylabel('XY','FontWeight','bold');
    elseif i==NPulse
        ylabel('YZ','FontWeight','bold');
    elseif i==2*NPulse-1
        ylabel('XZ','FontWeight','bold');
    end
    xticks(1:(Mat-1)/2:Mat);xticklabels({int2str(-FOV/2),'0',int2str(FOV/2)});
    yticks(1:(Mat-1)/2:Mat);yticklabels({int2str(-FOV/2),'0',int2str(FOV/2)});
    colormap gray;
end

if saveimage
    if ConventionalSpoiler
        print(gcf, '-dpng', '-r150', ['conventional_mat',mat2str(Mat),'_I.png']);
    else
        print(gcf, '-dpng', '-r150', ['new_mat',mat2str(Mat),'_I.png']);
    end
end

return
