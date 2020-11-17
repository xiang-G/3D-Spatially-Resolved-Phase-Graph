% function [logData]=clean_signal(logData,noise_threshold)
%
function [slogData]=clean_signal(slogData,noise_threshold)

disp(['Before clean, num of spins = ',int2str(length(slogData.signal))]);
sig_mag = abs(slogData.signal(:,1)+1i*slogData.signal(:,2));

NE=length(slogData.NKperEcho)-1;
NKperEcho=slogData.NKperEcho;
erased = sig_mag<noise_threshold;
num_erase(1)=0;
for i=1:NE
    num_erase(1+i) = num_erase(i)+sum(erased(1+NKperEcho(i):NKperEcho(i+1)));
end
slogData.signal(erased,:)=[];
slogData.kvector(erased,:)=[];
slogData.NKperEcho = slogData.NKperEcho-int32(num_erase)';

disp(['After clean, num of spins = ',int2str(length(slogData.signal))]);

end