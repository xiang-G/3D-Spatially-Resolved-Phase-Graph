% function [logData]=clean_signal(logData,noise_threshold)
%
function [logData]=clean_signal(logData,noise_threshold)

disp(['Before clean, num of spins = ',int2str(length(logData.signal))]);
sig_mag = abs(logData.signal(:,1)+1i*logData.signal(:,2));
logData.signal(sig_mag<noise_threshold,:)=0;
disp(['After clean, num of spins = ',int2str(length(logData.signal))]);

end