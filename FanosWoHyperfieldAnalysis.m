function FanosWoHyperfieldAnalysis
%Date: 19 of Jan, 2016. Updated by Rebekkah.

dbstop if error

cd('\\192.114.21.198\Dori_Data\data\rebekkah\rescaling data info\G3MD25PF3')
load('3sets G3MD25PF3 data and results.mat','peak_rates_sm_all')

len= length(peak_rates_sm_all); 

fanos_wo_hyper=nan(1,len);

for i =1:len
    
 peak_rates=peak_rates_sm_all{i};
 peak_rates(peak_rates==max(peak_rates))=[];   
 
    % Find Fano factor
    fanos_wo_hyper(i)= var(peak_rates)/mean(peak_rates);
    
end

save('fanos withoth hyperfield', 'fanos_wo_hyper')


disp('');
