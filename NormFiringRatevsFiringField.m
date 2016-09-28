
clearvars

%for real data
% cd('\\192.114.21.198\Dori_Docs\users\rebekkah\results and info of analysis');
% load('variability and border distances.mat', 'peak_rates_all_sm')
% load('variability and border distances.mat', 'peak_rates_all')

%for simulated data
 load('simulated results with fields 10.mat')

if ~exist('peak_rates_all_sm','var')
peak_rates_all_sm=sorted_rates;
end

figure;

all_norm_orders=[];
all_sorted_rates=[];

% take firing rates and field ranks and normalize to max
norm_rates=cell(1,length(peak_rates_all_sm));
norm_ranks=cell(1,length(peak_rates_all_sm));
all_norm_ranks=[];
all_norm_rates=[];
for i=1:length(peak_rates_all_sm)
    
     peak_rates= sort(peak_rates_all_sm{i});
     norm_rates{i}=peak_rates;
     % norm_rates{i}= peak_rates/ max(peak_rates);
     
     norm_ranks{i}= (1:length(peak_rates))/length(peak_rates);
     
    %convert from cell array to vector format 
    all_norm_ranks(end+1:end+length(norm_ranks{i}))=norm_ranks{i};  
    all_norm_rates(end+1:end+length(norm_rates{i}))=norm_rates{i};
end    
  
lower_bins=0:0.1:1;
upper_bins=0.1:0.1:1.1;

bin_upper_value=nan(1,length(all_norm_ranks)); 
for i=1:length(all_norm_ranks)
    
    ind=find(lower_bins<=all_norm_ranks(i)& upper_bins>all_norm_ranks(i));
    ind=ind(end);
    
    bin_upper_value(i)=upper_bins(ind);
end

mean_rates=nan(1,length(upper_bins));
mean_ranks=nan(1,length(upper_bins));
stderrs=nan(1,length(upper_bins));
for i=1:length(upper_bins);
    value= find(bin_upper_value==upper_bins(i)); 
mean_rates(i)=mean(all_norm_rates(value));
mean_ranks(i)=mean(all_norm_ranks(value));
stderrs(i)= std(all_norm_rates(value))/sqrt(length(value));
end

length(all_norm_ranks)

figure; scatter(all_norm_ranks, all_norm_rates, 'o');

% figure;
% all_norm_ranks= all_norm_ranks';
% all_norm_rates= all_norm_rates';
% f= fit(all_norm_ranks,all_norm_rates,'poly1');
% p= confint(f,0.95);
% 
% p1_1= p(1,1);
% p1_2= p(2,1);
% p2_1= p(1,2);
% p2_2= p(2,2);
% 
% x=0:2;
% y=(f.p1*x +f.p2);
% y1= (p1_1*x +p2_1);
% y2= (p1_2*x +p2_2);
% 
% h= plot(x,y, 'Color', [0.85 0.16 0], 'LineWidth', 2); hold on;
% h1= plot(x,y1,'--', 'Color', [0.85 0.16 0],'LineWidth', 2);
% h2= plot(x,y2,'--','Color', [0.85 0.16 0],'LineWidth', 2);
% text(-0.25,1.05,'C','Units', 'Normalized', 'VerticalAlignment', 'Top', 'fontw', 'b', 'fontsize', 14)
% 
% errorbar(mean_ranks, mean_rates, stderrs, 'x');hold on;
% 
% xlim([0 1.1])
% ylim([0 1.1])

% save('norm variability analysis results', 'all_norm_orders', ...
%     'all_sorted_rates', 'order', 'means_by_bins', 'stderr')


disp('');