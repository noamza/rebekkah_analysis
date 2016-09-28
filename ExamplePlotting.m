figure;
for h=1:6
subplot(2,3,h);

sorted_peak_rates_b=[];
sorted_peak_rates_b=sort(peak_rates_b{h});

 for r=1:length(sorted_peak_rates_b)
        sort_inds(r)= find(sorted_peak_rates_b(r) == peak_rates_b{h});
 end
    
 sorted_peak_rates_e=[];   
sorted_peak_rates_e=peak_rates_e{h(sort_ind)};

plot(1:length(peak_rates_b{inds(h)}), sorted_peak_rates_b, 'o-'); 
hold on;
plot(1:length(peak_rates_e{inds(h)}), sorted_peak_rates_e, 'ro-');
end