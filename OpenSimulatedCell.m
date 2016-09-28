load('simulated results 250cells 500x per cell.mat')

for sim=1:500
    
     fanos=[];
    dists=[];
    
    for c= 1:250
    
    all_sims= all_fano_factors{c};
    all_sim_dists= all_simulated_dists{c};
    
    fanos(c)=all_sims(sim);
    dists(c)=all_sim_dists(sim);
    end
    
    median_fano(sim)=median(fanos);
    sum_dists(sim)=sum(dists<=0.1) ;
       
end
    