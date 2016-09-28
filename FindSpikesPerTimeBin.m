function spikes_t = FindSpikesPerTimeBin(pos_t, spk_t) 

spikes_t=zeros(1, length(pos_t));

 for h= 1:length(spk_t)
            tmp=abs(pos_t-spk_t(h));
            [~, idx] = min(tmp);
            closest=pos_t(idx);
            spikes_t(pos_t==closest)= spikes_t(pos_t==closest)+1;
 end