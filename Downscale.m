function [zero_inds, non_zero_inds] = Downscale(zero_inds, non_zero_inds, number_of_PF, place_field_num)    

    size_zero_inds= length(zero_inds);
    size_nonzero_inds= length(non_zero_inds);
    
    PF_num_fire= place_field_num;
    PF_num_fire(zero_inds)= nan;
    
    PF_num_nonfire= place_field_num;
    PF_num_nonfire(non_zero_inds)= nan;
    
    PF_num_fire_orig= PF_num_fire;
    PF_num_nonfire_orig=PF_num_nonfire;
    
    to_be_removed=[];
    
    for h=1:number_of_PF
        fire_inds=[];
        nonfire_inds=[];
        fire_inds= find(PF_num_fire==h);
        nonfire_inds=find(PF_num_nonfire==h);
        PF_fire(h)= length(fire_inds);
        PF_nonfire(h)= length(nonfire_inds);
        if PF_fire(h)< PF_nonfire(h)
            place_field_inds=[];
            downscale_inds=[];
            diff= length(nonfire_inds)-length(fire_inds);
            downscale_inds = randperm(length(nonfire_inds), diff);
            to_be_removed=nonfire_inds(downscale_inds);
            PF_num_nonfire(to_be_removed)=nan;
        elseif PF_fire(h)>PF_nonfire(h)
            place_field_inds=[];
            downscale_inds=[];
            diff= length(fire_inds)-length(nonfire_inds);
            downscale_inds = randperm(length(fire_inds), diff);
            to_be_removed=fire_inds(downscale_inds);
            PF_num_fire(to_be_removed)=nan;
        end
    end
    
    removed_inds_fire=[];
    removed_inds_nonfire=[];
    
    count=1;
    count2=1;
    for h=1:length(PF_num_fire)
        if isnan(PF_num_fire(h)) && ~isnan(PF_num_fire_orig(h))
            removed_inds_fire(count)=h;
            count=count+1;
        end
    end
    
    count=1;
    count2=1;
    for h=1:length(PF_num_nonfire)
        if isnan(PF_num_nonfire(h)) && ~isnan(PF_num_nonfire_orig(h))
            removed_inds_nonfire(count2)=h;
            count2=count2+1;
        end
    end
    
    for h=1:length(removed_inds_fire);
        non_zero_inds(non_zero_inds==removed_inds_fire(h))= nan;
    end
    non_zero_inds(isnan(non_zero_inds))= [];
    
    for h=1:length(removed_inds_nonfire);
        zero_inds(zero_inds==removed_inds_nonfire(h))= nan;
    end
    zero_inds(isnan(zero_inds))= [];

