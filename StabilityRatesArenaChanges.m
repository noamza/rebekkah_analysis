parms.dir_load_data = 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\RESULTS\max';
% parms.dir_save_data = 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\stability between fields removed last arena';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};

count_all=1;

for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    load(file_name);
    
    file_count=1;
    max_inds=[{'max_inds1'};{'max_inds2'};{'max_inds3'}; {'max_inds4'}; {'max_inds5'}];
    
    for h=1:5
        
        name= sprintf('arena%d', h);
        arena_num= sprintf('%d', h);
        
        name2= sprintf('max_inds%s', arena_num);
        %[max_inds] = FindMaxIndsRateMap(rate_mats.(name));
       
        % stretch arena to size of original
        [orig_x orig_y]= size(rate_mats.arena1);
        [new_x new_y]= size(rate_mats.(name));
        
        factorx= orig_x/new_x;
        factory= orig_y/new_y;
            
        tform = maketform('affine',[factory 0 0; 0 factorx 0; 0 0 1]);
        stretch_arena = imtransform(rate_mats.(name),tform);
        
        stretch_arenas{h}=stretch_arena;
        
        max_inds{h}= FindMaxIndsRateMap(stretch_arena);
        autocorr= Cross_Correlation(stretch_arena,stretch_arena);
        autoMaxInds= FindAutoMaxInds(autocorr);
        PF_radius= findPlaceFieldRadius(autocorr, autoMaxInds);
        max_inds{h} = RemoveTooCloseMaxInds(max_inds{h}, PF_radius, stretch_arena);
        
        
        len= size(max_inds{h});
        len=len(1);
        
        rates=nan(1,length(len));
        
        for l= 1:len
            rate_mat=stretch_arenas{h};
            current_max_inds= max_inds{h};
            rates(l)= rate_mat(current_max_inds(l,1), current_max_inds(l,2));
        end
        
        all_rates{h}= rates;
        
    end
    
     max_inds_orig= max_inds{1};
    
    for h=1:5
        lenn= size(max_inds{h});
        lens(h)= lenn(1);
    end
    
    %sprintf
    
    
    zones1=length(max_inds{1});
    
    count=1;
    
    for k=5
        
        zones2= size(max_inds{k});
        zones2=zones2(1);
        
        if abs(zones1-zones2)<= 3 & zones2 ~=1  %%% change zones2 ~=1 for regular and >2 for max removed
            
            
            len= min(lens(1), lens(k));
            
            zone_val= nan([1 length(max_inds{1})]);
            rates_orig= nan([1 length(max_inds{1})]);
            rates_new= nan([1 length(max_inds{1})]);
            
            max_inds1= max_inds{1};
            max_inds2= max_inds{k};
            
            new_max_inds2=nan([length(max_inds2) 2]);
            
            for h=1:len
                
                [min_i, min_j] = FindTransfPtToPt2(max_inds1, max_inds2);
                
                new_max_inds2(min_j,1)= max_inds1(min_i,1);
                new_max_inds2(min_j,2)= max_inds1(min_i,2);
                
                max_inds1(min_i,:)= nan;
                max_inds2(min_j,:)= nan;
                
                zone_val(min_i)= min_j;
                
                
               %  rates(h)= new_max_inds2;
               % zone_vals(h) = zone_val;
               % count=count+1; 
            end
            
           % orig_rate_mat= rate_mats.arena1;
           % current_rate_mat= rate_mats.(name);
            
             rates_orig= all_rates{1}; 
             rates_orig(isnan(zone_val))=nan;
             rates_next= all_rates{k};
             rates_new= nan([1 len]);
             zone_val_nan= zone_val;
             zone_val_nan(isnan(zone_val_nan))=[];
             rates_new= rates_next(zone_val_nan);
             rates_orig(isnan(rates_orig))=[];
             
             %% TO REMOVE MAX FIELD 
             
             max_ind= find(rates_orig==max(rates_orig));
             rates_orig(max_ind)=[];
             rates_new(max_ind)=[];
             
%              %% TO REMOVE MAX FIELD //
%             
             if length(rates_orig) > 1
             rates_corr= corrcoef(rates_orig, rates_new);
             rates_corr=rates_corr(2);
             else
                 rates_corr=nan;
             end
            
            
             
%              sorted_orig=[];
%              sort_inds=[];
%              sorted_new=[];
%              
%              
%              
%              if rates_corr >= 0.85
%              sorted_orig= sort(rates_orig);
%              
%              
%              for r=1:length(rates_orig)
%                     sort_inds(r)= find(sorted_orig(r) == rates_orig);
%               end
%              
%              sorted_new=rates_new(sort_inds);
%              
%             figure; subplot(1,2,1)
%             plot(1:length(rates_orig),sorted_orig, 'o-'); hold on;
%             plot(1:length(rates_orig),sorted_new, 'ro-'); hold on;
%            title(sprintf('%s%d', file_name, file_count));
%             subplot(1,2,2); scatter(rates_orig, rates_new);
%               title(sprintf('%0.2f', rates_corr));
%              
%              end
              
%         cd(parms.dir_save_data)
%         save(sprintf('%s_file%d.mat', file_name, file_count), 'zone_val', 'new_max_inds2', ...
%             'max_inds_orig', 'orig_rate_mat', 'current_rate_mat', ...
%             'all_rates', 'rates_orig', 'rates_new', 'rates_corr');
%         cd(parms.dir_load_data);
        
        file_count=file_count+1; 
        
        rates_corr_all(count_all)= rates_corr;
            count_all=count_all+1;
            
        end
        
        
    end
end

nanmean(rates_corr_all)
sum(rates_corr_all>=0.5)
sum(rates_corr_all>=0.7)


