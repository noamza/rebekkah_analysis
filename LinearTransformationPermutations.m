
% parms.dir_load_data = 'N:\users\rebekkah\barry jeffery pivot point\FULL DATA SET\RESULTS\max';
% 
% dir_name=parms.dir_load_data;
% dir_list = dir(strcat(dir_name,'\*.mat'));
% file_names = {dir_list.name};
% count=1;
% 
% % enumerate on cells
% for j =1:length(file_names)
%     cd(parms.dir_load_data);
%     file_name = file_names{j};
%     load(file_name);
%     
%     
%     max_inds1_orig= FindMaxIndsRateMap(rate_mats.arena1);
%     autocorr1= Cross_Correlation(rate_mats.arena1,rate_mats.arena1);
%     autoMaxInds1= FindAutoMaxInds(autocorr1);
%     PF_radius1= findPlaceFieldRadius(autocorr1, autoMaxInds1);
%     max_inds1_orig = RemoveTooCloseMaxIndsRescaling(max_inds1_orig, PF_radius1, rate_mats.arena1);
%     
%     max_inds2_orig= FindMaxIndsRateMap(rate_mats.arena2);
%     autocorr2= Cross_Correlation(rate_mats.arena2,rate_mats.arena2);
%     autoMaxInds2= FindAutoMaxInds(autocorr2);
%     PF_radius2= findPlaceFieldRadius(autocorr2, autoMaxInds2);
%     max_inds2_orig = RemoveTooCloseMaxIndsRescaling(max_inds2_orig, PF_radius2, rate_mats.arena2);
%    
       %testing below%
    max_inds1_orig= [1 1; 15 15; 30 30;];
    max_inds2_orig= [1 0.5; 15 7.5; 30 15];   
    
    %testing above%
    
    min_len= min(length(max_inds1_orig), length(max_inds2_orig));
   
    diff=abs(length(max_inds1_orig)-length(max_inds2_orig));
    if length(max_inds1_orig) <= length(max_inds2_orig)
        perms_n = 1:length(max_inds1_orig);
        perms_m = perms(1:length(max_inds2_orig));
        perms_m(:,end-diff+1:end)=[];
        perms_m = unique(perms_m,'rows');
    elseif length(max_inds1_orig) > length(max_inds2_orig)
        perms_n = perms(1:length(max_inds1_orig));
        perms_n(:,end-diff+1:end)=[];
        perms_n = unique(perms_n,'rows');
        perms_m = 1:length(max_inds2_orig);
    end
    
    
    A= nan(1,factorial(min_len));
    B= nan(1,factorial(min_len));
    C= nan(1,factorial(min_len));
    D= nan(1,factorial(min_len));
    
    [permsn_len, ~] =size(perms_n);
    [permsm_len, ~] =size(perms_m);
    
    count=1;
    
    min_dist=inf;
    for n = 1:permsn_len  %% add pair permutation here
        
        for m= 1:permsm_len
            
            max_inds1= max_inds1_orig(perms_n(n,:),:);
            max_inds2= max_inds2_orig(perms_m(m,:),:);
            
            x= max_inds1(:,1);
            y=max_inds1(:,2);
            
            x1= max_inds2(:,1);
            y1= max_inds2(:,2);
            
      
    
    %size_max= length(rate_mats.arena1);
    %try all As from 1 to 2
    %try all Bs from 0 to size of arena
   
    size_max=50;
    
    
    for A=1:0.2:2
        for B=0:10:size_max
            for C=1:0.2:2
                for D=0:10:size_max
                    
                    dist=0;
                    for i= 1:length(x)
                        
                    dist= dist + (x(i) - A*x1(i) -B)^2 + (y(i) - C*y1(i) -D)^2;
                    
                    
                    end
                    dist= dist/length(x);
                    
                    if dist< min_dist;
                    min_dist= dist;
                    A_min=A;
                    B_min=B;
                    C_min=C;
                    D_min=D;
                    
                    n_min=n;
                    m_min=m;
                    
                    end
                    
                end
            end
        end
    end
    
    
        end
    end
    
  
            
%             a= 1:min_len;
%             A(count)= (sum(x1(a) .* x(a) -mean(x))) / ...
%                         (sum(x(a) .* x1(a) -mean(x1)));
%             
%             B(count)= mean(x) - A(count)*mean(x1); 
%             
%             C(count)= (sum(y1(a) .* y(a) -mean(y))) / ...
%                         (sum(y(a) .* y1(a) -mean(y1)));
%             
%             D(count)= mean(y) - C(count)*mean(y1);
%             
%             count=count+1;
%             
%         end
%         
%     end
    
            max_inds1= max_inds1_orig(perms_n(n_min,:),:);
            max_inds2= max_inds2_orig(perms_m(m_min,:),:);

    x= max_inds1_orig(:,1);
    y= max_inds1_orig(:,2);
    
    x1= max_inds2_orig(:,1);
    y1= max_inds2_orig(:,2);
%     
%     
%     count=1;
%     
%     S=nan(1,length(A));
%     for c= 1:length(A)
%         S=0;
%         for b= 1:length(x)
%             for a=1:length(x1)
%                 S(c)= S+ (x(b) - A(c)*x1(a) -B(c))^2 + (y(b) - C(c)*y1(a) -D(c))^2;
%             end
%         end   
%     end
%     
%     ind= find(S == min(S));
%     ind=ind(1);
%     
%     A_min= A(ind);  %stretch in x dimension
%     B_min= B(ind);  %shift in x dimension
%     C_min= C(ind);  %stretch in y dimension
%     D_min= D(ind);  %shift in y dimension
    
    for i=1:length(x1)
        final_x1(i)= A_min*(x1(i))-B_min;
        final_y1(i)= C_min*(y1(i))-D_min;
    end
    
    max_inds_final= nan(size(max_inds2_orig));
    max_inds_final(:,1)= final_x1;
    max_inds_final(:,2)= final_y1;
    
    figure;
    plot(max_inds1_orig(:,1), max_inds1_orig(:,2), 'o', 'MarkerFaceColor', 'k'); hold on;
    plot(max_inds_final(:,1), max_inds_final(:,2), 'o', 'MarkerFaceColor', 'r');
    
%end

disp('');
