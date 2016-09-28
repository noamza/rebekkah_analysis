function fourier_mat_new=CreateFourierTransformMat(fourier_mat,x1,x2,y1,y2)

% smooth to find peaks
sm_fourier_mat= abs(fourier_mat);
parms.sigma=3;
sm_fourier_mat=SmoothRateMat(sm_fourier_mat,parms);

% remove outside of the main square
fm_inds=FindMaxIndsRateMap(sm_fourier_mat);
fm_inds_remove=(fm_inds(:,1)<x1 | fm_inds(:,1)>x2 ...
    |fm_inds(:,2)<y1 | fm_inds(:,2)>y2);
fm_inds(fm_inds_remove,:)=[];

% remove too close to the center
mp=ceil(size(fourier_mat)/2);
mp=mp+[-9 9];

fm_inds_remove=(fm_inds(:,1)>mp(1) & fm_inds(:,1)<mp(2) ...
    & fm_inds(:,2)>mp(1) & fm_inds(:,2)<mp(2));
fm_inds(fm_inds_remove,:)=[];

abs_fourier=abs(fourier_mat);

% figure; imagesc(abs_fourier)
% hold on;
% plot(fm_inds(:,2),fm_inds(:,1),'kx')

fm_inds=RemoveTooCloseMaxInds(fm_inds, 8, abs_fourier, 1.9);

rates=nan(1,length(fm_inds));
for h=1:length(fm_inds)
    rates(h)=abs_fourier(fm_inds(h,1), fm_inds(h,2));
end

%
%     figure; imagesc(abs_fourier)
% hold on;
% plot(fm_inds(:,2),fm_inds(:,1),'kx')

% remove inds too close to center according to max ind distance
middle_pt=size(fourier_mat)/2;
max_rate_ind=find(rates==max(rates));
max_rate_ind=max_rate_ind(1);
circle_dist=Distance(middle_pt(1), middle_pt(2),fm_inds(max_rate_ind,1), fm_inds(max_rate_ind,2));

for h=1:length(fm_inds)
    dist=Distance(middle_pt(1), middle_pt(2),fm_inds(h,1), fm_inds(h,2));
    
    if dist < circle_dist*0.7
        fm_inds(h,:)=[nan nan];
        rates(h)=nan;
    end
    
end

remove=isnan(fm_inds(:,1));
fm_inds(remove,:)=[];
rates(isnan(rates))=[];

orig_max_rate=max(rates);
idxs=[];
for h=1:3
    max_rate_ind=find(rates==max(rates));
    
    if length(max_rate_ind)==2
        
        if max(rates) > orig_max_rate*0.5
            
            idxs(end+1)=max_rate_ind(1);
            idxs(end+1)=max_rate_ind(2);
            
            rates(idxs([end-1,end]))=0;
            
        end
        
    else
        rates(max_rate_ind)=0;
        h=h-1;
        
    end
end

len=length(idxs);
size_fm=size(fourier_mat);
rad=8;

fourier_mat_new=zeros(size_fm);

% if distance to max pt is less than PF radius, assign it value at max pt
for cen=1:len;
    for fig_i =1:size_fm(1)
        for j =1:size_fm(2)
            if Distance(fig_i, j, fm_inds(idxs(cen),1), fm_inds(idxs(cen),2)) < rad  %change this depending on how large you want fields to be
                fourier_mat_new(fig_i,j)= fourier_mat(fig_i,j);
            end
        end
    end
end


