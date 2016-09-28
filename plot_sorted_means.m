parms.dir_load_data= 'C:\Users\Dori\Desktop\sorted rates';
%parms.dir_save_images= 'N:\users\rebekkah\bin size 6 nonsmooth\suppl info all cells basic parameters';

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;

fig=figure;

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat=load(file_name);              % opens results for rate mats
    A=dat.S;
    name= file_name(1:12); %9-19
    name= strrep(name, '_', ' ');
    
    n=4;
    m=2;
%     
%     pos_mean_x=(A.pos.x + A.pos.x2)/2;
%     pos_mean_y=(A.pos.y + A.pos.y2)/2;
%         
%     spk_x=interp1(A.pos.t,pos_mean_x,A.spk.t);
%     spk_y=interp1(A.pos.t,pos_mean_y,A.spk.t);
%        
%     a= subplot(n,m,count);
%     plot(pos_mean_x,pos_mean_y,'k','MarkerSize', 5);hold on;
%     plot(spk_x,spk_y,'.r', 'MarkerSize', 5);
%     axis equal; axis off;
%     axis ij
%     name= file_name(9:19); %9-19 
%     name2=file_name(21:end-4);
%   %  name3=file_name(29:end-4)
%     name= strrep(name, '_', ' ');
%     name2= strrep(name2, '_', ' ');
%    % name3= strrep(name3, '_', ' ');
%     title(sprintf('%s %s',name, name2), 'Interpreter', 'none', 'FontName', 'calibri', 'fontsize', 8, 'horizontalAlignment', 'right');
%     % set(gca,'xtick',[]);
%      % set(gca,'ytick',[])
%     
%     count= count+1;
%     
%     max_freq= max(max(A.rate_mat));
%     
%     b= subplot(n,m,count);
%     imagesc(A.rate_mat)
%     title(sprintf('%0.1f Hz', max_freq), 'FontName', 'calibri', 'horizontalAlignment', 'left', 'fontsize', 8);
%     axis square; axis off; hold on;
%    
%      plot(A.max_index(2), A.max_index(1), 'x', 'MarkerSize', 17, 'LineWidth', 5, 'color', 'k');
%     
%     count= count+1;
%     
%     c=subplot(n,m,count);
%     imagesc(A.zone_mat)
%     axis square; axis off;
%     title('');
%     
%     count=count+1;
%     
%     p_1 = get(a, 'pos');
%     p_2 = get(b, 'pos');
%     p_3 = get(c, 'pos');
% 
%     
%     p_1(1)= 0.2;
%     p_2(1)= 0.28;
%     p_3(1)= 0.36;
% 
%     
%     set(a, 'pos', p_1)
%     set(b, 'pos', p_2)
%     set(c, 'pos', p_3)
% 
%     
%     disp('');
%     
    
%     subplot(n,m,count)
%     imagesc(A.rate_mats.arena1);
%     axis equal; axis off;
%     title(sprintf('%s',name), 'Interpreter', 'none', 'FontName', 'calibri', 'fontsize', 9, 'horizontalalignment', 'right')
%     
%     
%     count=count+1;
%     subplot(n,m,count)
%     imagesc(A.rate_mats.arena2);
%     axis equal; axis off;
%     
%     count=count+1;
%     subplot(n,m,count)
%     h=imagesc(A.rate_mats.arena3);
%     axis equal; axis off;
%     
%     
%     count=count+1;
%     subplot(n,m,count)
%     h=imagesc(A.rate_mats.arena4);
%     axis equal; axis off;
%     
%     
%     count=count+1;
%     subplot(n,m,count)
%     h= imagesc(A.rate_mats.arena5);
%     axis equal; axis off;
%     
%     
%     count=count+1;
%     %...........................................
%     
%     
%     
%     subplot(n,m,count)
%     h=imagesc(A.rate_mats.arena1);
%     axis equal; axis off;
%     
%     set(h,'AlphaData',0.6);
%     
%     S= get_zones(A.rate_mats.arena1, A.S);
%     plot(S.max_index(2), S.max_index(1), 'x', 'MarkerSize', 20, 'LineWidth', 5, 'color', 'k');
%     first_index= S.max_index/100;
%     
%     
%     count=count+1;
%     subplot(n,m,count)
%     h= imagesc(A.rate_mats.arena2);
%     axis equal; axis off;
%     set(h,'AlphaData',0.6);
%     
%     S= get_zones(A.rate_mats.arena2, A.S);
%     plot(S.max_index(2), S.max_index(1), 'x', 'MarkerSize', 20, 'LineWidth', 5, 'color', 'k');
%     second_index= S.max_index/100;
%     
%     count=count+1;
%     subplot(n,m,count)
%     h=imagesc(A.rate_mats.arena3);
%     axis equal; axis off;
%     set(h,'AlphaData',0.6);
%     
%     S= get_zones(A.rate_mats.arena3, A.S);
%     plot(S.max_index(2), S.max_index(1), 'x', 'MarkerSize', 20, 'LineWidth', 5, 'color', 'k');
%     third_index= S.max_index/100;
%     
%     count=count+1;
%     subplot(n,m,count)
%     h=imagesc(A.rate_mats.arena4);
%     axis equal; axis off;
%     set(h,'AlphaData',0.6);
%     
%     S= get_zones(A.rate_mats.arena4, A.S);
%     plot(S.max_index(2), S.max_index(1), 'x', 'MarkerSize', 20, 'LineWidth', 5, 'color', 'k');
%     fourth_index= S.max_index/100;
%     
%     count=count+1;
%     subplot(n,m,count)
%     h= imagesc(A.rate_mats.arena5);
%     axis equal; axis off;
%     set(h,'AlphaData',0.6);
%     
%     S= get_zones(A.rate_mats.arena5, A.S);
%     plot(S.max_index(2), S.max_index(1), 'x', 'MarkerSize', 20, 'LineWidth', 5, 'color', 'k');
%     fifth_index= S.max_index/100;
%     
%     count=count+1;
%     
%     disp('')
        fano_factor= var(A.sorted_means)/mean(A.sorted_means);
    
        subplot(n,m,count);
        plot(A.sorted_means, 'o-')
        ymax= max(A.sorted_means);
        xmax= length(A.sorted_means)+1;
        axis ([0 xmax 0 ymax]);  axis square; %axis off;
        hold on;
        title(sprintf('Fano factor= %0.2f', fano_factor), 'FontName', 'calibri', 'fontsize', 11);
        ylabel('firing rate', 'FontName', 'calibri', 'fontsize', 8);
        set(gca,'xtick',[])
    
        count= count+1;
    
end

% p = get(fig, 'pos', 'units', 'inches');
% p(3)= 7;
% p(4)= 11;
% 
% set(fig, 'pos', p);