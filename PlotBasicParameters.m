

% plot spike mat, rate mat, zone mat, and sorted rates plot into one
% figure.

% spike mat title: Cell name, axis off.
% rate_mat: axis off.
%
% sorted rates title: Fano factor. axis_x off, axis_y is firing rate.

parms.dir_load_data= 'N:\users\rebekkah\final data smoothed\ROTATED ARENA';
parms.dir_save_images= 'N:\users\rebekkah\bin size 6 nonsmooth\suppl info all cells basic parameters'; 

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
count=1;

fig=figure;
page=0;

% enumerate on cells
for i =1:length(file_names)
    cd(parms.dir_load_data);
    file_name = file_names{i};
    dat=load(file_name);              % opens results for rate mats
    A=dat.S;
    
    
    pos_mean_x=(A.pos.x + A.pos.x2)/2;
    pos_mean_y=(A.pos.y + A.pos.y2)/2;
        
    spk_x=interp1(A.pos.t,pos_mean_x,A.spk.t);
    spk_y=interp1(A.pos.t,pos_mean_y,A.spk.t);
    
    n=6;
    m=4;
    
    if mod(i,6) - 1 ==0
        close all
        fig=figure;
        
        count=1;
    end
    
    a= subplot(n,m,count);
    plot(pos_mean_x,pos_mean_y,'k','MarkerSize', 5);hold on;
    plot(spk_x,spk_y,'.r', 'MarkerSize', 5);
    axis square;%axis off;
    axis ij
    name= file_name(13:19); %9-19 
    name2=file_name(21:end-4);
  %  name3=file_name(29:end-4)
    name= strrep(name, '_', ' ');
    name2= strrep(name2, '_', ' ');
   % name3= strrep(name3, '_', ' ');
    ylabel(sprintf('%s %s',name, name2), 'Interpreter', 'none', 'FontName', 'calibri', 'fontsize', 8); %'horizontalAlignment', 'right');
     set(gca,'xtick',[]);
      set(gca,'ytick',[])
    
    count= count+1;
    
    max_freq= max(max(A.rate_mat));
    
    b= subplot(n,m,count);
    imagesc(A.rate_mat)
    title(sprintf('%0.1f Hz', max_freq), 'FontName', 'calibri', 'horizontalAlignment', 'left', 'fontsize', 8);
    axis square; axis off;
    
    count= count+1;
    
    c=subplot(n,m,count);
    imagesc(A.zone_mat)
    axis square; axis off;
    
    count=count+1;
    
    fano_factor= var(A.sorted_means)/mean(A.sorted_means);
    
    d=subplot(n,m,count);
    plot(A.sorted_means, 'o-')
    ymax= max(A.sorted_means);
    xmax= length(A.sorted_means)+1;
    axis ([0 xmax 0 ymax]);  axis square; %axis off;
    hold on;
    title(sprintf('Fano factor= %0.2f', fano_factor), 'FontName', 'calibri', 'fontsize', 8);
    ylabel('firing rate', 'FontName', 'calibri', 'fontsize', 8);
    set(gca,'xtick',[])

    
    
    p_1 = get(a, 'pos');
    p_2 = get(b, 'pos');
    p_3 = get(c, 'pos');
    p_4 = get(d, 'pos');
    
    p_1(1)= 0.2;
    p_2(1)= 0.28;
    p_3(1)= 0.36;
    p_4(1)= 0.44;
    
    set(a, 'pos', p_1)
    set(b, 'pos', p_2)
    set(c, 'pos', p_3)
    set(d, 'pos', p_4)
    
     count=count+1;
       
    
    
     
    if mod(count,24) ==0
        
        fig.Units= 'inches'
        fig.Position= [0 0 9 6]  
    
        
        page=page+1;
        cd(parms.dir_save_images)
        saveas(fig,sprintf('suppl_info%d.fig',page));
        saveas(fig,sprintf('suppl_info%d.jpg',page));
    end
    
end