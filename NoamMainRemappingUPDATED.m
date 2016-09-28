function NoamMainRemappingUPDATED
%Date: 05 of Apr, 2016. Updated by Rebekkah.
%WHY 5?
dbstop if error
%{ 
parms.dir_load_data = '\\192.114.21.198\Dori_Data\data\rebekkah\All rats\all cells with no pos data fixing';
parms.dir_save_pictures= '\\192.114.21.198\Dori_Data\data\rebekkah\All rats\all cells images';
parms.dir_save_images= '\\192.114.21.198\Dori_Data\data\rebekkah\All rats\images original analysis';
 
cd(parms.dir_load_data);

dir_name=parms.dir_load_data;
dir_list = dir(strcat(dir_name,'\*.mat'));
file_names = {dir_list.name};
%}
binsize = 10;
fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_a.mat',binsize);
tic; cells = load(fn); cellsNoam = cells.cells; toc;
cellsNoam([288])=[];%([36,95,229,288]) = [];%,259, 262, 266, 298]) = []; %ERROR 36:95:229: PF_radius not found %288
count = 0;
% enumerate on cells
na = 0;
for i = 1 : length(cellsNoam)%1:375 %length(file_names)
    TheCell = cellsNoam{i};
    file_name = TheCell.ind;
    fprintf('i %d ind %d,\n', i, file_name);
    add = 0;
    arena = '';
    len = 1;
    after = false;
    for k = 1:length(TheCell.middle);
        if TheCell.middle{k}.exists
            len = len + 1;
        end
    end
    bins = len - 1;
    if TheCell.after.exists
        len = len + 1;
        after = true;
    end
        %initialize variables
    zone_mats = cell(1,len);
    arena_types= cell(1,len);
    max_inds_all= nan(len,2);
    all_max_indices=cell(1,len);
    rate_mats_all_arenas=cell(1,len);
    gridness_scores_all_arenas= nan(1,len);
    MD_scores_all_arenas= nan(1,len);
    autocorrs_all_arenas=cell(1,len);
    posx_aa = cell(1,len);
    posy_aa = cell(1,len);
    post_aa = cell(1,len);
    spkx_aa = cell(1,len);
    spky_aa = cell(1,len);
    spkt_aa = cell(1,len);
    %loop through cell
    for j = 1:len;
        if j == 1           
            Cell = TheCell.before;
            arena = 'before';
            Cell.exists = true;
        elseif j <= bins + 1
            Cell = TheCell.middle{j-1};
            arena = sprintf('mid%d',j-1);
        else
            Cell = TheCell.after;
            arena = 'after';
        end
        if Cell.exists
            add = add + 1;
            pos_x = Cell.px; %(Cell.pos_data.x1); %+Cell.pos_data.x2)/2;
            pos_y = Cell.py; %(Cell.pos_data.y1); %+Cell.pos_data.y2)/2;
            pos_t = Cell.pt'; %Cell.pos_data.t;
            spk_t = Cell.st'; %%Cell.spike_data.ts;
            spk_x = Cell.sx; %interp1(pos_t,pos_x,spk_t);
            spk_y = Cell.sy; %interp1(pos_t,pos_y,spk_t);
            valid{add} = true;
            %[~,MD_score,~]=ComputeMovingDirectionalityWithSpeedThreshHold(pos_t,pos_x,pos_y,spk_t,0);
            % Create Rate Maps (smoothed and unsmoothed)
            %parms.sigma=5; parms.bin_size=4;
            rate_mat = Cell.rm;%CreateRateMap(pos_x,pos_y,pos_t,spk_x,spk_y,spk_t,parms);
            % Create AutoCorr
            autocorr = Cell.ac; %Cross_Correlation(rate_mat, rate_mat);
            %Find AutoMaxInds
            auto_max_inds = FindAutoMaxInds(autocorr);
            %noam
            if ~isempty(auto_max_inds) && size(auto_max_inds,1) > 1
                PF_radius = findPlaceFieldRadius(autocorr, auto_max_inds);
            else
                PF_radius = 0.00001; %1.414213562373095;
                auto_max_inds = [0, 0];
                warning(sprintf('ERROR %d: PF_radius not found\n', i));
                valid{add} = false;
            end
            % Find Max_Inds of smoothed & nonsmoothed rate mat
            max_inds = FindMaxIndsRateMap(rate_mat);
            strength = 1.9;
            max_inds = RemoveTooCloseMaxInds(max_inds, PF_radius, rate_mat, strength);
            if max_inds == 0 %noam
                max_inds = [1, 1];
                valid{add} = false;
            end
            % Find peak firing rate at max_inds
            peak_rates = findPeakRates(max_inds, rate_mat);
            % Find location of max-firing field (hyperfield)- nonsmoothed data
            ind = find(peak_rates == max(peak_rates));
            ind = ind(1); %in cases where same rate, take first
            max_index = max_inds(ind,:);
            size_rate_mat = size(rate_mat);
            norm_max_index = max_index ./ size_rate_mat;
            % max_indices{i}= max_inds_ns; %save for future max peak shuffling
            max_indices = max_inds; %save for future max peak shuffling
            % create zone mats for images and number_zone_mat for future shuffling
            [zone_mat, ~] = CreateZoneMat(size(rate_mat), PF_radius, max_inds, peak_rates);
            % [zone_mat_ns, number_zone_mat{i}]= CreateZoneMat(rate_mat_ns, PF_radius_ns, max_inds_ns, peak_rates_ns);
            %find gridness score
            R_outer = FindROuter(autocorr);
            if isnan(R_outer)
                valid{add} = false;
                gridness = -2;
            else
                gridness = GridnessRadius(autocorr,R_outer,i);
            end       
            % COMMENTS HERE
            %save images
            zone_mats{add} = zone_mat;
            arena_types{add} = arena;%Cell.arena;          %NOAM
            max_inds_all(add,:) = norm_max_index;
            all_max_indices{add} = max_indices;
            peak_rates_all_arenas{add} = peak_rates;
            rate_mats_all_arenas{add} = rate_mat;
            gridness_scores_all_arenas(add) = gridness;
            MD_scores_all_arenas(add) = -9999; %MD_score;    %NOAM
            autocorrs_all_arenas{add} = autocorr;
            posx_aa{add} = pos_x;
            posy_aa{add} = pos_y;
            post_aa{add} = pos_t;
            spkx_aa{add} = spk_x;
            spky_aa{add} = spk_y;
            spkt_aa{add} = spk_t;

                                    %SPLIT BY GROUP NUMBER
            if add == len 
                %set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.4, 17.4], 'PaperUnits', 'centimeters', 'PaperSize', [17.4, 17.4])
                %cd(parms.dir_save_pictures); saveas(fig,sprintf('%s.jpg', file_name(1:end-4)));cd(parms.dir_load_data);
                if any(gridness_scores_all_arenas>0.3) && any(MD_scores_all_arenas<0.15)                   
                    count = count+1;
                    filenames{count} = file_name;
                    zone_mats_all{count} = zone_mats;
                    arena_types_all{count} = arena_types;
                    norm_max_index_all{count} = max_inds_all;
                    max_indices_all{count} = all_max_indices;
                    peak_rates_all{count} = peak_rates_all_arenas;
                    rate_mats_all{count} = rate_mats_all_arenas;
                    gridness_all{count} = gridness_scores_all_arenas;
                    MD_scores_all{count} = MD_scores_all_arenas;
                    autocorrs_all{count} = autocorrs_all_arenas;
                    pos_x_all{count} = posx_aa;
                    pos_y_all{count} = posy_aa;
                    pos_t_all{count} = post_aa;
                    spk_x_all{count} = spkx_aa;
                    spk_y_all{count} = spky_aa;
                    spk_t_all{count} = spkt_aa;
                    after_all{count} = after;
                    valid_all{count} = valid;
                else
                    na = na+1; %not added
                end
                
            end
        
        end
        %disp('processing: '); %i
    end  
    
end
fprintf('not added %d\n', na);

disp('C:\Noam\Output\rebekkah\');

%cd('\\192.114.21.198\Dori_Data\data\rebekkah\rerunning analysis');
cd('C:\Noam\Output\rebekkah\');
save('remapping data info noam mids', 'rate_mats_all', 'zone_mats_all', ...
    'max_indices_all', 'peak_rates_all', 'autocorrs_all', 'filenames', ...
    'gridness_all','MD_scores_all', 'arena_types_all',...
    'pos_x_all', 'pos_y_all', 'pos_t_all', ...
    'spk_x_all', 'spk_y_all', 'spk_t_all', 'after_all', 'valid_all')
cd('C:\Noam\Dropbox\GitTechnion\rebekkah')

%%%%%%%%%%%%%%%%%%%%%%%%%END








    %         image
    %                     if add==1
    %                         fig=figure;
    %                     end
    %
    %                     n=len;m=5;
    %                     subplot(n,m,1+(5*(add-1)))
    %                     plot(pos_mean_x,pos_mean_y,'k');hold on;
    %                     plot(spk_x,spk_y,'.r');
    %                     axis equal;axis off;
    %                     axis ij;
    %                     title_name= Cell.cut_file(1:end-4);
    %                     title(sprintf('%s %s %s %s', title_name, Cell.arena, Cell.tetrode, Cell.cell), 'Interpreter', 'none');
    %
    %                     subplot(n,m,2+(5*(add-1)))
    %                     imagesc(rate_mat);
    %                     axis equal; axis off;
    %                     title(sprintf('%0.1f Hz', max(peak_rates)), 'HorizontalAlignment', 'left');
    %
    %                     subplot(n,m,3+(5*(add-1)))
    %                     imagesc(zone_mat);
    %                     title(sprintf('cell = %d', count));
    %                     axis equal; axis off;
    %
    %                     subplot(n,m,4+(5*(add-1)))
    %                     plot(1:length(peak_rates),sort(peak_rates), 'ko-');
    %                     title(sprintf('Fano factor=%0.1f', fano_factor(i)));
    %                     axis square;
    %
    %                     subplot(n,m,5+(5*(add-1)))
    %                     imagesc(autocorr);
    %                     axis equal; axis off;
    %                     title(sprintf('gridness= %0.2f', gridness));
    %
    %                     subplot(n,m,6)
    %                     imagesc(rate_mat_ns);
    %                     axis equal; axis off;
    %                     title(sprintf('%0.1f Hz', max(peak_rates_ns)),'HorizontalAlignment', 'left');
    %
    %                     subplot(n,m,7)
    %                     imagesc(zone_mat_ns);
    %                     axis equal; axis off;
    %                     title(sprintf('dist=%0.2f', norm_dist(i)));
    %
    %                     subplot(n,m,8)
    %                     plot(1:length(peak_rates_ns),sort(peak_rates_ns), 'ko-');
    %                     title(sprintf('Fano factor=%0.1f', fano_factor_ns));
    

