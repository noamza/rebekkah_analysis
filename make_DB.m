function makeTestingDB(excel_file,start_num)
%
%
% makeHairpinDB
%
%
% read an excel file containing the file and position data, and generate a
% data-base containing this information.
%
% excel_file - name of file containing unit data
% start_num - starting from the <start_num> line in the excel_sheet
% (deafult = 1)
%
% 2nd version, which keeps the data structure simple
%
% Excel contains the following columns:
% rat       - number of rat
% date      - run date
% A1        - suffix of first set of files (open field)
% A2        - suffix of 4th set of files (open field)
% exp_type - info on the type of experiment
% food      - info on food administration
% tetrode   - tetrode number (i.e. t6)
% cell      - cell number (i.e. t6c1)
% arena     - name of arena (room10old, room9 etc...)
% histology 
% ... additional optional fields ...
%
% Created by Dori Derdikman, CBM, Feb 2008

if ~exist('start_num')
    start_num = 1;
end

excel_file= 'N:\users\rebekkah\AxDB\EC_cells';

% various parameters

parms.data_dir = 'N:\users\rebekkah\AxDB';  % root of data directory
parms.dim = 160; % dimensions of box + 10
parms.excel_file = excel_file;  % excel file containing the data-base data

% read excel sheet

excel_sheet = read_excel(excel_file);

% create db structure from excel sheet

db_all = build_db_from_excel(excel_sheet,parms);

% read spike data and pos data into db structure

ncells = length(db_all);

hProg = progbar( 'start' ); Progress = 0;

for cell_num = start_num: ncells

    
    rec = db_all(cell_num);
    
    % read pos data
    
    rec = HP_read_pos_data(rec,parms);
               
    % read spike data
 
    for AB = {'A1' 'A2'}
         
        AB = AB{:};
        rec = HP_read_spike_data(AB,rec,parms);
    end
    
%     db_all(cell_num) = rec;
%     save('HP_db','db_all');
    
    
    ind = find(excel_file == '.');
    excel_prefix = excel_file(1:ind-1);
    
    num_str = sprintf('%03d',cell_num);
    
    key = [excel_prefix '_' num_str '_' rec.rat '_' rec.date '_t'  ...
           num2str(rec.tetrode) 'c' num2str(rec.cell) ];
    db = rec;
    file_name = ['OUT_' key];
    save(file_name,'db','parms','db_all');
    disp(['===> save cell in file ' file_name ' <===']);
    
    Progress = cell_num/ncells; progbar( 'update', hProg, Progress );
 
    pack;
end               

Progress = 1; progbar( 'close', hProg, Progress );

disp('finished');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function db = build_db_from_excel(excel_sheet,parms)
%
% copy fields from excel sheet to DB
%
for i = 1:length(excel_sheet)
    
    db(i) = excel_sheet(i);
    db(i).rat = excel_sheet(i).rat;
    db(i).date = excel_sheet(i).date;
    db(i).exp_type = excel_sheet(i).exp_type;
    db(i).food = excel_sheet(i).food;
    db(i).tetrode = excel_sheet(i).tetrode;
    db(i).cell = excel_sheet(i).cell;
    db(i).arena = excel_sheet(i).arena;
    
    db(i).tetrode = str2num(db(i).tetrode(2)); % this is done because of the strange way it is kept in the excel
    db(i).cell = str2num(db(i).cell(4)); % assume < 10 cells (which is ok for my data).
    
    %
    % determine files, cut files, sessions, and cut_file in structure(AB)
    % later pos_data and spike_data will also be put in this structure.
    %
    
    for AB = {'A1' 'A2'}
        AB = AB{:};
        db(i).(AB).sessions = excel_sheet(i).(AB); % A1 in DB is a structure etc.
        db(i).(AB).cut_file = [db(i).date db(i).(AB).sessions(1:2) ...
            '_' num2str(db(i).tetrode) '.cut'];
        if length(db(i).(AB).sessions) == 5
            db(i).(AB).files{1} = [db(i).date db(i).(AB).sessions(1:2)];
            db(i).(AB).files{2} = [db(i).date db(i).(AB).sessions(4:5)]; % assuming sessions is of the
            % form 03,04
        elseif length(db(i).(AB).sessions) == 2 % whitlock
            db(i).(AB).files{1} = [db(i).date db(i).(AB).sessions(1:2)];
        end
    end

end % enumeration on excel sheet



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function excel_sheet = read_excel(excel_file)
%
% read excel sheet into structure
%

% read excel sheet

[num, txt, raw] = xlsread(excel_file);

fields = deblank(raw(1,:));
raw = raw(2:end,:); 

% % get rid of empty rows % does not work because of a matlab bug
% 
% full_rows = find(~isnan([raw{:,1}])); % find non-empty rows (in which first cell is a NaN after reading)
% raw = raw(full_rows,:);

% get rid of empty columns

keep_columns = find(cellfun('isclass',fields,'char'));
fields = fields(keep_columns);
raw = raw(:,keep_columns);
       
% turn all numeric values into text

num_vals = cellfun('isclass',raw,'double');
num_inds = find(num_vals);

for ind = num_inds'
        if isnan(raw{ind})
            raw{ind} = '';
        else
            raw{ind} = num2str(raw{ind});
        end
end

excel_sheet = cell2struct(deblank(raw),fields,2);
