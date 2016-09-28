
cd('\\192.114.21.198\Dori_Data\data\rebekkah\data sets');

filename = 'list_of_cells_and_layers.xlsx';

len=631;

[ndata, text, alldata] = xlsread(filename);
Rat=cell2mat(alldata(2:len,1));
tet=cell2mat(alldata(2:len,10));
cell=cell2mat(alldata(2:len,11));
layer=alldata(2:len,4);
name=alldata(2:len,13);

for i=1:length(name)


tmp=name{i};
session(i,1:2)=tmp(end-1:end);
date(i,1:6)=tmp(end-7:end-2);

end


 
for i=1:length(name)

tmp_name=sprintf('Cell_r%d_d%s_s%s_t%d_c%d.mat',Rat(i),date(i,:),...
       session(i,:),tet(i,:),cell(i,:));
cell_list(i,1:length(tmp_name))=tmp_name;


end   
  
for i=1:length(name)

tmp_layer=layer{i};

switch(tmp_layer)
   
    case 'MEC LII'
       layer_list(i)=2;
    case 'MEC LIII'
       layer_list(i)=3;
    case 'MEC LIV'
       layer_list(i)=4;
    case 'MEC LV'
       layer_list(i)=5;
    case 'MEC LVI'
       layer_list(i)=6;
       
       
end


end   
  

save('layer_cells_list.mat','cell_list','layer_list')