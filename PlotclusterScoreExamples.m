%% open the examples

function PlotclusterScoreExamples

figure;
x=4;
y=5;

num=11;
%dat=load('214_06-03-10.mat_cell1.mat'); %example cluster 0
dat=load('214_06-04-21.mat_cell1.mat'); %example cluster 0
PlotFiveRateMaps(dat.rate_mats.arena1, dat.rate_mats.arena2, dat.rate_mats.arena3, dat.rate_mats.arena4, dat.rate_mats.arena5,x,y,num, dat.anchor_indices)


num=1;
%dat=load('217_06-05-23.mat_cell2.mat'); %example cluster 0.1
dat=load('217_06-07-05.mat_cell1.mat'); %example cluster 0
PlotFiveRateMaps(dat.rate_mats.arena1, dat.rate_mats.arena2, dat.rate_mats.arena3, dat.rate_mats.arena4, dat.rate_mats.arena5,x,y,num, dat.anchor_indices)


figure;

num=11;
%dat=load('236_06-09-10.mat_cell1.mat') ;%example cluster 0.3
dat=load('217_06-08-01.mat_cell1.mat'); %example cluster 0
PlotFiveRateMaps(dat.rate_mats.arena1, dat.rate_mats.arena2, dat.rate_mats.arena3, dat.rate_mats.arena4, dat.rate_mats.arena5,x,y,num, dat.anchor_indices)


num=1;
dat=load('216_06-08-30.mat_cell2.mat') ;%example cluster 0.4
dat=load('214_06-03-10.mat_cell1.mat'); %example cluster 0
PlotFiveRateMaps(dat.rate_mats.arena1, dat.rate_mats.arena2, dat.rate_mats.arena3, dat.rate_mats.arena4, dat.rate_mats.arena5,x,y,num, dat.anchor_indices)


figure;
num=11;
dat=load('214_06-04-04.mat_cell5.mat') ;%example cluster 0.6
PlotFiveRateMaps(dat.rate_mats.arena1, dat.rate_mats.arena2, dat.rate_mats.arena3, dat.rate_mats.arena4, dat.rate_mats.arena5,x,y,num, dat.anchor_indices)


num=1;
dat=load('216_06-08-21.mat_cell1.mat') ;%example cluster 1
PlotFiveRateMaps(dat.rate_mats.arena1, dat.rate_mats.arena2, dat.rate_mats.arena3, dat.rate_mats.arena4, dat.rate_mats.arena5,x,y,num, dat.anchor_indices)


%% PLOT a 0 cluster score example


function PlotFiveRateMaps(first, second, third, fourth, fifth,x,y, num, anchor_indices)

subplot(x,y,num)
imagesc(first); axis equal; axis off;

subplot(x,y,num+1)
imagesc(second); axis equal; axis off;

subplot(x,y,num+2)
imagesc(third); axis equal; axis off;

subplot(x,y,num+3)
imagesc(fourth); axis equal; axis off;

subplot(x,y,num+4)
imagesc(fifth); axis equal; axis off;

% ......................................... plot transparent rate map with
% max peak crossed

subplot(x,y,num+5)
imagesc(first, 'AlphaData', 0.5)
axis equal; axis off;
hold on;
[sizex,sizey]= size(first);
index_x= anchor_indices(1,1)*sizex;
index_y=anchor_indices(1,2)*sizey;
plot(index_y,index_x, 'x', 'MarkerSize', 15, 'LineWidth', 3, 'color', 'k');

subplot(x,y,num+6)
imagesc(second,'AlphaData', 0.5)
axis equal; axis off;
hold on;
[sizex,sizey]= size(second);
index_x= anchor_indices(2,1)*sizex;
index_y=anchor_indices(2,2)*sizey;
plot(index_y,index_x, 'x', 'MarkerSize', 15, 'LineWidth', 3, 'color', 'k');


subplot(x,y,num+7)
imagesc(third,'AlphaData', 0.5)
axis equal; axis off;
hold on;
[sizex,sizey]= size(third);
index_x= anchor_indices(3,1)*sizex;
index_y=anchor_indices(3,2)*sizey;
plot(index_y,index_x, 'x', 'MarkerSize', 15, 'LineWidth', 3, 'color', 'k');

subplot(x,y,num+8)
imagesc(fourth,'AlphaData', 0.5)
axis equal; axis off;
hold on;
[sizex,sizey]= size(fourth);
index_x= anchor_indices(4,1)*sizex;
index_y=anchor_indices(4,2)*sizey;
plot(index_y,index_x, 'x', 'MarkerSize', 15, 'LineWidth', 3, 'color', 'k');

subplot(x,y,num+9)
imagesc(fifth,'AlphaData', 0.5)
axis equal; axis off;
hold on;
[sizex,sizey]= size(fifth);
index_x= anchor_indices(5,1)*sizex;
index_y=anchor_indices(5,2)*sizey;
plot(index_y,index_x, 'x', 'MarkerSize', 15, 'LineWidth', 3, 'color', 'k');

end

end