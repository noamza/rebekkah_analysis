


    stretch_one= [0.6 0.3 0.6 1 0.1 0.3 0.4 0.6 1 0.6 0.8 0.4 0.8 1 0 0.4 1 1 ];
    stretch_three=[0 0.1 0.2 0.3 0.2 0.1 0.3 1 0.2 0.4 1 0.6];
    stretch_four=[0.6 0.6 1 0.4 0.6 0.6 0.6 0.6 0.2 0.6 0.6 0.1 0.6 0.3 0.6 0.2 0.4 0.6 0.6 0.6 0.4 0.6 0.3 0.2 1 0.4 0.6 0.3 0.3 1 0.4 0.5 0.167];
    
   stretch_means= [mean(stretch_one), mean(stretch_three), mean(stretch_five)] ;
   
   
   
   figure;

   
std_one=std(stretch_one);
%std_two=std(stretch_two);
std_three=std(stretch_three);
%std_four=std(stretch_four);
std_five=std(stretch_five);
   
        plot(1,stretch_means(1),'*m'); hold on;
        %plot(2,stretch_means(2),'*r');
        plot(2,stretch_means(2),'*g');
        %plot(4,stretch_means(4),'*k');
        plot(3,stretch_means(3),'*y');
        errorbar(1,stretch_means(1),std_one,'.m','linewidth',1);
        %errorbar(2,stretch_means(2),std_two,'.r','linewidth',1);
        errorbar(2,stretch_means(2),std_three,'.g','linewidth',1);
        %errorbar(4,stretch_means(4),std_four,'.k','linewidth',1);
        errorbar(3,stretch_means(3),std_five,'.y','linewidth',1)