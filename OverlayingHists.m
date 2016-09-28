function OverlayingHists(data1,data2)


hist(data1)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on


hist(data2, 49)
h = findobj(gca,'Type','patch');
set(h,'facealpha',0.75);
