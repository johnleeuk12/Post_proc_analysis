
tracks = unique([D.unit_cf.Track]);



edges = logspace(3,4.5,50);

figure(1)
% for t = tracks
%     X = [D.unit_cf(find([D.unit_cf.Track] == t)).cf];
%     histogram(X,edges);
%     hold on
% end

histogram([D.unit_cf.cf],edges)
xticks([1e3 : 2.5*1e3:3*1e4])
% set(gca,'xticklabel',num2str(get(gca,'xtick')/1e3','%.1f'))
xticklabels(num2str(get(gca,'xtick')'/1e3))
xlabel('kHz')
set(gca,'xscale','log')