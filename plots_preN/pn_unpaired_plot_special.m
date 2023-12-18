function stats = unpaired_plot_special(a,b,atype,btype,xlabs,ylab,which_test)

markers = {'o','x','*','+'};

na = length(a);
nb = length(b);

assert(length(atype)==na);
assert(length(btype)==nb);

all_types = unique([atype;btype]);
ntypes = length(all_types);
assert(ntypes <=4);

[~,awhich] = (ismember(atype,all_types));
[~,bwhich] = (ismember(btype,all_types));

amarkers = markers(awhich);
bmarkers = markers(bwhich);

for i = 1:na
    plot(1+randn*0.05,a(i)+randn*0.05,'marker',amarkers{i},'linewidth',2,...
        'color',[0 0.4470 0.7410],'markersize',10)
    hold on
end
hold on
for i = 1:nb
    plot(2+randn*0.05,b(i)+randn*0.05,'marker',bmarkers{i},'linewidth',2,...
        'color',[0.8500 0.3250 0.0980],'markersize',10)
end
xticks([1 2])
xticklabels(xlabs)
yticks(1:5)
ylabel(ylab)

switch which_test
    case 'non_para'
        [p,~,stats_stuff] = ranksum(a,b);
        stats.p = p;
        stats.ns = [sum(~isnan(a)) sum(~isnan(b))];
        stats.ranksum = stats_stuff.ranksum;
    case 'para'
        [~,p,~,ostats] = ttest2(a,b);
        stats.p = p;
        stats.tstat = ostats.tstat;
        stats.df = ostats.df;
        stats.means = [nanmean(a) nanmean(b)];
        stats.sd = [nanstd(a) nanstd(b)];
end
max_all = max([max(a) max(b)]);
min_all = min([min(a) min(b)]);
span = max_all-min_all;
ylim([min_all - 0.1*span max_all+0.1*span])
yl = ylim;
ybar = yl(1)+(yl(2)-yl(1))*1.1;
ytext = yl(1)+(yl(2)-yl(1))*1.15;
ylnew = [yl(1) yl(1)+(yl(2)-yl(1))*1.4];
ylim(ylnew)
plot([1 2],[ybar,ybar],'k','linewidth',2)
text(1.5,ytext,get_p_text_el(p),'fontsize',15,'horizontalalignment','center')
text(2,yl(1)+(yl(2)-yl(1))*1.35,sprintf('%s = %s',markers{1},all_types{1}),...
    'fontsize',15,'horizontalalignment','right')
text(2,yl(1)+(yl(2)-yl(1))*1.27,sprintf('%s = %s',markers{2},all_types{2}),...
    'fontsize',15,'horizontalalignment','right')
set(gca,'fontsize',15)

end