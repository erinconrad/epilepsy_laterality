function stats = unpaired_plot_tle(a,b,xlabs,ylab,which_test,no_plot)

if ~exist('no_plot','var')
    no_plot = 0;
end

na = length(a);
nb = length(b);
if ~no_plot
    plot(1+randn(na,1)*0.05,a,'o','linewidth',2)
    hold on
    plot(2+randn(nb,1)*0.05,b,'o','linewidth',2)
    xticks([1 2])
    xticklabels(xlabs)
    ylabel(ylab)

end

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
        stats.ns = [sum(~isnan(a)) sum(~isnan(b))];
        pooled_sd = sqrt(((stats.ns(1)-1)*stats.sd(1)^2+...
            (stats.ns(2)-1)*stats.sd(2)^2)/...
            (stats.ns(1)+stats.ns(2)-2));
        cohen_d = (stats.means(1)-stats.means(2))/pooled_sd;
        stats.cohen_d = cohen_d;

end
max_all = max([max(a) max(b)]);
min_all = min([min(a) min(b)]);
span = max_all-min_all;

if ~no_plot
    ylim([min_all - 0.1*span max_all+0.1*span])
    yl = ylim;
    ybar = yl(1)+(yl(2)-yl(1))*1.1;
    ytext = yl(1)+(yl(2)-yl(1))*1.15;
    ylnew = [yl(1) yl(1)+(yl(2)-yl(1))*1.2];
    ylim(ylnew)
    plot([1 2],[ybar,ybar],'k','linewidth',2)
    text(1.5,ytext,sprintf('%s, \\itd\\rm = %1.2f',get_p_text_el(p),cohen_d),...
        'fontsize',25,'horizontalalignment','center')
    set(gca,'fontsize',25)
end

end