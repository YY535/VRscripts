function [p2, pg, pp]= StatePlot(bprd,ObjTypes, titleS, savepath,filename,AnimalN,varargin)%
xx = [];
c= [];
x = [];
L = [];
p = zeros(1,3);
thrP = 480; % thereshold of the period length in samples. 240 smaples/s
clear g
lb = {'pre','on','post'};
klb = {'V:', 'R:'};
p2 = nan(3);
pp = zeros(2,3);
% sAnimalN = AnimalN(AnimalN>0);
uSAni = setdiff(unique(AnimalN),0);
mAni = zeros(length(uSAni), length(AnimalN));
for k = 1:length(uSAni)
mAni(k, AnimalN == uSAni(k)) = 1;
end
pg = zeros(2,3);
for k = 1:2
    L = [L,sprintf('\n')];
    L = [L , klb{k}];
        ccc = [sq((bprd(:,1,1)- bprd(:,2,1))./sum(bprd(:,:,1),2)),...
            sq((bprd(:,1,2)-bprd(:,2,2))./sum(bprd(:,:,2),2)),...
            sq((bprd(:,1,3)-bprd(:,2,3))./sum(bprd(:,:,3),2))];
        if strcmp(titleS,'traj')
            dcc = [sum(bprd(:,:,1),2), sum(bprd(:,:,2),2), sum(bprd(:,:,3),2)];
        else
            dcc = [sum(bprd(:,:,1),2), sum(bprd(:,:,2),2), sum(bprd(:,:,3),2)];
            figure(3);clf;% sum(sq(<720),2)<1  sum(<1,2)<1
            for nnn = 1:3
                subplot(3,1,nnn)
                plot(dcc(:,nnn)/240, ccc(:,nnn),'.')
                hold on
                plot([1 1]*thrP/240, [-1 1], 'r--')
            end
            title(titleS)
            pause(.01)
        end
        slc = ObjTypes==(k-1) & AnimalN>0 & sum(sum(dcc<thrP,2),2)<1;
        % select sessions according to the Object, animal number and any
        % period included should be larger than the thrP. 
        ccc = ccc(slc,:);
        tmpA = mAni(:,slc);
        % ccc(AnimalN(ObjTypes==(k-1))==0,:)=[];% ccc(sum(~isnan(ccc),2)<1,:)=[];% 
        ccc(isnan(ccc))=0;
        if strcmp(filename(end), 'G') % "Grouped" by animals 
            ccc = bsxfun(@rdivide,tmpA, sum(tmpA,2))*ccc;
        end
        if k ==2
            for kk = 1:3
                % compare real and VR
                [pg(1,kk),~,stat] = ranksum(c(:,kk),ccc(:,kk));
                pg(2,kk) = stat.zval;%ranksum;
            end
        end
        c = [c;ccc];
        for kk = 1:3
            % test DR ? 0
            [p(kk), ~, stat] = signrank(ccc(:,kk));
            if isfield(stat, 'zval')
            zz(k,kk) = stat.zval;
            else
            zz(k,kk) = nan;
            end
            L = [ L , sprintf(' %s: %.5f(z:%.5f,n:%d),', lb{kk}, p(kk), zz(k,kk),length(ccc(:,kk)))];
            pp(k,kk) = p(kk);
        end
        xx = [xx; repmat({'1','2','3'},size(ccc,1),1)];
        x = [x;ones(size(ccc))*k + 0*repmat([-.2 0 .2],size(ccc,1),1)];
        if 0 % nargin<7 test difference between different periods. two sided. 
            if k ==1
                p2(1,2) = signrank(ccc(:,2) - ccc(:,1));
                p2(1,3) = signrank(ccc(:,3) - ccc(:,1));
                p2(2,3) = signrank(ccc(:,2) - ccc(:,3));
            else
                p2(2,1) = signrank(ccc(:,2) - ccc(:,1));
                p2(3,1) = signrank(ccc(:,3) - ccc(:,1));
                p2(3,2) = signrank(ccc(:,2) - ccc(:,3));
            end
        else
            if k ==1
                [p2(1,2), ~, stat] = signrank(ccc(:,2) - ccc(:,1));% , 0, 'tail','right'
                
                if ~isfield(stat, 'zval')
                    stat.zval = nan;
                end
                    z2(1,2) = stat.zval;
                
                [p2(1,3), ~, stat] = signrank(ccc(:,3) - ccc(:,1));% , 0, 'tail','right'
                if ~isfield(stat, 'zval')
                    stat.zval = nan;
                end
                 z2(1,3)= stat.zval;
                [p2(2,3), ~, stat] = signrank(ccc(:,2) - ccc(:,3));% , 0, 'tail','right'
                if ~isfield(stat, 'zval')
                    stat.zval = nan;
                end
                 z2(2,3)= stat.zval;
            else
                [p2(2,1), ~, stat] = signrank(ccc(:,2) - ccc(:,1));% , 0, 'tail','right'
                if ~isfield(stat, 'zval')
                    stat.zval = nan;
                end
                 z2(2,1)= stat.zval;
                [p2(3,1), ~, stat] = signrank(ccc(:,3) - ccc(:,1));% , 0, 'tail','right'
                if ~isfield(stat, 'zval')
                    stat.zval = nan;
                end
                 z2(3,1)= stat.zval;
                [p2(3,2), ~, stat] = signrank(ccc(:,2) - ccc(:,3));% , 0, 'tail','right'
                if ~isfield(stat, 'zval')
                    stat.zval = nan;
                end
                 z2(3,2)= stat.zval;
            end
        end
end
if 1 % plot use gramm (ggplot). 
    clear g
    figure(1)
    clf
    g(1,1)=gramm('x',[x(:,1);x(:,2);x(:,3)],'y',c(:),'color',[xx(:,1);xx(:,2);xx(:,3)]);% ;, 'label',L(:)
    %  g.geom_label('color','k','dodge',0.7,'VerticalAlignment','bottom','HorizontalAlignment','center');
    g(1,2) = copy(g(1,1));
    g(1,1).geom_jitter('width',0.4,'height',0);
    g(1,1).stat_summary('geom',{'bar','black_errorbar'}, 'type', 'sem');% stat_boxplot();
    g(1,2).stat_boxplot();
    g(1,2).stat_violin('fill','transparent');
    g.set_limit_extra([.1,.1], [.05 .05]);
    % g.stat_boxplot('geom',{'bar','black_errorbar'});
    g.set_names('color','Phase','x',[], 'y',[]);
    g.set_title([titleS, ' pval:', L]);
    g.draw();
    pos = get(1,'Position');
    kf = figure(1);
    set(kf, 'Position',  [ pos(1:2)  960   420])
    for k =1:2
        subplot(2,2,2*k)
        plot(repmat(1:3, sum(x(:,1) == k),1)', c(x(:,1) == k, :)', 'LineWidth', 2, 'Color',[.4, .4, .4])
        hold on
        plot(repmat(1:3, sum(x(:,1) == k),1)', c(x(:,1) == k, :)', 'k.', 'LineWidth', 5)
        if k ==1
            xlabel(sprintf('%s 1-2:%.5f, 2-3:%.5f, 1-3:%.5f\n z:1-2:%.5f, 2-3:%.5f, 1-3:%.5f', klb{k}, p2(1,2), p2(2,3), p2(1,3), z2(1,2), z2(2,3), z2(1,3)))
        else
            xlabel(sprintf('%s 1-2:%.5f, 2-3:%.5f, 1-3:%.5f\n z:1-2:%.5f, 2-3:%.5f, 1-3:%.5f', klb{k}, p2(2,1), p2(3,2), p2(3,1), z2(2,1), z2(3,2), z2(3,1)))
            for nn = 1:3
                if pg(nn)<.05
                    plot(nn, 1, 'rx', 'LineWidth', 2)
                    text(nn+.05,.9,[ sprintf(' %.5f', pg(nn))])% '\leftarrow',
                end
            end
        end
        xlim([.5 3.5])
        ylim([-1 1])
        grid on
    end
    print(1, [savepath, filename, '.eps'], '-depsc')
    savefig(1,[savepath, filename, '.fig'])
end
save([savepath, filename,'.p2p.mat'], 'p2')