% speed change w.r.t. geodisic distance. 
figure(225);clf
nn = {'obj','sham'};
nbins = 20;%
cdistance = .05;
limT = 4*240;
bg = 4;
for m = 1:2 % two backgrounds
    for n =1:2 % 2 phases
        for k = 1:2 % objects
            if k ==1
                ll = (DefDis(Names.VRSpatNovelc == m, ObjSet, n));
                ll = cell2mat(ll(:));
                ll2 = (dDefDis(Names.VRSpatNovelc == m, ObjSet, n));
                ll2 = cell2mat(ll2(:));
                Dislb = (Deflect(Names.VRSpatNovelc == m, ObjSet, n));
                Dislb = cell2mat(Dislb(:));
            else
                ll = (DefDis(Names.VRSpatNovelc == m, setdiff(1:Onn,ObjSet), n));
                ll = cell2mat(ll(:));
                ll2 = (dDefDis(Names.VRSpatNovelc == m, setdiff(1:Onn,ObjSet), n));
                ll2 = cell2mat(ll2(:));
                Dislb = (Deflect(Names.VRSpatNovelc == m, setdiff(1:Onn,ObjSet), n));
                Dislb = cell2mat(Dislb(:));
            end
            xv = linspace(0,2,nbins);
            crit = Dislb(:,1)<cdistance & Dislb(:,2)>=(pi/2) & Dislb(:,3)<limT;
            h = histc(bsxfun(@rdivide, ll(crit,:), ll(crit,1)),xv);% (abs(Dislb(:,2))>(pi/2),:)
            h = h(:,bg:end);
            hy = linspace(0,1,nDefDis);
            hy = hy(bg:end);
            subplot(4,4,(m-1)*4 + (n-1)*2 +k)
            imagesc(xv,hy, bsxfun(@rdivide,h,sum(h,1))')
            title(sprintf('%s in phase.%d', nn{k}, n))
            xv2 = linspace(-.015,.015,nbins);
            h2 = histc(ll2(crit,:),xv2);% (abs(Dislb(:,2))>(pi/2),:)
            subplot(4,4,8+(m-1)*4 + (n-1)*2 +k)
            imagesc(xv2,linspace(0,1,nDefDis), bsxfun(@rdivide,h2,sum(h2,1))')
            title(sprintf('%s in phase.%d', nn{k}, n))
            hold on; plot( [0 0],[0 1],'r:')
            plot(median(ll2(crit,:)),linspace(0,1,nDefDis),'k')
            
        end
    end
end
subplot(4,4,1);ylabel('VR')
subplot(4,4,5);ylabel('Black')
subplot(4,4,9);ylabel('VR_d')
subplot(4,4,13);ylabel('Black_d')
xlabel('pass')
colormap default
% savefig(225, 'VRSpatialNovel.Fix.traj.pass.fig')

figure(226);clf
nn = {'obj','sham'};
for m = 1:2 % two backgrounds
    for n =1:2 % 2 phases
        for k = 1:2 % objects
            if k ==1
                ll = (DefDis(Names.VRSpatNovelc == m, ObjSet, n));
                ll = cell2mat(ll(:));
                ll2 = (dDefDis(Names.VRSpatNovelc == m, ObjSet, n));
                ll2 = cell2mat(ll2(:));
                Dislb = (Deflect(Names.VRSpatNovelc == m, ObjSet, n));
                Dislb = cell2mat(Dislb(:));
            else
                ll = (DefDis(Names.VRSpatNovelc == m, setdiff(1:Onn,ObjSet), n));
                ll = cell2mat(ll(:));
                ll2 = (dDefDis(Names.VRSpatNovelc == m, setdiff(1:Onn,ObjSet), n));
                ll2 = cell2mat(ll2(:));
                Dislb = (Deflect(Names.VRSpatNovelc == m, setdiff(1:Onn,ObjSet), n));
                Dislb = cell2mat(Dislb(:));
            end
            crit = Dislb(:,1)<cdistance  & Dislb(:,3)<limT;%& Dislb(:,2)<(pi/2)
            h = histc(bsxfun(@rdivide, ll(crit,:), ll(crit,1)),xv);% (abs(Dislb(:,2))>(pi/2),:)
            h = h(:,bg:end);
            hy = linspace(0,1,nDefDis);
            hy = hy(bg:end);
            subplot(4,4,(m-1)*4 + (n-1)*2 +k)
            imagesc(xv,hy, bsxfun(@rdivide,h,sum(h,1))')
            title(sprintf('%s in phase.%d', nn{k}, n))
            h2 = histc(ll2(crit,:),xv2);% (abs(Dislb(:,2))>(pi/2),:)
            subplot(4,4,8+(m-1)*4 + (n-1)*2 +k)
            imagesc(xv2,linspace(0,1,nDefDis), bsxfun(@rdivide,h2,sum(h2,1))')
            title(sprintf('%s in phase.%d', nn{k}, n))
            hold on; plot([0 0],[0 1], 'r:')
            plot(median(ll2(crit,:)),linspace(0,1,nDefDis),'k')
           
        end
    end
end
subplot(4,4,1);ylabel('VR')
subplot(4,4,5);ylabel('Black')
subplot(4,4,9);ylabel('VR_d')
subplot(4,4,13);ylabel('Black_d')
xlabel('close')
colormap default
% savefig(226, 'VRSpatialNovel.Fix.traj.close.fig')




figure(215);clf
% nn = {'obj','sham'};
% nbins = 20;%
% cdistance = .04;
% limT = 600;
for m = 1:2 % two backgrounds
    for n =1:2 % 2 phases
        for k = 1:2 % objects
            if k ==1
                ll = (hDefDis(Names.VRSpatNovelc == m, ObjSet, n));
                ll = cell2mat(ll(:));
                ll2 = (hdDefDis(Names.VRSpatNovelc == m, ObjSet, n));
                ll2 = cell2mat(ll2(:));
                Dislb = (Deflect(Names.VRSpatNovelc == m, ObjSet, n));
                Dislb = cell2mat(Dislb(:));
            else
                ll = (hDefDis(Names.VRSpatNovelc == m, setdiff(1:Onn,ObjSet), n));
                ll = cell2mat(ll(:));
                ll2 = (hdDefDis(Names.VRSpatNovelc == m, setdiff(1:Onn,ObjSet), n));
                ll2 = cell2mat(ll2(:));
                Dislb = (Deflect(Names.VRSpatNovelc == m, setdiff(1:Onn,ObjSet), n));
                Dislb = cell2mat(Dislb(:));
            end
            crit = Dislb(:,1)<cdistance & Dislb(:,2)>=(pi/2) & Dislb(:,3)<limT;
            h = histc(bsxfun(@rdivide, ll(crit,:), ll(crit,1)),xv);% (abs(Dislb(:,2))>(pi/2),:)
            h = h(:,bg:end);
            hy = linspace(0,1,nDefDis);
            hy = hy(bg:end);
            subplot(4,4,(m-1)*4 + (n-1)*2 +k)
            imagesc(xv,hy, bsxfun(@rdivide,h,sum(h,1))')
            title(sprintf('%s in phase.%d', nn{k}, n))
            %         ll =
%             xv2 = linspace(0,.3,nbins);
            xv2 = linspace(-.015,.015,nbins);
            h2 = histc(ll2(crit,:),xv2);% (abs(Dislb(:,2))>(pi/2),:)
%             [h2, xv2] = hist(ll2(Dislb(:,1)<.03 & Dislb(:,2)>=(pi/2),:),nbins);% (abs(Dislb(:,2))>(pi/2),:)
            subplot(4,4,8+(m-1)*4 + (n-1)*2 +k)
            imagesc(xv2,linspace(0,1,nDefDis), bsxfun(@rdivide,h2,sum(h2,1))')
            title(sprintf('%s in phase.%d', nn{k}, n))
            hold on; plot( [0 0],[0 1],'r:')
            plot(median(ll2(crit,:)),linspace(0,1,nDefDis),'k')
            
        end
    end
end
subplot(4,4,1);ylabel('VR')
subplot(4,4,5);ylabel('Black')
subplot(4,4,9);ylabel('VR_d')
subplot(4,4,13);ylabel('Black_d')
xlabel('pass')
colormap default
% savefig(215, 'VRSpatialNovel.Fix.traj.pass.ph.fig')

figure(216);clf
nn = {'obj','sham'};
for m = 1:2 % two backgrounds
    for n =1:2 % 2 phases
        for k = 1:2 % objects
            if k ==1
                ll = (hDefDis(Names.VRSpatNovelc == m, ObjSet, n));
                ll = cell2mat(ll(:));
                ll2 = (hdDefDis(Names.VRSpatNovelc == m, ObjSet, n));
                ll2 = cell2mat(ll2(:));
                Dislb = (Deflect(Names.VRSpatNovelc == m, ObjSet, n));
                Dislb = cell2mat(Dislb(:));
            else
                ll = (hDefDis(Names.VRSpatNovelc == m, setdiff(1:Onn,ObjSet), n));
                ll = cell2mat(ll(:));
                ll2 = (hdDefDis(Names.VRSpatNovelc == m, setdiff(1:Onn,ObjSet), n));
                ll2 = cell2mat(ll2(:));
                Dislb = (Deflect(Names.VRSpatNovelc == m, setdiff(1:Onn,ObjSet), n));
                Dislb = cell2mat(Dislb(:));
            end
            crit = Dislb(:,1)<cdistance & Dislb(:,3)<limT; % & Dislb(:,2)<(pi/2)
            h = histc(bsxfun(@rdivide, ll(crit,:), ll(crit,1)),xv);% (abs(Dislb(:,2))>(pi/2),:)
            h = h(:,bg:end);
            hy = linspace(0,1,nDefDis);
            hy = hy(bg:end);
            subplot(4,4,(m-1)*4 + (n-1)*2 +k)
            imagesc(xv,hy, bsxfun(@rdivide,h,sum(h,1))')
            title(sprintf('%s in phase.%d', nn{k}, n))
            %         ll =
            h2 = histc(ll2(crit,:),xv2);% (abs(Dislb(:,2))>(pi/2),:)
%             [h2, xv2] = hist(ll2(Dislb(:,1)<.03 & Dislb(:,2)<(pi/2),:),nbins);% (abs(Dislb(:,2))>(pi/2),:)
            subplot(4,4,8+(m-1)*4 + (n-1)*2 +k)
            imagesc(xv2,linspace(0,1,nDefDis), bsxfun(@rdivide,h2,sum(h2,1))')
            title(sprintf('%s in phase.%d', nn{k}, n))
            hold on; plot([0 0],[0 1], 'r:')
            plot(median(ll2(crit,:)),linspace(0,1,nDefDis),'k')
           
        end
    end
end
subplot(4,4,1);ylabel('VR')
subplot(4,4,5);ylabel('Black')
subplot(4,4,9);ylabel('VR_d')
subplot(4,4,13);ylabel('Black_d')
xlabel('close')
colormap default
% savefig(216, 'VRSpatialNovel.Fix.traj.close.ph.fig')
return
%% 
% assign labels
SpeedCmp = [];
lbDef = cell(size(Deflect));
for k = 1:size(Deflect,1)
    for n = 1:length(Deflect(k,:))
    lbDef{k,n} = ones(size(Deflect{k,n},1),1)*k;
    end
end
%
bgn = 1;%[];
SpeedCmpS = cell(2,1);
for phs = 1:2
for k = 1:Onn
    if isempty(bgn)
        tmp = [cell2mat(Deflect(:,k,phs)), cell2mat(lbDef(:,k,phs))];%
    else
        tmp = [cell2mat(Deflect(Names.VRSpatNovelc==bgn,k,phs)), cell2mat(lbDef(Names.VRSpatNovelc==bgn,k,phs))];%
    end
    tmp(tmp(:,1)>.05 |tmp(:,3)>limT ,:)=[];% tmp(:,1)>.06|tmp(:,2)<(pi/2)
    SpeedCmp = [SpeedCmp; [tmp(:,[4 5]), ones(size(tmp,1),1)*k,tmp(:,6) ]];
end
token = SpeedCmp(:,3)==1;
snum = [accumarray(SpeedCmp(token,4),1, [nfile,1]), accumarray(SpeedCmp(~token,4),1, [nfile,1])];
SpeedCmpS{phs} = [accumarray(SpeedCmp(token,4), diff(SpeedCmp(token,1:2),1,2), [nfile,1]), accumarray(SpeedCmp(~token,4),diff(SpeedCmp(~token,1:2),1,2), [nfile,1])]./snum; 

%%
end
%
figure(221);clf
clear g
g = gramm('x',reshape(repmat({'obj', 'sham'}, 2*nfile,1),[],1), 'y', reshape(cell2mat(SpeedCmpS),[],1), 'color', reshape([ones(nfile,2);2*ones(nfile,2)],[],1));
g.geom_jitter('width',0.4,'height',0);
g.stat_violin('fill','transparent');%
g.stat_boxplot('width',.2);
[pp(1,1), ~, stat] = signrank(diff(SpeedCmpS{1},1,2));
zz(1,1) = stat.zval;
n(1,1) = sum(~isnan(diff(SpeedCmpS{1},1,2)));
[pp(1,2), ~, stat] = signrank(diff(SpeedCmpS{2},1,2));
zz(1,2) = stat.zval;
n(1,2) = sum(~isnan(diff(SpeedCmpS{2},1,2)));

[pp(2,1), ~, stat] = signrank(SpeedCmpS{1}(:,1));
zz(2,1) = stat.zval;
n(2,1) = sum(~isnan(SpeedCmpS{1}(:,1)));
[pp(2,2), ~, stat] = signrank(SpeedCmpS{1}(:,2));
zz(2,2) = stat.zval;
n(2,2) = sum(~isnan(SpeedCmpS{1}(:,2)));
[pp(3,1), ~, stat] = signrank(SpeedCmpS{2}(:,1));
zz(3,1) = stat.zval;
n(3,1) = sum(~isnan(SpeedCmpS{2}(:,1)));
[pp(3,2), ~, stat] = signrank(SpeedCmpS{2}(:,2));
zz(3,2) = stat.zval;
n(3,2) = sum(~isnan(SpeedCmpS{2}(:,2)));
g.set_title(sprintf('Speed Compare: obj vs sham: ph1: %.2f, ph2: %.2f\n ph1: obj: %.2f, sham: %.2f\n ph2: obj: %.2f, sham: %.2f', signrank(diff(SpeedCmpS{1},1,2)),signrank(diff(SpeedCmpS{2},1,2)), ...
    signrank(SpeedCmpS{1}(:,1)),  signrank(SpeedCmpS{1}(:,2)), ...
    signrank(SpeedCmpS{2}(:,1)),  signrank(SpeedCmpS{2}(:,2))))
g.set_names('x', [],'y', 'speed difference (m/s)', 'color', 'phase')
g.draw()
drawnow
figure(221);
grid on
return
%% check the time spending 
figure(220);clf
clear g
phn = {'ph.1', 'ph.2'};
objn = {'Obj', 'sham'};
for kk = 1:3:13
    spdtime = cell(2,2);
    for k = 1:2 % environment
        for n = 1:2 % phase
            spdtime{k,n} = [sq(lPeriod(Names.VRSpatNovelc==k, 1, n, kk)), ...
                sq(mean(lPeriod(Names.VRSpatNovelc==k, 2:end, n, kk),2))];
            spdtime{k,n} = [spdtime{k,n}, k*ones(size(spdtime{k,n},1),1), ...
                n*ones(size(spdtime{k,n},1),1)];
            
        end
    end
    %  i take only3:3:15 cm
    spdtime = cell2mat(spdtime(:));
    spdtime(spdtime(:,1) ==0 |spdtime(:,2) ==0,:)=[];
    for n = 1:2
        g((kk-1)/3 +1,n) = gramm('x', phn(spdtime(:,4)), 'y', spdtime(:,n)/240, 'color', spdtime(:,3));
        g((kk-1)/3 +1,n).stat_violin('fill','transparent');
        g((kk-1)/3 +1,n).stat_boxplot('width',.2);
        g((kk-1)/3 +1,n).set_title(sprintf('%d cm %s', fix(Dms(kk)*100), objn{n}));
    end
end
g.set_names('color','env', 'x',[], 'y','time(sec)');
g.draw();
%%
% close 220
meantrajs  = zeros(20, nfile, Onn, 2);
startSpeed = cell(nfile,Onn,2);
tmp(:,1)>.05 |tmp(:,3)>limT
for k =1:nfile
    for n = 1:Onn
        for m = 1:2
            if isempty(DefDis{k,n,m})
                continue
            end
            tmp = Deflect{k,n,m}(:,1)<.05 & Deflect{k,n,m}(:,3)<limT;
            if isempty(tmp)
                continue
            end
            startSpeed{k,n,m} = DefDis{k,n,m}(tmp,1);
            meantrajs(:,k,n,m) = mean(bsxfun(@rdivide, DefDis{k,n,m}(tmp,:), startSpeed{k,n,m}),1);
        end
    end
end
%
x = linspace(0, 1, 20);
yo = sq(meantrajs(:,:,:,2));
% yo(:,:,2) = mean(tmp(:,:,tmp()),3);
y = cell(nfile,2);
c = cell(nfile,2);
% figure
emptymats = false(nfile,1);
for k = 1:nfile
    if isempty(DefDis{k,1,2})
        emptymats(k) = true;
        continue
    end
    y{k,1} = yo(:,k,1);
    y{k,2} = mean(yo(:,k,sq(yo(1,k,:)>0)),3);
    c{k,1} = 'obj';
    c{k,2} = 'sham';
end
y = y(Names.VRSpatNovelc==1 & ~emptymats,:);
c = c(Names.VRSpatNovelc==1 & ~emptymats,:);
g = gramm('x', x, 'y', y(:), 'color', c(:))
g.stat_summary('type','std');
g.set_title('speed compare');
figure('Position',[100 100 800 550]);
g.draw();
%%
c1 = cell2mat(startSpeed(Names.VRSpatNovelc==1,1,2));
c2 = cell2mat(reshape(startSpeed(Names.VRSpatNovelc==1,2:end,2),[],1));
clear g
a = ['obj', 'shame']
x = [repmat(1,length(c2),1);repmat(2, length(c2),1)];
y = [c1;c2]
g = gramm('x', x, 'y', y)
g.geom_jitter();
g.stat_boxplot()
g.set_title(sprintf('init speed compare p = %.3f', ranksum(c1,c2)));
figure
g.draw()