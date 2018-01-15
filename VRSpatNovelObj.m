% VRSpatNovel VR object comparison 
preVRnew
% VRObj:
k =2;
Names.([keyWords{k}, 'c']) = ones(size(Names.(keyWords{k})));
for n = 1:length(Names.(keyWords{k}))
    subDir = Names.(keyWords{k}){n};
    addnum = 0;
    if strcmp(subDir(16:20), 'Black')
        Names.([keyWords{k}, 'c'])(n)= 2;
        addnum = 11;
    end
    if subDir(22)=='6'
        if str2num(subDir([24:25]+addnum))<27
            Names.([keyWords{k}, 'c'])(n)= false;
            continue
        end
    end
    if ~strcmp(subDir([42:43]+addnum), '3D')
            Names.([keyWords{k}, 'c'])(n)= false;
            continue
    end
    try
        xyz = h5read([Workdir, subDir, '/', subDir, '.h5'], '/preprocessed/Rigid Body/Rat/Position');
    catch
        Names.([keyWords{k}, 'c'])(n)= 0;
        fprintf('\n%s', subDir)
    end
end
fprintf('\n\n the prs for %s is: %.2i', keyWords{k}, sum(Names.([keyWords{k}, 'c'])>0)/length(Names.([keyWords{k}, 'c'])))

%%
F_Pos = h5readatt([Workdir, subDir, '/', subDir, '.h5'], '/', 'VR_SPATIAL_FIXED_OBJECT_POSITIONS');
N_Pos = h5readatt([Workdir, subDir, '/', subDir, '.h5'], '/', 'VR_SPATIAL_NOVELTY_OBJECT_POSITIONS');
disP = @(x,y)(sqrt(sum(bsxfun(@minus, x, y).^2,2)));
%  [nfile, nposition, nphase, nradium]
nfile = length(Names.VRSpatNovelc);
Onn = 403;
nposition = Onn;
nphase = 4;
Dms = .05;%[3:15]/100;%.07+.05*[0:2];
nradium = length(Dms);
lPeriod = zeros(nfile, nposition, nphase, nradium);% low period 
hPeriod = zeros(nfile, nposition, nphase, nradium);% high period
lwPeriod = zeros(nfile, nposition, nphase, nradium);% lowwalk period
hwPeriod = zeros(nfile, nposition, nphase, nradium);% highwalk period
TCount = zeros(nfile, nposition, nphase);% lowwalk period
AnimalN = zeros(length(Names.VRSpatNovelc),1);
figure(224); clf;caxis([0 .3])
Deflect = cell(nfile,nposition,2);
NF_Pos = [-.3 0 .3;0 0 0; 0 0 0];
ObjTypes = zeros(nfile,1);
figure(224)
nDefDis = 20;
DefDis = cell(nfile,Onn,2);
dDefDis = cell(nfile,Onn,2);
hDefDis = cell(nfile,Onn,2);
hdDefDis = cell(nfile,Onn,2);
ObjSet = 1;%:2;
ObjPositions = cell(nfile,1);
%%
for n = 26:nfile %length(Names.VRSpatNovelc)% 24
    if ~Names.VRSpatNovelc(n) || n==43 || n==91 || n==26
        continue
    elseif Names.VRSpatNovelc(n)==2
        addnum = 11;
        ObjTypes = 1;
    else
        addnum = 0;
        ObjTypes = 0;
    end
    %% data loading
    clear evs xyz
    subDir = Names.VRSpatNovel{n};
    xyz=h5read([Workdir, subDir, '/', subDir, '.h5'], '/preprocessed/Rigid Body/Rat/Position');
    samplingrate = ceil(1/diff(xyz.Time(1:2)));
    evs.evlog=h5read([Workdir, subDir, '/', subDir, '.h5'], '/events/eventlog');
    evs.evFram=h5read([Workdir, subDir, '/', subDir, '.h5'], '/events/phaseStartFrameNum');
    evs.evName=h5read([Workdir, subDir, '/', subDir, '.h5'], '/events/eventNames');
    evs.evArg=h5read([Workdir, subDir, '/', subDir, '.h5'], '/events/eventArguments');
    
    AnimalN(n) = find(strcmp(uniRatN, subDir([39:40]+addnum)));
    ObjN = [str2num(subDir(45 +addnum)), str2num(subDir(47 +addnum)),str2num(subDir(49 +addnum))];
    ObjP = zeros(Onn,2);
    ObjP(1,:) = F_Pos([1 3],  [ObjN(3)])';
    if length(ObjSet)>1
        ObjP(2,:) = N_Pos([1 3],  [ObjN(1)])';
    end
    ObjAVO = [N_Pos([1 3],  ObjN(1)), F_Pos([1 3],  ObjN(3))]';
    randcan = bsxfun(@times,rand(Onn*100,2)-.5, 2*[35 35])/100;% [25 5]
    randcan(abs(randcan(:,2))>.1,:)=[];
    %   find enough random samples.
    dd = @(x,y)(sqrt(sum(bsxfun(@minus, x, y).^2,2)));
    drcan = find(dd(randcan, ObjAVO(1,:))>=.1 & dd(randcan, ObjAVO(2,:))>=.1);% , Onn-length(ObjSet),'first'.2
    if length(drcan)< (Onn-length(ObjSet))
        drcan = find(dd(randcan, ObjAVO(1,:))>=.15 & dd(randcan, ObjAVO(2,:))>=.15);
    end
    if length(drcan) > (Onn-length(ObjSet))
        [~, dummy] = sort(randcan(drcan,1));
        dummy = dummy(fix(linspace(1,length(drcan),(Onn-length(ObjSet)))));
        drcan = drcan(dummy);
    end
    
    ObjP((length(ObjSet)+1):end,:) = randcan(drcan,:);
    ObjPositions{n} = ObjP;
    FR = 10; %3;%  lowpassed frequency
    xyz.Time(isnan(xyz.X)) = [];
    xyz.Frame(isnan(xyz.X)) = [];
    
    xyz.X = ButterFilter(xyz.X(~isnan(xyz.Z)), FR, samplingrate, 'low');
    xyz.Y = ButterFilter(xyz.Y(~isnan(xyz.Z)), FR, samplingrate, 'low');
    xyz.Z = ButterFilter(xyz.Z(~isnan(xyz.Z)), FR, samplingrate, 'low');
    
    %% walking. 
    
    lprd = find(diff(evs.evlog.Frame,1,1)>10000);
    Periods = zeros(length(lprd),2);
    tmpY = xyz.Y;
    
    tmp_out = abs(xyz.X-F_Pos(1,2))>.5 | abs(xyz.Z-F_Pos(3,2))>.25;
    nint = 40;
    spd = SmoothSpeed([xyz.X, xyz.Z], 1, 40);
    v = sqrt(sum(SmoothSpeed([xyz.X, xyz.Z], xyz.Time, 40).^2,2));
    dv = SmoothSpeed(v, 1, 2);
    disc = zeros(length(tmpY),size(ObjP,1));
    headv = zeros(length(tmpY),size(ObjP,1));
    headang = zeros(length(tmpY),size(ObjP,1));
    for k = 1:size(ObjP,1)
        disc(:,k) = disP([xyz.X xyz.Z],ObjP(k,:));
        headv(:,k) = sum(bsxfun(@minus,ObjP(k,:), [xyz.X, xyz.Z]).*spd,2)>0;
        headang(:,k) = angle(sqrt(-1)*(-xyz.Z+ObjP(k,2))+(-xyz.X+ObjP(k,1)));
    end
    disc(tmp_out,:) = 10;% exclude all at once
        fprintf('\n nfile:%d', n)
    for kk = 1:2%length(lprd)
        Periods(kk,:) = [evs.evlog.Frame(lprd(kk)), evs.evlog.Frame(lprd(kk)+1)];
        if kk == 1
            Periods(kk,1) = Periods(kk,1) + fix((Periods(kk,2) - Periods(kk,1))/3*2);% took the last 2 min in phase 1
        end
        % I took only 'wait_duration' event for each event. 
        tt = xyz.Frame>max(Periods(kk,1)-10,1) & xyz.Frame<=Periods(kk,2);
        tmp_t = xyz.Frame(tt);
        tmp_v = v(tt);
        tmp_dv = dv(tt);
        tmp_l = tmpY(tt)<.15;
        tmp_h = tmpY(tt)>=.15 ;% & tmpY(tt)<=.28;
        tmp_x = xyz.X(tt);
        tmp_z = xyz.Z(tt);
        % about trajectory
        tmp_disc1 = disc(tt,:)<=.1;%(.03 + .05
        disc_tmp = disc(tt,:);
        ang_tmp = headang(tt,:);
        
        for k =1:Onn
        if 0
            trajs = StartEnding(tmp_disc1(:,k));
            if ~isempty(trajs)% && Names.VRSpatNovelc(n)==1
                trajs((trajs(:,2)- trajs(:,1))<60,:) = [];
                ntraj = size(trajs,1);
                Deflect{n,k,kk} = zeros(ntraj, 5);% [closes_distance, angle]
                DefDis{n,k,kk} = zeros(ntraj, nDefDis);
                dDefDis{n,k,kk} = zeros(ntraj, nDefDis);
                hDefDis{n,k,kk} = zeros(ntraj, nDefDis);
                hdDefDis{n,k,kk} = zeros(ntraj, nDefDis);
                for nn = 1:ntraj
                    [dummy, dummyid] = min(disc_tmp(trajs(nn,1): trajs(nn,2),k));
                   
                    Deflect{n,k,kk}(nn,:) = [dummy, ...
                        angle(exp(sqrt(-1)*(2*pi- ang_tmp(trajs(nn,2),k) + ang_tmp(trajs(nn,1),k)))), ...
                        (trajs(nn,2)- trajs(nn,1)), ...
                        mean(tmp_v(max(trajs(nn,1)-5,1):min(trajs(nn,1),length(tmp_v)))), ...( trajs(nn,1)+dummyid-1)
                        mean(tmp_v(max(trajs(nn,1)+dummyid,1):min(trajs(nn,1)+dummyid+5,length(tmp_v))))];%end<[maxtime, ang_diff, duration]
                    DefDis{n,k,kk}(nn,:) = tmp_v(floor(linspace(trajs(nn,1), trajs(nn,2), nDefDis)));% [min(disc_tmp(trajs(nn,1): trajs(nn,2),k)),angle(exp(sqrt(-1)*(2*pi- ang_tmp(trajs(nn,2),k) + ang_tmp(trajs(nn,1),k))))];%<
                    dDefDis{n,k,kk}(nn,:) = tmp_dv(floor(linspace(trajs(nn,1), trajs(nn,2), nDefDis)));% [min(disc_tmp(trajs(nn,1): trajs(nn,2),k)),angle(exp(sqrt(-1)*(2*pi- ang_tmp(trajs(nn,2),k) + ang_tmp(trajs(nn,1),k))))];%<
                    hDefDis{n,k,kk}(nn,:) = tmp_v(floor(linspace(trajs(nn,1), trajs(nn,1)+dummyid-1, nDefDis)));% [min(disc_tmp(trajs(nn,1): trajs(nn,2),k)),angle(exp(sqrt(-1)*(2*pi- ang_tmp(trajs(nn,2),k) + ang_tmp(trajs(nn,1),k))))];%<
                    hdDefDis{n,k,kk}(nn,:) = tmp_dv(floor(linspace(trajs(nn,1), trajs(nn,1)+dummyid-1, nDefDis)));% [min(disc_tmp(trajs(nn,1): trajs(nn,2),k)),angle(exp(sqrt(-1)*(2*pi- ang_tmp(trajs(nn,2),k) + ang_tmp(trajs(nn,1),k))))];%<
                    
                end
                TCount(n,k,kk) = ntraj;%<.15
            else
                TCount(n,k,kk) = 0;%<.15
            end
        end
            for nn = 1:nradium % .5 + .5*nn
                tmp_disc = disc(tt,:)<=Dms(nn);%(.03 + .05*(nn-1));
                lPeriod(n,k,kk, nn) = sumdiff(tmp_t(StartEnding(tmp_v>.03 &tmp_l & tmp_disc(:,k))),1,2);%<.15
                hPeriod(n,k,kk, nn) = sumdiff(tmp_t(StartEnding(tmp_h & tmp_disc(:,k))),1,2);% .15~.28
                
                lwPeriod(n,k,kk, nn) = sumdiff(tmp_t(StartEnding(tmp_v>.03 & headv(tt)>.03 & tmp_l & tmp_disc(:,k))),1,2);%<.15
                hwPeriod(n,k,kk, nn) = sumdiff(tmp_t(StartEnding(headv(tt)>.03 & tmp_h & tmp_disc(:,k))),1,2);% .15~.28
            end
        end
%         pause
%         figure(224);clf
    end
end
fprintf('\n finish \n')
% VRSpatialNoveltyTrajV
% return
if 0
%% statistics of deflection. 
% I'd need to compare the statistics from sham location vs the statistics
% from the acture object location. 
% here I could compute the comparison of the obj/sham and distorted
% obj/sham 
limT = 4*240;
limD = .03;
for bgs = 1:2
figure(218+bgs);clf
for phs = 1:2;%1;
    %
AngComp = cell(2,1);
token = Names.VRSpatNovelc==bgs;
AngComp{1} = cell2mat(reshape(sq(Deflect(token,ObjSet,phs)),[],1));
AngComp{1}(AngComp{1}(:,1)>limD |AngComp{1}(:,3)>limT ,:) = [];% 
AngComp{1} = [AngComp{1}, ones(size(AngComp{1},1),1)];

% AngComp{2} = cell2mat(reshape(sq(Deflect(:,:,1)),[],1));
AngComp{2} = cell2mat(reshape(sq(Deflect(token,setdiff(1:Onn,ObjSet),phs)),[],1));
AngComp{2}(AngComp{2}(:,1)>limD |AngComp{2}(:,3)>limT,:) = [];
AngComp{2} = [AngComp{2}, 2*ones(size(AngComp{2},1),1)];
sizeAngComp = [size(AngComp{1},1), size(AngComp{2},1)];
sizeAngComp(3) = sum(sizeAngComp(1:2))
AngComp = cell2mat(AngComp);
nshuff = 1000;
DefRatio = zeros(nshuff+1,1);
Rconst = sizeAngComp(2)/sizeAngComp(1);
AngComp(:,2) = abs(AngComp(:,2));
for k = 1:nshuff
    [~,rids] = sort(rand(sizeAngComp(3),1));
    DefRatio(k) = Rconst*sum(AngComp(rids(1:sizeAngComp(1)), 2)<(pi/2))/ ...
        sum(AngComp(rids((1+sizeAngComp(1)):end), 2)<(pi/2));
%     DefRatio(k) = sum(AngComp(rids(1:sizeAngComp(1)), 2)<(pi/2));
end
DefRatio(end)=Rconst*sum(AngComp(1:sizeAngComp(1), 2)<(pi/2))/...
    sum(AngComp((1+sizeAngComp(1)):end, 2)<(pi/2));
DefRatio(isinf(DefRatio))=0;
EmpiricalPval = sum(DefRatio(1:(end-1))>DefRatio(end))/nshuff
figure(218+bgs);
subplot(1,2,phs)
hold on
for k = 1:2
    if k ==1
        tokens =cell2mat(reshape(sq(Deflect(token,ObjSet,phs)),[],1));
        tokens(tokens(:,3)>limT ,:) = [];% tokens(:,1)>.03 |
    else
        tokens =cell2mat(reshape(sq(Deflect(token,setdiff(1:Onn,ObjSet),phs)),[],1));
        tokens(tokens(:,3)>limT ,:) = [];% tokens(:,1)>.03 |
    end
    plot(abs(tokens(:,2)), tokens(:,1),'o')
end
set(gca,'XTick', [0:.25:1]*pi, 'XTicklabel', [0:.25:1]*180)
xlabel('theta (degree)')
ylabel('closest distance (m)')
grid on; axis tight
legend('Obj', 'sham')
title(sprintf('Phase %d: eP-Val is %.3f', phs, EmpiricalPval))
end
end
end
return

%%
save(sprintf('%sVRSpatialNovelObj.Obj.%d-%d.mat',Savedir, Dms([1 3])*100),'Dms','Deflect','DefDis', 'lPeriod','hPeriod','TCount', 'lwPeriod','hwPeriod','AnimalN')

return
%%
spdist = zeros(nradium,2,2);% nrad, 2obj, phases, 2 behavior
 figure(2);clf;hold on; 
 
Sspdist = zeros(nradium,2,2);% nrad, 2obj, phases, 2 behavior
ranksumSig = @(x,y)((median(x)-median(y)));
 for n = 1:nradium
     for m =1:2
     spdist(n,:,m,1) = [ranksum(sq(lPeriod(Names.VRSpatNovelc==1, 1,m,n)), sq(lPeriod(Names.VRSpatNovelc==2, 1,m,n))), ...
         ranksum(sq(sum(lPeriod(Names.VRSpatNovelc==1, 2:3,m,n),2)), sq(sum(lPeriod(Names.VRSpatNovelc==2, 2:3,m,n),2)))];
     spdist(n,:,m,2) = [ranksum(sq(lwPeriod(Names.VRSpatNovelc==1, 1,m,n)), sq(lwPeriod(Names.VRSpatNovelc==2, 1,m,n))), ...
         ranksum(sq(sum(lwPeriod(Names.VRSpatNovelc==1, 2:3,m,n),2)), sq(sum(lwPeriod(Names.VRSpatNovelc==2, 2:3,m,n),2)))];
     Sspdist(n,:,m,1) = [ranksumSig(sq(lPeriod(Names.VRSpatNovelc==1, 1,m,n)), sq(lPeriod(Names.VRSpatNovelc==2, 1,m,n))), ...
         ranksumSig(sq(sum(lPeriod(Names.VRSpatNovelc==1, 2:3,m,n),2)), sq(sum(lPeriod(Names.VRSpatNovelc==2, 2:3,m,n),2)))];
     Sspdist(n,:,m,2) = [ranksumSig(sq(lwPeriod(Names.VRSpatNovelc==1, 1,m,n)), sq(lwPeriod(Names.VRSpatNovelc==2, 1,m,n))), ...
         ranksumSig(sq(sum(lwPeriod(Names.VRSpatNovelc==1, 2:3,m,n),2)), sq(sum(lwPeriod(Names.VRSpatNovelc==2, 2:3,m,n),2)))];
     end
 for k = 1:2;
     subplot(nradium,2,n*2-1);
     plot(sq(sum(lPeriod(Names.VRSpatNovelc==k, :,1,n),2)), sq(sum(lPeriod(Names.VRSpatNovelc==k, :,2,n),2)), '.');hold on; 
     grid on; axis tight
%      xlabel('phase 1')
     ylabel([num2str(Dms(n)*100),'cm'])% 'phase 2:', 
%      title('low for all ')
     subplot(nradium,2,n*2);
     plot(sq(sum(lwPeriod(Names.VRSpatNovelc==k, :,1,n),2)), sq(sum(lwPeriod(Names.VRSpatNovelc==k, :,2,n),2)), '.');hold on; 
     grid on; axis tight
%      xlabel('phase 1')
%      ylabel('phase 2')
%      title('lowFacing for all')
     
 end
 end
 for k = 1:2;
     figure(20+k);
     for n =1:2;
         subplot(2,2,n);
         imagesc(Dms*100, [], sq(spdist(:,:,n,k))', [0 .1]);
         xlabel(['phs.', num2str(n)]);
         set(gca,'YTick',1:2,'YTicklabel', {'obj','sham'})
         subplot(2,2,2+n);
         imagesc(Dms*100, [], sq(Sspdist(:,:,n,k))', [-.1 .1]);
         xlabel(['phs.', num2str(n)]);
         ylabel('V-B')
         set(gca,'YTick',1:2,'YTicklabel', {'obj','sham'})
     end
     colormap jet
 end
 savefig(2,'VRSN.explorationtime.F_Pos.VBcmp.fig')
 savefig(21,'VRSN.explorationtime.F_Pos.VBcmp.lPeriod.fig')
 savefig(22,'VRSN.explorationtime.F_Pos.VBcmp.lwPeriod.fig')
print(21,'VRSN.explorationtime.F_Pos.VBcmp.lPeriod.jpg', '-djpeg')
print(22,'VRSN.explorationtime.F_Pos.VBcmp.lwPeriod.jpg', '-djpeg')
 

return
%%
figure(1);clf
clear g
totime = sum(sq(sum(lPeriod(:,[1 2 4:end],1:2)>0,3)),2);
secs = Names.VRSpatNovelc==1 & totime >0;
x = reshape(repmat([1 2],sum(secs),1),[],1);
yo = sq(mean(lPeriod(secs,1:2,1:2),2));
ys = sq(mean(lPeriod(secs,4:end,1:2),2));
y = reshape((yo-ys)./(yo+ys),[],1);
% x = x(yo(:)>0 & ys(:)>0);
% y = y(yo(:)>0 & ys(:)>0);
p = zeros(2,1);
zs = zeros(2,1);
for k = 1:2
    [p(k), ~, stat] = signrank(y(x==k));% , 0, 'tail','right'
    zs(k)=stat.zval;
end
g = gramm('x',x(:),'y',y(:))
g.geom_jitter('width',0.4,'height',0);
g.stat_summary('geom',{'bar','black_errorbar'}, 'type', 'sem');% stat_boxplot();
    g.set_title(sprintf(' ph1:pval:%.5f, zval:%.5f\nph2:pval:%.5f, zval:%.5f', p(1),zs(1),p(2),zs(2)));
    g.draw();
    figure(1)
    set(gca,'XTick',1:2)
    grid on
    savefig(1,'VRspatialnovelty.fig')
     [p,~,stat]=signrank(y(1:54)-y(55:end))
% 
% p =
% 
%     0.0183
% 
% 
% stat = 
% 
%           zval: -2.3593
%     signedrank: 449
