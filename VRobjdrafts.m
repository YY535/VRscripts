preVRnew
% VRObj:
for k =1% 1 for VRobj, 2 for VRspatialNov
    Names.([keyWords{k}, 'c']) = true(size(Names.(keyWords{k})));
    for n = 1:length(Names.(keyWords{k}))
        subDir = Names.(keyWords{k}){n};
        try
            xyz = h5read([Workdir, subDir, '/', subDir, '.h5'], '/preprocessed/Rigid Body/Rat/Position');
        catch
            Names.([keyWords{k}, 'c'])(n)= false;
            fprintf('\n%s', subDir)
        end
    end
    fprintf('\n\n the utility for %s is: %.2i', keyWords{k}, sum(Names.([keyWords{k}, 'c']))/length(Names.([keyWords{k}, 'c'])))
end

%
load('/storage/weiwei/data/VRnew/VRobjPos.mat')
aa = ObjPos;
% figure(1);clf;hold on for visualize results
disP = @(x,y)(sqrt(sum(bsxfun(@minus, x, y).^2,2)));
C_R = @(x,A,c)(bsxfun(@minus,x,c)*A');
%  [nfile, nposition, nphase, nradium]
nfile = length(Names.VRObjc);
nposition = 2;
nphase = 3;
Dms = [6:15]/100;%.07+.05*[0:2]; size of the distance. 
nradium = length(Dms);
lPeriod = zeros(nfile, nposition, nphase, nradium);% low period 
hPeriod = zeros(nfile, nposition, nphase, nradium);% high period
lwPeriod = zeros(nfile, nposition, nphase, nradium);% lowwalk period
hwPeriod = zeros(nfile, nposition, nphase, nradium);% highwalk period
TCount = zeros(nfile, nposition, nphase, nradium);% lowwalk period
% we used only the low or low  facing period in figures.  
ObjTypes = zeros(length(Names.VRObjc),1);
AnimalN = zeros(length(Names.VRObjc),1);
figure(224);clf;hold on
figure(225);clf;hold on
collectTraj = cell(2,1);
collectTraj2 = cell(2,1);
isplottraj = false;% true;%
spds = @(x,n)(sqrt(sum(([x((n+1):end,:);zeros(n,size(x,2))]- x).^2,2))/n);
spdst  = @(x,n,t)(sqrt(sum(([x((n+1):end,:);zeros(n,size(x,2))]- x).^2,2))./([t((n+1):end);(t(end)+t(1:n))] - t));
%%
for n = 24:length(Names.VRObjc)% no proper data before session 24
    if ~Names.VRObjc(n) || n==32 || n==71 ||n==75 || n==103 || n==92 || n==113 || n==132
        continue
    end
    %%
    subDir = Names.VRObj{n};
    xyz=h5read([Workdir, subDir, '/', subDir, '.h5'], '/preprocessed/Rigid Body/Rat/Position');
    samplingrate = ceil(1/diff(xyz.Time(1:2)));
    evs.evlog=h5read([Workdir, subDir, '/', subDir, '.h5'], '/events/eventlog');
    evs.evFram=h5read([Workdir, subDir, '/', subDir, '.h5'], '/events/phaseStartFrameNum');
    evs.evName=h5read([Workdir, subDir, '/', subDir, '.h5'], '/events/eventNames');
    evs.evArg=h5read([Workdir, subDir, '/', subDir, '.h5'], '/events/eventArguments');
    evs.ObjType = h5readatt([Workdir, subDir, '/', subDir, '.h5'], '/','VR_OBJECT_TYPE');
    
    if strcmp(evs.ObjType, 'VR')
        ObjTypes(n) = 0;
    else
        ObjTypes(n) = 1;
    end
    AnimalN(n) = find(strcmp(uniRatN, subDir(36:37)));
%     figure(224);clf;
    %% behavior segmentation
    FR = 10; % lowpassed frequency3;%
    xyz.Time(isnan(xyz.X)) = [];
    xyz.Frame(isnan(xyz.X)) = [];
    
    xyz.X = ButterFilter(xyz.X(~isnan(xyz.Z)), FR, samplingrate, 'low');
    xyz.Y = ButterFilter(xyz.Y(~isnan(xyz.Z)), FR, samplingrate, 'low');
    xyz.Z = ButterFilter(xyz.Z(~isnan(xyz.Z)), FR, samplingrate, 'low');
    vz = spdst(xyz.Y,1,xyz.Time);% diff([0;xyz.Y],1,1)./diff([0;xyz.Time],1,1);
    v2 = sqrt(sum(diff([0 0;[xyz.X xyz.Z]],1,1).^2,2))./diff([0;xyz.Time],1,1);
    vxz = sqrt(sum(SmoothSpeed([xyz.X, xyz.Z],  xyz.Time, 40).^2,2));
    %% climb-jump behavior 
    % detect and then exclude. 
    % Detected periods are stored for consistency. 
    % idea: detect beginning by "center-crossing" event, then detect ending
    % after the beginning as point when animal leave the cylinder. 
    
    if 0% ~exist([Savedir, '/', subDir, '.jumpPrd.mat'])
        try
            % jumping climbing periods prd
            tmp_j = conv(1*(vz<-.3), [1 1 1 1]', 'same')>0;
            prd = StartEnding(tmp_j(:));
            stP = [xyz.X(prd(:,1)), xyz.Z(prd(:,1))];
            [tmpD, tmpI]= min([disP(stP, aa(1,[1 3])), disP(stP, aa(2,[1 3]))],[],2);
            zz = tmpD < .2 & tmpD >= .08;
            PosE = prd(zz,:);% & v2(prd(:,1))>.25
            tmpD = tmpD(zz,:);
            tmpI = tmpI(zz,:);
            stP = [xyz.X(prd(zz,2)), xyz.Z(prd(zz,2))];
            tmpD2 = disP(stP, aa(tmpI,[1 3]));
            %% center crossing as the begining
            prdC = [StartEnding([0;disP([xyz.X, xyz.Z], aa(1,[1 3]))<.035]); ...
                StartEnding([0;disP([xyz.X, xyz.Z], aa(2,[1 3]))<.035])]; 
            % start and ending point of entering into the r = 3.5cm areas
            % around the two objects, respectively.  
            prdC(prdC(:,1)<200,:)=[];
            % not at very first when rat just entered into the arena
            cprdC = fix(mean(prdC,2));
            % center of the crossing points. 
            cprdC(cprdC>(length(xyz.X)-300)) = [];
            dists =sqrt( (xyz.X(cprdC-240) - xyz.X(cprdC+240)).^2 + ...
                (xyz.Z(cprdC-240) - xyz.Z(cprdC+240)).^2);
            % also condition on the distance during this "crossing" event. 
            % which should at least longer than 8 cm. 
            % to makesure this is a proper "crossing", instead of a
            % "rearing" when animal is sitting on the top. 
            prdC = prdC(dists >= .08,:);
            cprdC = cprdC(dists >= .08,:);
            prdB = zeros(size(prdC,1),1);
            prdE = zeros(size(prdC,1),1);
            %%
            ee = [];
            for kk = 1:size(prdC,1)
                try
                    if sum(prdC(kk)>prdB & prdC(kk)<prdE)
                        continue
                    end
                    tmp_center = [xyz.X(cprdC(kk)), xyz.Z(cprdC(kk))];
                    [~, tmp_cid] = min(disP(aa(1:2,[1 3]),tmp_center));
                    if (aa(tmp_cid,[1 3]) - [xyz.X(cprdC(kk)-60), xyz.Z(cprdC(kk)-60)])*(aa(tmp_cid,[1 3]) - [xyz.X(cprdC(kk)+60), xyz.Z(cprdC(kk)+60)])'>0
                        continue
                    end
                    try
                        prdB(kk) = prdC(kk,1) - find(vz(prdC(kk,1)-[0:600])>.15 & disP([xyz.X(prdC(kk,1)-[0:600]), xyz.Z(prdC(kk,1)-[0:600])],aa(tmp_cid,[1 3]))>=.05, 1,'first');
                    catch
                        prdB(kk) = prdC(kk,1) - find(vz(prdC(kk,1)-[0:900])>.15 & disP([xyz.X(prdC(kk,1)-[0:900]), xyz.Z(prdC(kk,1)-[0:900])],aa(tmp_cid,[1 3]))>=.05, 1,'first');
                    end
                    PosEid = find(PosE(:,1)>prdB(kk) & tmpD2>.15,1,'first');
                    tmp2 = PosE(PosEid,2);
                    tmp2 = min(tmp2, prdB(kk)+100+find(xyz.Y((prdB(kk)+100):tmp2)<.15, 1,'first'));
                    prdE(kk) = tmp2 + find(vz(tmp2+[1:200])>-.1 ,1,'first');
                    
                    
                    if isplottraj
                        %%
                        figure(224)
                        hold on
                        [traj, A]= CenterRotate([xyz.X(prdB(kk):10:prdE(kk)), xyz.Z(prdB(kk):10:prdE(kk))], aa(tmp_cid,[1 3]), .05); 
                        plot(traj(1,1), traj(1,2),'rx')
                        plot(traj(end,1), traj(end,2),'bx')
                        jp = C_R(tmp_center, A, aa(tmp_cid,[1 3]));
                        plot(jp(1),jp(2),'k.')
                        jp = C_R([xyz.X(PosE(PosEid,1):10:PosE(PosEid,2)), xyz.Z(PosE(PosEid,1):10:PosE(PosEid,2))], A, aa(tmp_cid,[1 3]));
                        plot(traj(:,1), traj(:,2),'Color', [.5 .5 .5])
                        plot(jp(:,1), jp(:,2), 'r')
                    end
                catch
                end
            end
            figure(224);         grid on; axis tight
            %%
            if sum(prdB>0)
                prdE(prdB==0)=[];
                prdB(prdB==0)=[];
                jumpPrd = [xyz.Frame(prdB), xyz.Frame(prdE), sign(xyz.X(prdB))];
                save([Savedir, '/', subDir, '.jumpPrd.mat'],'jumpPrd')
            else
                jumpPrd = [];
            end
        catch
            jumpPrd = [];
        end
    else
        if exist([Savedir, '/', subDir, '.jumpPrd.mat'])
            load([Savedir, '/', subDir, '.jumpPrd.mat'])
        else
            jumpPrd = [];
        end
    end
    %% walking. 
    %
    Periods = zeros(length(evs.evFram),2);
    tmpY = xyz.Y;
    nint = 40;
    spd = [[zeros(nint,1);xyz.X]-[xyz.X;zeros(nint,1)], ...
       [zeros(nint,1);xyz.Z]-[xyz.Z;zeros(nint,1)]];
   spd = spd((nint+1):end,:);
   if subDir(39) == 'L'
       disc = [disP([xyz.X xyz.Z],aa(1,[1 3])), disP([xyz.X xyz.Z],aa(2,[1 3]))];
       headv = [sum(bsxfun(@minus,aa(1,[1 3]), [xyz.X, xyz.Z]).*spd,2)>0, sum(bsxfun(@minus,aa(2,[1 3]), [xyz.X, xyz.Z]).*spd,2)>0];
       
       cpos = aa(1:2,[1 3]);
       % object go first.
   else
       disc = [disP([xyz.X xyz.Z],aa(2,[1 3])), disP([xyz.X xyz.Z],aa(1,[1 3]))];
       headv = [sum(bsxfun(@minus,aa(2,[1 3]), xyz.X).*spd,2)>0, sum(bsxfun(@minus,aa(1,[1 3]), xyz.X).*spd,2)>0];
       cpos = aa([2 1],[1 3]);
   end
    for k = 1:size(jumpPrd,1)
        tmpY(xyz.Frame>=jumpPrd(k,1) & xyz.Frame<=jumpPrd(k,2)) = 3;
        disc(xyz.Frame>=jumpPrd(k,1) & xyz.Frame<=jumpPrd(k,2),:) = 10;
        spd(xyz.Frame>=jumpPrd(k,1) & xyz.Frame<=jumpPrd(k,2)) = -1;
        headv(xyz.Frame>=jumpPrd(k,1) & xyz.Frame<=jumpPrd(k,2),:) = false;
        % exclude jumping periods. 
    end
    for kk = 1:length(evs.evFram)
        nn = find(evs.evlog.Frame == evs.evFram(kk));
        Periods(kk,:) = [evs.evlog.Frame(nn+2), evs.evlog.Frame(nn+3)];
        % I took only 'wait_duration' event for each event. 
        tt = xyz.Frame>Periods(kk,1) & xyz.Frame<=Periods(kk,2);
        tmp_t = xyz.Frame(tt);
        tmp_l = tmpY(tt)<.15 ;%& vxz(tt)>.01;
        tmp_h = tmpY(tt)>=.15;% & tmpY(tt)<=.28;
        
        for nn = 1:nradium % .5 + .5*nn
            tmp_disc = disc(tt,:)<=Dms(nn);% (.05 + .04*nn);
            TCount(n,:,kk, nn) = [size(StartEnding(tmp_disc(:,1)),1), size(StartEnding( tmp_disc(:,2)),1)];%<.15
            
            lPeriod(n,:,kk, nn) = [sumdiff(tmp_t(StartEnding(tmp_l & tmp_disc(:,1))),1,2), sumdiff(tmp_t(StartEnding(tmp_l & tmp_disc(:,2))),1,2)];%<.15
            hPeriod(n,:,kk, nn) = [sumdiff(tmp_t(StartEnding(tmp_h & tmp_disc(:,1))),1,2), sumdiff(tmp_t(StartEnding(tmp_h & tmp_disc(:,2))),1,2)];% .15~.28
            
            lwPeriod(n,:,kk, nn) = [sumdiff(tmp_t(StartEnding(headv(tt)>.03 & tmp_l & tmp_disc(:,1))),1,2), sumdiff(tmp_t(StartEnding(headv(tt)>.03 & tmp_l & tmp_disc(:,2))),1,2)];%<.15
            hwPeriod(n,:,kk, nn) = [sumdiff(tmp_t(StartEnding(headv(tt)>.03 & tmp_h & tmp_disc(:,1))),1,2), sumdiff(tmp_t(StartEnding(headv(tt)>.03 & tmp_h & tmp_disc(:,2))),1,2)];% .15~.28
        end
    end
    
end
%% 
save([Savedir, 'VRObj.Occupancy.v.mat'], 'lPeriod','hPeriod','ObjTypes','TCount', 'lwPeriod','hwPeriod','AnimalN','Dms')
sumPlot
return
