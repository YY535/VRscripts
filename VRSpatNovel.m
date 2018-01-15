% VRSpatNovel
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
nposition = 3;
nphase = 4;
Dms = [3:15]/100;%.07+.05*[0:2];
nradium = length(Dms);
lPeriod = zeros(nfile, nposition, nphase, nradium);% low period 
hPeriod = zeros(nfile, nposition, nphase, nradium);% high period
lwPeriod = zeros(nfile, nposition, nphase, nradium);% lowwalk period
hwPeriod = zeros(nfile, nposition, nphase, nradium);% highwalk period
TCount = zeros(nfile, nposition, nphase, nradium);% lowwalk period
AnimalN = zeros(length(Names.VRSpatNovelc),1);
%% 
% period extraction, similar to VRobj
for n = 26:nfile %length(Names.VRSpatNovelc)% 24
    if ~Names.VRSpatNovelc(n) || n==43 || n==91 || n==26
        continue
    elseif Names.VRSpatNovelc(n)==2
        addnum = 11;
    else
        addnum = 0;
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
    % init new fix
    ObjP = [N_Pos(1,ObjN(1)), N_Pos(3,ObjN(1)); ...
        N_Pos(1,ObjN(2)), N_Pos(3,ObjN(2));...
        F_Pos(1,ObjN(3)), F_Pos(3,ObjN(3))];
    
    FR = 10; % 3;% lowpassed frequency
    xyz.Time(isnan(xyz.X)) = [];
    xyz.Frame(isnan(xyz.X)) = [];
    
    xyz.X = ButterFilter(xyz.X(~isnan(xyz.Z)), FR, samplingrate, 'low');
    xyz.Y = ButterFilter(xyz.Y(~isnan(xyz.Z)), FR, samplingrate, 'low');
    xyz.Z = ButterFilter(xyz.Z(~isnan(xyz.Z)), FR, samplingrate, 'low');
    
    %% walking. 
    %
    
    lprd = find(diff(evs.evlog.Frame,1,1)>10000);
    Periods = zeros(length(lprd),2);
    tmpY = xyz.Y;
    
    tmp_out = abs(xyz.X-F_Pos(3,2))>.5 | abs(xyz.Z-F_Pos(3,2))>.25;
    nint = 40;
    spd = SmoothSpeed([xyz.X, xyz.Z], 1, 40);
    v = sqrt(sum(SmoothSpeed([xyz.X, xyz.Z.^2], xyz.Time, 40).^2,2));
    disc = zeros(length(tmpY),size(ObjP,1));
    headv = zeros(length(tmpY),size(ObjP,1));
    for k = 1:size(ObjP,1)
        disc(:,k) = disP([xyz.X xyz.Z],ObjP(k,:));
        headv(:,k) = sum(bsxfun(@minus,ObjP(k,:), [xyz.X, xyz.Z]).*spd,2)>0;
    end
    disc(tmp_out,:) = 10;% exclude all at once
      
    for kk = 1:length(lprd)
        Periods(kk,:) = [evs.evlog.Frame(lprd(kk)), evs.evlog.Frame(lprd(kk)+1)];
        if kk == 1
            Periods(kk,1) = Periods(kk,1) + fix((Periods(kk,2) - Periods(kk,1))/3*2);% took the last 2 min in phase 1
        end
        % I took only 'wait_duration' event for each event. 
        tt = xyz.Frame>Periods(kk,1) & xyz.Frame<=Periods(kk,2);
        tmp_t = xyz.Frame(tt);
        tmp_v = v(tt);
        tmp_l = tmpY(tt)<.15;
        tmp_h = tmpY(tt)>=.15 ;% & tmpY(tt)<=.28;
        for nn = 1:nradium % .5 + .5*nn
            tmp_disc = disc(tt,:)<=Dms(nn);%(.03 + .05*(nn-1));
            for k =1:3
                TCount(n,k,kk, nn) = size(StartEnding(tmp_disc(:,k)),1);%<.15
                lPeriod(n,k,kk, nn) = sumdiff(tmp_t(StartEnding(tmp_v>.005 &tmp_l & tmp_disc(:,k))),1,2);%<.15
                hPeriod(n,k,kk, nn) = sumdiff(tmp_t(StartEnding(tmp_h & tmp_disc(:,k))),1,2);% .15~.28
                
                lwPeriod(n,k,kk, nn) = sumdiff(tmp_t(StartEnding(tmp_v>.005 & headv(tt)>.03 & tmp_l & tmp_disc(:,k))),1,2);%<.15
                hwPeriod(n,k,kk, nn) = sumdiff(tmp_t(StartEnding(headv(tt)>.03 & tmp_h & tmp_disc(:,k))),1,2);% .15~.28
            end
        end
    end
end
%%
save(sprintf('%sVRSpatialNovel.Occupancy.%d-%d.mat',Savedir, Dms([1 3])*100),'Dms', 'lPeriod','hPeriod','TCount', 'lwPeriod','hwPeriod','AnimalN')
sumPlotSN
return
