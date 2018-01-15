% load VRSpatialNovel.Occupancy.3-5.mat
clear g
p2 = cell(2,1);
p2g = cell(2,1);
figure(2)
colormap gray
figure(3)
colormap gray
titleS = {'low';'high';'traj';'Lfacing';'Hfacing';};
pp = zeros(2,2,5,nradium);
pp2 = zeros(2,2,5,nradium);
for k = 1%1:5
    switch k
        case 1
            bprd = lPeriod;
            tk = titleS{k};
        case 2
            bprd = hPeriod;
            tk = titleS{k};
            continue
        case 3
            bprd = TCount*100;
            tk = titleS{k};
            continue
        case 4
            bprd = lwPeriod;
            tk = titleS{k};
        case 5
            bprd = hwPeriod;
            tk = titleS{k};
            continue
    end
    for nn = 3%1:nradium
        Dm =  Dms(nn);%.05+ .04*nn; %
        try
            [p2{1}(k,nn), p2{2}(k,nn), pp(:,:,k,nn)] = StatePlotSN(sq(bprd(:,:,:,nn))/(Dm*10)^2,Names.VRSpatNovelc, sprintf('SN%s at %d cm',tk, fix(Dm*100)), Savedir,['VRSpatNovelExp.statSum.', tk, '.oneside.', num2str( Dm*100)],AnimalN,1);
        catch
            p2{1}(k,nn) = nan;
            p2{2}(k,nn) = nan;
        end
        return
        try
        [p2g{1}(k,nn), p2g{2}(k,nn), pp2(:,:,k,nn)] = StatePlotSN(sq(bprd(:,:,:,nn))/(Dm*10)^2,Names.VRSpatNovelc, sprintf('SN%s at %d cm',tk, fix(Dm*100)), Savedir,['VRSpatNovelExp.statSum.', tk, '.oneside.', num2str( Dm*100), '.G'],AnimalN,1);
        catch
            p2g{1}(k,nn) = nan;
            p2g{2}(k,nn) = nan;
        end
    end
end
return
save('VRSpacialNovelty.p.oneside.mat', 'pp', 'pp2')
save('VRSpacialNovelty.p2p.oneside.mat', 'p2', 'p2g')
% save('VRSpacialNovelty.p.mat', 'pp', 'pp2')
% save('VRSpacialNovelty.p2p.mat', 'p2', 'p2g')
