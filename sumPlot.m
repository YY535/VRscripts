% load VRObj.Occupancy.mat
nradium = length(Dms);
clear g
p2 = cell(5,nradium);
p2g = cell(5,nradium);
pp = zeros(2,3,5,nradium);
ppg = zeros(2,3,5,nradium);
pg = zeros(3,5,nradium);
pgg = zeros(3,5,nradium);
figure(2);clf
colormap gray

k = 4;
bprd = lwPeriod;
tk = 'Lfacing';

% k = 1;
% bprd = lPeriod;
% tk = 'low';

nn = 7; % select which radium
Dm =  Dms(nn);% radium
[p2{k,nn}, pg(:,k,nn), pp(:,:,k,nn)] = StatePlot(sq(bprd(:,:,:,nn)),ObjTypes, sprintf('%s at %d cm',tk, fix(Dm*100)), Savedir,['VRobj.statSum.', tk, '.', num2str( Dm*100)],AnimalN,1);
[p2g{k,nn}, pgg(:,k,nn), ppg(:,:,k,nn)] = StatePlot(sq(bprd(:,:,:,nn)),ObjTypes, sprintf('%s at %d cm',tk, fix(Dm*100)), Savedir,['VRobj.statSum.', tk, '.', num2str( Dm*100), '.G'],AnimalN,1);