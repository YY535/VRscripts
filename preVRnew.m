
Workdir ='/storage2/nickdg/data/VR_Experiments_Round_2/processed_data/';%
%'/storage/share/Projects/VRExperiments_Round2/Data/motion_tracking/processed/';% '/storage/share/Projects/VRExperiments_Round2/Data/motion_tracking/processed/';
% '/storage/share/Projects/VRExperiments_Round2/Data/motion_tracking/processed/';
% '/storage/nickdg/data/VR_Experiments_Round_2/processed_data/';
cd(Workdir)
keyWords = {'VRObj', 'VRSpatNovel'};
Names = [];
% expdir = {'/storage/nickdg/data/VR_Experiments_Round_2/processed_data_by_experiment/VRObjectExp/';...
%     '/storage/nickdg/data/VR_Experiments_Round_2/processed_data_by_experiment/VRSpatNovelExp/'};
for kk = 1:2
%     cd(expdir{kk})
a = dir([keyWords{kk}, '*']);
gg = false(size(a));
for k =1:length(a)
    gg(k)=exist([Workdir, a(k).name,'/', a(k).name, '.h5'],'file');
end
a = a(gg);
FileB = cell(size(a));
for k =1:length(a)
    FileB{k} = a(k).name;
end
Names.(keyWords{kk})=FileB;
end
Savedir = '/storage/weiwei/data/VRnew/PreProcess/';
cd(Savedir)
%%
RatNames = cell(length(Names.VRObj),1);
for k = 1:length(Names.VRObj)
    RatNames{k} = Names.VRObj{k}(36:37);
end
uniRatN = unique(RatNames);
