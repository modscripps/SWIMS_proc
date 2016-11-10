DD='D:\swims\BS03\data_mat\CTD';
IX='D:\swims\BS03\indexes\CTD_BS03_matfiles.mat';
CFs = dir([DD '\CTD*']);

lf = length(CFs);

ydays=[]; prs=[];

for i=ib:ie
    clear SWraw
    load([DD '\' CFs(i).name])
    ydays = [ydays SWraw.SWIMStime];
    prs = [prs SWraw.Pr];
end