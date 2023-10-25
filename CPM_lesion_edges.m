function CPM_lesion_edges(all_mats,atlas,network1,network2,prefix,dataset)

% lesion edges from a 3D set of matrices for a given network-pair, or select those edges only
% possible network names are on lines 20-21 (Shen atlas) and 25-26 (Schaefer300-7Network atlas)
% if network1 is same as network2, then that network's nodes alone will be targeted 
% INPUTS:
% 1. FC matrices (edges x edges x subjects/trials) - load these in
% 2. atlas (1=Shen268, 2=Schaefer300 - 7 network version)
% 3. network1= name of network 1 of network-network lesion pair
% 4. network2= name of network 2 of network-network lesion pair
% 5. prefix= prefix for output file names
% 6. name of dataset folder name
% OUTPUTS:
% saves two output matrices (lesioned and network-only versions) in Group/Lesioned_mats

global globalDataDir;
datapath=[globalDataDir filesep dataset];
if atlas==1
atlas_name='Shen268';
net_labels=importdata('Lookup_shen268');
net_labels=net_labels(:,7);
net_names={'SMN';'[]'; 'CO'; 'AUD'; 'DMN';'[]'; 'VIS'; 'FPN';...
    'SAL'; 'SUB'; 'VAN'; 'DAN'; 'Unknown'};
elseif atlas==2
net_labels=importdata(['Schaefer_label300_7networks.mat']);
net_names={'DMN'; 'DAN'; 'SAL'; 'FPCN'; 'LIM'; 'VIS'; 'SMN'};
atlas_name='Schaefer300';
end

% get network numbers
net1_num=strmatch(network1,net_names,'exact');
net2_num=strmatch(network2,net_names,'exact');

% loop through matrices
all_mats_lesioned=[]; all_mats_networksOnly_allsubs=[];
for i=1:size(all_mats,3)
    curr_mat=[]; curr_mat_lesioned=[]; networksOnly_mat=[];
    % lesion network-pair
    curr_mat=all_mats(:,:,i); networksOnly_mat=NaN(size(curr_mat,1),size(curr_mat,2));
    curr_mat_lesioned=curr_mat;
    curr_mat_lesioned(net_labels==net1_num,net_labels==net2_num)=NaN;
    curr_mat_lesioned(net_labels==net2_num,net_labels==net1_num)=NaN;
    all_mats_lesioned=cat(3,all_mats_lesioned,curr_mat_lesioned);
    % select network-pair only (lesion everything else)
    networksOnly_mat(net_labels==net1_num,net_labels==net2_num)=curr_mat(net_labels==net1_num,net_labels==net2_num);
    networksOnly(net_labels==net2_num,net_labels==net1_num)=curr_mat(net_labels==net2_num,net_labels==net1_num);
    all_mats_networksOnly_allsubs=cat(3,all_mats_networksOnly_allsubs,networksOnly_mat);
end

% save
mkdir([datapath filesep 'Group']);
mkdir([datapath filesep 'Group' filesep 'Lesioned_mats']);
outpath=[datapath filesep 'Group' filesep 'Lesioned_mats' filesep];
%cd([datapath filesep 'Group' filesep 'Lesioned_mats']);
if atlas==1
save([outpath 'Shen268_' prefix '_' network1 '-' network2 '_lesioned'],'all_mats_lesioned');
save([outpath 'Shen268_' prefix '_' network1 '-' network2 '_only'],'all_mats_networksOnly_allsubs');
elseif atlas==2
save([outpath 'Schaefer300_' prefix '_' network1 '-' network2 '_lesioned'],'all_mats_lesioned');
save([outpath 'Schaefer300_' prefix '_' network1 '-' network2 '_only'],'all_mats_networksOnly_allsubs');
end
display(['done - outputs are in ' outpath]); 


