function data = getData()

tab = readtable('data/reliability_exp1.csv');

names = table2cell(unique(tab(:,'subject_name')));

data = zeros(0,5);
for i = 1:length(names)
    name = names{i};
    tabSub = tab(strcmp(tab.subject_name,name),:);
    tabSub = tabSub(strcmp(tabSub.block_type,'Test'),:);
    tabSub = tabSub(strcmp(tabSub.task,'B'),:);
    reliabilities = unique(tabSub.stim_reliability);
    rels = zeros(size(tabSub.stim_reliability));
    for iRel = 1:length(reliabilities)
        rels(tabSub.stim_reliability==reliabilities(iRel)) = iRel;
    end
    data = cat(1,data,[repmat(i,height(tabSub),1),rels,...
        tabSub.stim_orientation,tabSub.resp_category,tabSub.stim_category]);
end