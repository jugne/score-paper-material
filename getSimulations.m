clear

currdir = cd;

workingdir = {'FullDatesInference'};

% delete all xmls in all working dirs

for i = 1 : length(workingdir)
    convertXmlToSimulation(workingdir{i}, {'Anseriformes', 'Galliformes'})
    cd(currdir)
    convertXmlToSimulation_rootFilteringMax(workingdir{i}, {'Anseriformes', 'Galliformes'})
    cd(currdir)
    convertXmlToSimulation_segmentRootFilteringMax(workingdir{i}, {'Anseriformes', 'Galliformes'})
    cd(currdir)
end