%     
% virus = 'h5n1';
% directory = 'h5n1_no_reject';
% template = 'template_no_rejection.xml';
% n_runs = 5;
% nrsequences = 200;
% from = 1996;
% to = 2018;
% temperature = 0.01;
% seeds = [505231315, 188090728, 66658260, 32382753, 894204728];

function [] = getXMLallSegmentsFullDates(virus, directory, fasta_file, template, runs, nrsequences, from, to, temperature, seeds)
for run = 1 : runs
    rng(seeds(run));
    mkdir (sprintf('%s',['FullDatesInference/' directory]),sprintf('%s',['run_' num2str(run) '/xmls_coupled']));
    
    % get all fasta files
    data_dir = ['/Users/jugne/Documents/SCORE-paper/StructCoalRe-Dataset_/h5n1/' fasta_file '/run_' num2str(run) '/'];
    segs_files = dir([data_dir '*.fasta']);
    % segments = cell('ha','na');

    seqs = cell(0,0);
    seq_seqs = cell(0,0);
    segments = cell(0,0);
    for i = 1 : length(segs_files)
        fasta = fastaread([data_dir segs_files(i).name]);
        seq_seqs{i} = cell(0,0);
        tmp = strsplit(segs_files(i).name, '_');
        tmp = strsplit(tmp{3}, '.');
        segments{i} = tmp{1};
        for j = 1 : length(fasta)
            seqs{end+1} = fasta(j).Header;
            seq_seqs{i}{j} = fasta(j).Header;
        end
    end

    % get all unique sequences
    unique_seqs = unique(seqs);
    % get the frequency of all sequences
    for i = length(unique_seqs):-1:1
        if sum(ismember(seqs, unique_seqs{i}))<8
            % check if the year is within the limits
            tmp = strsplit(unique_seqs{i}, '|');
            tmp2 = strsplit(tmp{2}, '-');
            if str2double(tmp2{1})<from || str2double(tmp2{1})>to
                unique_seqs(i) = [];
            end              
        end
    end

%     % delete any sequence for which the data isn't specified to the day
%     year = zeros(0,0);
%     for i = length(unique_seqs) : -1 : 1
%         tmp = strsplit(unique_seqs{i}, '|');
%         tmp = strsplit(tmp{2}, '-');
%         year(i) = str2double(tmp{1});
%         if length(tmp) < 3
%             unique_seqs(i) = [];   
%             year(i) = [];
%         end        
%     end
%     % check if the dimension of the year vector is correct
%     if length(year)~=length(unique_seqs)
%         error('error in the definition of the number of years');
%     end
    % get the weights of each sequence as the inverse number of samples per
    % year if the year is within from to
%     weights = zeros(length(year),1);
%     for i = 1 : length(weights)
%         if year(i) >= from && year(i) <= to
%             weights(i) = 1/(sum(year==year(i)));
%         end
%     end
%     
%     use_seqs = zeros(min(sum(weights>0), nrsequences),1);
%     for j = 1 : length(use_seqs)    
%         add_seq = randsample(length(unique_seqs), 1, true, weights);
%         year(add_seq) = -1;
%         
%         % recompute weights
%         weights = zeros(length(year),1);
%         for i = 1 : length(weights)
%             if year(i) >= from && year(i) <= to
%                 weights(i) = 1/(sum(year==year(i)));
%             end
%         end
%         use_seqs(j) = add_seq;
%     end
%     
%     use_seqs = sort(use_seqs);

    use_seqs = zeros(nrsequences,1);
    for j = 1:length(use_seqs)
        use_seqs(j) = j;
    end

    use_segs = segments;

    for r = 0 : 2
        f = fopen(template);
        g = fopen(['FullDatesInference/' directory '/run_' num2str(run) '/xmls_coupled/' virus '_rep' num2str(r) '.xml'], 'w');
        
        est_tip_time = cell(0,0);
        while ~feof(f)
            line = fgets(f);
            if contains(line, 'insert_types')
                for j = 1 : length(use_seqs)
                    type = strsplit(unique_seqs{use_seqs(j)}, '|');
                    type = type{3};
                    fprintf(g, '\t\t\t%s=%s', unique_seqs{use_seqs(j)}, type);
                    if j < length(use_seqs)
                        fprintf(g, ',\n');
                    else
                        fprintf(g, '\n');
                    end
                end
            elseif contains(line, 'insert_alignment')
                for i = 1 : length(segments)
                    fprintf(g, '\t<data id="%s">\n',segments{i});
                    fasta = fastaread([data_dir segs_files(i).name]);

                    seq_length(i) = length(fasta(1).Sequence);

                    for j = 1 : length(use_seqs)
%                         disp(unique_seqs{use_seqs(j)})
%                         if i==2
%                             fds
%                         end
                        ind = find(ismember(seq_seqs{i}, unique_seqs{use_seqs(j)}));
                        if isempty(ind)
                            disp(unique_seqs{use_seqs(j)})
                            disp(segments{i})
                        end
                        ind = ind(1);
                        if (segments{i} == "ha")
                            fprintf(g, '\t\t<sequence id="%s.%s" taxon="%s" totalcount="4" value="%s"/>\n',...
                                segments{i}, unique_seqs{use_seqs(j)},...
                                unique_seqs{use_seqs(j)}, strrep(fasta(ind).Sequence(1:1700), 'U','T'));
                        else
                            fprintf(g, '\t\t<sequence id="%s.%s" taxon="%s" totalcount="4" value="%s"/>\n',...
                                segments{i}, unique_seqs{use_seqs(j)},...
                                unique_seqs{use_seqs(j)}, strrep(fasta(ind).Sequence(1:1350), 'U','T'));
                        end
                            
                    end
                    fprintf(g, '\t</data>\n');

                end
            elseif contains(line, 'insert_run_header')
%                  fprintf(g, '\t\t<run spec="MCMC" chainLength="150000000" storeEvery="1000" preBurnin="100000">\n');            

                 fprintf(g, '\t\t<run id="mcmc" spec="beast.coupledMCMC.CoupledMCMC" chainLength="150000000" storeEvery="1000" deltaTemperature="%.4f" chains="2" resampleEvery="1000" preBurnin="100000">\n', temperature);            

            elseif contains(line, 'insert_taxa')
                for j = 1 : length(use_seqs)
                    fprintf(g, '\t\t\t<taxon spec="Taxon" id="%s"/>\n', unique_seqs{use_seqs(j)});
                end
            elseif contains(line, 'insert_sampling_times')
                for j = 1 : length(use_seqs)
                    time = strsplit(unique_seqs{use_seqs(j)}, '|');
                    tmp = strsplit(time{4}, '-');
                    % put the sampling time mid month or mid day
                    if length(tmp) == 2
                        error('sampling data no accurate enough');
                    end
                    if ~isempty(strfind(time{4}, 'XX-XX'))
                        est_tip_time{end+1} = unique_seqs{use_seqs(j)};
                        time{4} = strrep(time{4}, 'XX-XX', '07-01');
                    elseif ~isempty(strfind(time{4}, 'XX'))
                        est_tip_time{end+1} = unique_seqs{use_seqs(j)};
                        time{4} = strrep(time{4}, 'XX', '15');
                    end
                        
                    deztime = (datenum(time{4},'yyyy-mm-dd')- datenum(tmp{1},'yyyy'))...
                        /(datenum(num2str(str2double(tmp{1})+1),'yyyy')-datenum(tmp{1},'yyyy'))...
                        +str2double(tmp{1});
                    dezstring = sprintf('%.8f', deztime);
                    fprintf(g, '\t\t\t\t%s=',unique_seqs{use_seqs(j)});
                    if j < length(use_seqs)
                        fprintf(g, '%s,\n', dezstring);
                    else
                        fprintf(g, '%s\n', dezstring);
                    end
                end  
            elseif contains(line, 'insert_nr_segments')
                fprintf(g, strrep(line, 'insert_nr_segments', num2str(length(use_segs))));
            elseif contains(line, '<parameter id="clockRate.c" name="stateNode">')
                 if ~isempty(est_tip_time)
                     fprintf(g, strrep(line, 'name="stateNode">', 'lower="0.0015" name="stateNode">'));
                 else
                     fprintf(g, line);
                 end
            elseif contains(line, 'insert_parameters')
                if ~isempty(est_tip_time)
                    fprintf(g, '\t\t\t\t\t\t<parameter id="dateOffset" name="stateNode">0.0</parameter>\n');
                end

                for s = 1 : length(use_segs)
                    fprintf(g, '\t\t\t\t\t\t<parameter id="kappa.s:%s_1" lower="0.0" name="stateNode">%f</parameter>\n',use_segs{s}, exprnd(1));
                    fprintf(g, '\t\t\t\t\t\t<parameter id="kappa.s:%s_3" lower="0.0" name="stateNode">%f</parameter>\n',use_segs{s}, exprnd(1));
                    fprintf(g, '\t\t\t\t\t\t<parameter id="mutationRate.s:%s_1" name="stateNode">1</parameter>\n',use_segs{s});
                    fprintf(g, '\t\t\t\t\t\t<parameter id="mutationRate.s:%s_3" name="stateNode">1</parameter>\n',use_segs{s});
                    fprintf(g, '\t\t\t\t\t\t<parameter id="gammaShape.s:%s_1" name="stateNode">%f</parameter>\n',use_segs{s}, exprnd(1));
                    fprintf(g, '\t\t\t\t\t\t<parameter id="gammaShape.s:%s_3" name="stateNode">%f</parameter>\n',use_segs{s}, exprnd(1));
                    fprintf(g, '\t\t\t\t\t\t<parameter id="freqParameter.s:%s_1" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>\n',use_segs{s});
                    fprintf(g, '\t\t\t\t\t\t<parameter id="freqParameter.s:%s_3" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>\n',use_segs{s});
                end
            elseif contains(line, '<init spec="SegmentTreeInitializer"')
                if ~isempty(est_tip_time)
                    fprintf(g, '\t\t\t\t<init spec="coalre.util.DateOffsetInitializer" dateOffset="@dateOffset">\n');
                    for i = 1 : length(use_segs)
                        fprintf(g, '\t\t\t\t\t<segmentTree idref="%s.tree"/>\n', use_segs{i}); 
                    end
%                     fprintf(g, '\t\t\t\t\t<segmentTree idref="dummy.tree"/>\n');
                    fprintf(g, '\t\t\t\t</init>\n');
                end
                fprintf(g, line)
            elseif contains(line, 'insert_priors')
                % insert sampling time priors
                if ~isempty(est_tip_time)
                    for s = 1 : length(est_tip_time)
                        fprintf(g, '\t\t\t\t\t\t\t\t<distribution id="tipprior.%s" spec="coalre.distribution.TipPrior" dateOffset="@dateOffset" network="@network">\n', est_tip_time{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t<taxonset id="tip.%s" spec="TaxonSet">\n', est_tip_time{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t\t<taxon idref="%s"/>\n', est_tip_time{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t</taxonset>\n');
                        tmp = strsplit(est_tip_time{s}, '|');
                        tmp_ = strsplit(tmp{4}, '-');
                        if ~isempty(strfind(tmp{4}, 'XX-XX'))
%                         fprintf(g, '\t\t\t\t\t\t\t\t\t<distr spec="coalre.util.WeightedSumDistribution">\n');
                        fprintf(g, '\t\t\t\t\t\t\t\t\t\t<Uniform id="Unform.%s" name="distr" lower="%s" upper="%d"/>\n', est_tip_time{s}, tmp_{1}, str2double(tmp_{1})+1);
%                         fprintf(g, '\t\t\t\t\t\t\t\t\t\t<Uniform id="Unform.%s.2" name="distr" lower="1950" upper="1970"/>\n', est_tip_time{s});
%                         fprintf(g, '\t\t\t\t\t\t\t\t\t\t<parameter id="weights.s:%s_3" lower="0.0" name="weights" upper="1.0">0.9 0.05</parameter>\n', est_tip_time{s});
%                         fprintf(g, '\t\t\t\t\t\t\t\t\t</distr>\n');
                        
                        else
                            low = strrep(tmp{4}, 'XX', '01');
                            lower = (datenum(low,'yyyy-mm-dd')- datenum(tmp_{1},'yyyy'))...
                                        /(datenum(num2str(str2double(tmp_{1})+1),'yyyy')-datenum(tmp_{1},'yyyy'))...
                                        +str2double(tmp_{1});
                                    
                            if ~isempty(strfind(tmp_{2}, '04')) || ~isempty(strfind(tmp_{2}, '06')) || ~isempty(strfind(tmp_{2}, '09')) || ~isempty(strfind(tmp_{2}, '11'))
                                up = strrep(tmp{4}, 'XX', '30');
                                upper = (datenum(up,'yyyy-mm-dd')- datenum(tmp_{1},'yyyy'))...
                                        /(datenum(num2str(str2double(tmp_{1})+1),'yyyy')-datenum(tmp_{1},'yyyy'))...
                                        +str2double(tmp_{1});
                            elseif ~isempty(strfind(tmp_{2}, '02')) && (~isempty(strfind(tmp_{1}, '1996')) || ~isempty(strfind(tmp_{1}, '2000')) || ~isempty(strfind(tmp_{1}, '2004')) || ~isempty(strfind(tmp_{1}, '2008')) || ~isempty(strfind(tmp_{1}, '2012')) || ~isempty(strfind(tmp_{1}, '2016')))
                                up = strrep(tmp{4}, 'XX', '29');
                                upper = (datenum(up,'yyyy-mm-dd')- datenum(tmp_{1},'yyyy'))...
                                        /(datenum(num2str(str2double(tmp_{1})+1),'yyyy')-datenum(tmp_{1},'yyyy'))...
                                        +str2double(tmp_{1});
                            elseif ~isempty(strfind(tmp_{2}, '02')) 
                                up = strrep(tmp{4}, 'XX', '28');
                                upper = (datenum(up,'yyyy-mm-dd')- datenum(tmp_{1},'yyyy'))...
                                        /(datenum(num2str(str2double(tmp_{1})+1),'yyyy')-datenum(tmp_{1},'yyyy'))...
                                        +str2double(tmp_{1});
                            else 
                                up = strrep(tmp{4}, 'XX', '31');
                                upper = (datenum(up,'yyyy-mm-dd')- datenum(tmp_{1},'yyyy'))...
                                        /(datenum(num2str(str2double(tmp_{1})+1),'yyyy')-datenum(tmp_{1},'yyyy'))...
                                        +str2double(tmp_{1});
                            end
                            
                                %                         fprintf(g, '\t\t\t\t\t\t\t\t\t<distr spec="coalre.util.WeightedSumDistribution">\n');
                        fprintf(g, '\t\t\t\t\t\t\t\t\t\t<Uniform id="Unform.%s" name="distr" lower="%d" upper="%d"/>\n', est_tip_time{s}, lower, upper);
%                         fprintf(g, '\t\t\t\t\t\t\t\t\t\t<Uniform id="Unform.%s.2" name="distr" lower="1950" upper="1970"/>\n', est_tip_time{s});
%                         fprintf(g, '\t\t\t\t\t\t\t\t\t\t<parameter id="weights.s:%s_3" lower="0.0" name="weights" upper="1.0">0.9 0.05</parameter>\n', est_tip_time{s});
%                         fprintf(g, '\t\t\t\t\t\t\t\t\t</distr>\n');
                           
                        end
                            

                        fprintf(g, '\t\t\t\t\t\t\t\t</distribution>\n');                    
                    end
                end
                for s = 1 : length(use_segs)

                    fprintf(g, '\t\t\t\t\t\t\t\t<prior id="KappaPrior.s:%s_1" name="distribution" x="@kappa.s:%s_1">\n', use_segs{s}, use_segs{s});
                    fprintf(g, '\t\t\t\t\t\t\t\t\t<LogNormal id="LogNormalDistributionModel.%s.1_1" name="distr" M="1.0" S="1.25"/>\n', use_segs{s});
                    fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');
                    fprintf(g, '\t\t\t\t\t\t\t\t<prior id="KappaPrior.s:%s_3" name="distribution" x="@kappa.s:%s_3">\n', use_segs{s}, use_segs{s});
                    fprintf(g, '\t\t\t\t\t\t\t\t\t<LogNormal id="LogNormalDistributionModel.%s.1_3" name="distr" M="1.0" S="1.25"/>\n', use_segs{s});
                    fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');
                    fprintf(g, '\t\t\t\t\t\t\t\t<prior id="GammaPrior.s:%s_1" name="distribution" x="@gammaShape.s:%s_1">\n', use_segs{s}, use_segs{s});
                    fprintf(g, '\t\t\t\t\t\t\t\t\t<Exponential id="ExponentialDistribution.%s_3" name="distr"/>\n', use_segs{s});
                    fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');
                    fprintf(g, '\t\t\t\t\t\t\t\t<prior id="GammaPrior.s:%s_3" name="distribution" x="@gammaShape.s:%s_3">\n', use_segs{s}, use_segs{s});
                    fprintf(g, '\t\t\t\t\t\t\t\t\t<Exponential id="ExponentialDistribution.%s_30" name="distr"/>\n', use_segs{s});
                    fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');

                end
            elseif contains(line, 'insert_operators')
                 
                if ~isempty(est_tip_time)
                    for s = 1 : length(est_tip_time)
                        fprintf(g, '\t\t\t\t<operator spec="MultiTipDatesRandomWalker" useGaussian="true" windowSize="0.1" weight="5">\n');
                        fprintf(g, '\t\t\t\t\t<taxonset idref="tip.%s"/>\n', est_tip_time{s});
                        fprintf(g, '\t\t\t\t\t<tree idref="%s.tree"/>\n', use_segs{1});                 
                        for i = 1 : length(use_segs)
                            fprintf(g, '\t\t\t\t\t<trees idref="%s.tree"/>\n', use_segs{i});                 
                        end
                        fprintf(g, '\t\t\t\t</operator>\n');                    
                    end
                end

                
                for s = 1 : length(use_segs)
                    fprintf(g, '\t\t\t\t<operator id="KappaScaler.s:%s_1" spec="ScaleOperator" parameter="@kappa.s:%s_1" scaleFactor="0.5" weight="0.1"/>\n', use_segs{s}, use_segs{s});
                    fprintf(g, '\t\t\t\t<operator id="KappaScaler.s:%s_3" spec="ScaleOperator" parameter="@kappa.s:%s_3" scaleFactor="0.5" weight="0.1"/>\n', use_segs{s}, use_segs{s});
                    fprintf(g, '\t\t\t\t<operator id="FrequenciesExchanger.s:%s_1" spec="DeltaExchangeOperator" delta="0.01" weight="0.1" parameter="@freqParameter.s:%s_1"/>\n', use_segs{s}, use_segs{s});
                    fprintf(g, '\t\t\t\t<operator id="FrequenciesExchanger.s:%s_3" spec="DeltaExchangeOperator" delta="0.01" weight="0.1" parameter="@freqParameter.s:%s_3"/>\n', use_segs{s}, use_segs{s});
                    fprintf(g, '\t\t\t\t<operator id="alpha_scaler_1.%s" spec="ScaleOperator" parameter="@gammaShape.s:%s_1" scaleFactor="0.75" weight="0.1"/>\n', use_segs{s}, use_segs{s});
                    fprintf(g, '\t\t\t\t<operator id="alpha_scaler_3.%s" spec="ScaleOperator" parameter="@gammaShape.s:%s_3" scaleFactor="0.75" weight="0.1"/>\n', use_segs{s}, use_segs{s});
                end
             elseif contains(line, 'insert_mut_par')
                for s = 1 : length(use_segs)                   
                    fprintf(g, '\t\t\t\t\t<parameter idref="mutationRate.s:%s_1"/>\n', use_segs{s});
                    fprintf(g, '\t\t\t\t\t<parameter idref="mutationRate.s:%s_3"/>\n', use_segs{s});
                end
             elseif contains(line, 'insert_param_log')
                 if ~isempty(est_tip_time)
                    for s = 1 : length(est_tip_time)
                        fprintf(g, '\t\t\t\t<log idref="tipprior.%s"/>\n', est_tip_time{s});
                    end
                 end

                for s = 1 : length(use_segs)                   
                    fprintf(g, '\t\t\t\t<log idref="kappa.s:%s_1"/>\n', use_segs{s});
                    fprintf(g, '\t\t\t\t<log idref="kappa.s:%s_3"/>\n', use_segs{s});
                    fprintf(g, '\t\t\t\t<log idref="mutationRate.s:%s_1"/>\n', use_segs{s});
                    fprintf(g, '\t\t\t\t<log idref="mutationRate.s:%s_3"/>\n', use_segs{s});
                    fprintf(g, '\t\t\t\t<log idref="gammaShape.s:%s_1"/>\n', use_segs{s});
                    fprintf(g, '\t\t\t\t<log idref="gammaShape.s:%s_3"/>\n', use_segs{s});
                    fprintf(g, '\t\t\t\t<log idref="freqParameter.s:%s_1"/>\n', use_segs{s});
                    fprintf(g, '\t\t\t\t<log idref="freqParameter.s:%s_3"/>\n', use_segs{s});
                end
            elseif contains(line, 'insert_tree_likelihood')
                for i = 1 : length(use_segs)
                    fprintf(g, '\t\t\t\t\t\t\t<distribution id="treeLikelihood.%s_1" spec="ThreadedTreeLikelihood" tree="@%s.tree" useAmbiguities="true">\n',use_segs{i},use_segs{i});
                    fprintf(g, '\t\t\t\t\t\t\t\t<data id="%s_1" spec="FilteredAlignment" filter="1::3,2::3">\n',use_segs{i});
                    fprintf(g, '\t\t\t\t\t\t\t\t\t<data idref="%s"/>\n',use_segs{i});
                    fprintf(g, '\t\t\t\t\t\t\t\t</data>\n');
                    fprintf(g, '\t\t\t\t\t\t\t\t<siteModel id="SiteModel.s:%s_1" spec="SiteModel"  gammaCategoryCount="4" shape="@gammaShape.s:%s_1" mutationRate="@mutationRate.s:%s_1">\n',use_segs{i},use_segs{i},use_segs{i});
                    fprintf(g, '\t\t\t\t\t\t\t\t\t<substModel id="hky.s:%s_1" spec="HKY" kappa="@kappa.s:%s_1">\n',use_segs{i},use_segs{i});
                    fprintf(g, '\t\t\t\t\t\t\t\t\t\t<frequencies id="estimatedFreqs.s:%s_1" spec="Frequencies" frequencies="@freqParameter.s:%s_1"/>\n',use_segs{i},use_segs{i});
                    fprintf(g, '\t\t\t\t\t\t\t\t\t</substModel>\n');
                    fprintf(g, '\t\t\t\t\t\t\t\t</siteModel>\n');
                    fprintf(g, '\t\t\t\t\t\t\t\t<branchRateModel id="StrictClock.%s_1" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c"/>\n',use_segs{i});
                    fprintf(g, '\t\t\t\t\t\t\t</distribution>\n');
                    fprintf(g, '\t\t\t\t\t\t\t<distribution id="treeLikelihood.%s_3" spec="ThreadedTreeLikelihood" tree="@%s.tree" useAmbiguities="true">\n',use_segs{i},use_segs{i});
                    fprintf(g, '\t\t\t\t\t\t\t\t<data id="%s_3" spec="FilteredAlignment" filter="3::3">\n',use_segs{i});
                    fprintf(g, '\t\t\t\t\t\t\t\t\t<data idref="%s"/>\n',use_segs{i});
                    fprintf(g, '\t\t\t\t\t\t\t\t</data>\n');
                    fprintf(g, '\t\t\t\t\t\t\t\t<siteModel id="SiteModel.s:%s_3" spec="SiteModel"  gammaCategoryCount="4" shape="@gammaShape.s:%s_3" mutationRate="@mutationRate.s:%s_3">\n',use_segs{i},use_segs{i},use_segs{i});
                    fprintf(g, '\t\t\t\t\t\t\t\t\t<substModel id="hky.s:%s_3" spec="HKY" kappa="@kappa.s:%s_3">\n',use_segs{i},use_segs{i});
                    fprintf(g, '\t\t\t\t\t\t\t\t\t\t<frequencies id="estimatedFreqs.s:%s_3" spec="Frequencies" frequencies="@freqParameter.s:%s_3"/>\n',use_segs{i},use_segs{i});
                    fprintf(g, '\t\t\t\t\t\t\t\t\t</substModel>\n');
                    fprintf(g, '\t\t\t\t\t\t\t\t</siteModel>\n');
                    fprintf(g, '\t\t\t\t\t\t\t\t<branchRateModel id="StrictClock.%s_3" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c"/>\n',use_segs{i});
                    fprintf(g, '\t\t\t\t\t\t\t</distribution>\n');
                end
            elseif contains(line, 'insert_weights')
                for i = 1 : length(use_segs)
                    fprintf(g, '%d %d ', round(seq_length(i)*2/3),round(seq_length(i)/3));
                end
            elseif contains(line, 'insert_seg_tree_state')
                for i = 1 : length(use_segs)
                     fprintf(g, '\t\t\t\t\t\t<stateNode id="%s.tree" spec="Tree" taxonset="@taxonSet" trait="@traitSet"/>\n', use_segs{i});                 
                end
%                 fprintf(g, '\t\t\t\t\t\t<stateNode id="dummy.tree" spec="Tree" taxonset="@taxonSet" trait="@traitSet"/>\n');                 
            elseif contains(line, 'insert_seg_tree')
                for i = 1 : length(use_segs)
                     fprintf(g, '\t\t\t\t\t<segmentTree idref="%s.tree"/>\n', use_segs{i});                 
                end
%                 fprintf(g, '\t\t\t\t\t<segmentTree idref="dummy.tree"/>\n'); 
            elseif contains(line, 'insert_seg_logger')
                for i = 1 : length(use_segs)
                    fprintf(g, '\t\t\t<logger spec="Logger" logEvery="50000" mode="tree" fileName="$(filebase).%s.trees">\n', use_segs{i});           
                    fprintf(g, '\t\t\t\t<log idref="%s.tree"/>\n', use_segs{i});           
                    fprintf(g, '\t\t\t</logger>\n');
                end

            elseif contains(line, 'insert_stats_log')
                for i = 1 : length(use_segs)
                    fprintf(g, '\t\t\t\t<log spec="TreeStatLogger" tree="@%s.tree"/>\n', use_segs{i});      
                end
            elseif contains(line, 'insert_seg_tree')
                for i = 1 : length(use_segs)
                     fprintf(g, '\t\t\t\t\t<segmentTree idref="%s.tree"/>\n', use_segs{i});                 
                end

            else
                fprintf(g, line);
            end
        end
        fclose(f);fclose(g);
    end 
end
end
