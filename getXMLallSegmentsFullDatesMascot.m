function [] = getXMLallSegmentsFullDatesMascot(virus, directory, fasta_file, template, runs, nrsequences, from, to, temperature, seeds)
    for run = 1 : runs
    rng(seeds(run));
    mkdir (sprintf('%s',['FullDatesInference/' directory]),sprintf('%s',['run_' num2str(run) '/xmls']));
    
    data_dir = ['/Users/jugne/Documents/SCORE-paper/StructCoalRe-Dataset_/h5n1/' fasta_file '/run_' num2str(run) '/'];
    segs_files = dir([data_dir '*.fasta']);
    
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
    
    use_seqs = zeros(nrsequences,1);
    for j = 1:length(use_seqs)
        use_seqs(j) = j;
    end
    
    use_segs = segments;
    for r = 0 : 2
        f = fopen(template);
        g = fopen(['FullDatesInference/' directory '/run_' num2str(run) '/xmls/' virus '_rep' num2str(r) '.xml'], 'w');
        
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
            elseif contains(line, 'insert_parameters')
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
           elseif contains(line, 'tree_init')
                for s = 1 : length(use_segs)
                    fprintf(g, '\t\t\t\t<init id="RandomTree.t:%s" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@%s.tree" taxa="@%s">\n',use_segs{s},use_segs{s},use_segs{s});
                    fprintf(g, '\t\t\t\t\t<populationModel id="ConstantPopulation0.t:%s" spec="ConstantPopulation">\n',use_segs{s});
                    fprintf(g, '\t\t\t\t\t\t<parameter id="randomPopSize.t:%s" name="popSize">1.0</parameter>\n',use_segs{s});
                    fprintf(g, '\t\t\t\t\t</populationModel>\n');
                    fprintf(g, '\t\t\t\t</init>\n');
                end
           elseif contains(line, 'insert_mascot_prior')
                for s = 1 : length(use_segs)              
                    fprintf(g, '\t\t\t\t<distribution id="structuredCoalescent.t:%s" spec="beast.mascot.distribution.Mascot" tree="@%s.tree">\n',use_segs{s},use_segs{s});
                    fprintf(g, '\t\t\t\t\t<dynamics id="constant.t:%s" spec="beast.mascot.dynamics.Constant" Ne="@popSize.t" backwardsMigration="@migrationRate" dimension="0" typeTrait="@typeTrait"/>\n',use_segs{s});
	                fprintf(g, '\t\t\t\t\t<structuredTreeIntervals id="TreeIntervals.t:%s" spec="beast.mascot.distribution.StructuredTreeIntervals" tree="@%s.tree"/>\n',use_segs{s},use_segs{s});
                    fprintf(g, '\t\t\t\t</distribution>\n');
                end
          elseif contains(line, 'insert_priors')
                if ~isempty(est_tip_time)
                    for s = 1 : length(est_tip_time)
                            fprintf(g, '\t\t\t\t\t\t\t\t<distribution id="tipprior.%s" spec="beast.math.distributions.MRCAPrior" tipsonly="true" tree="@ha.tree">\n', est_tip_time{s});
                            fprintf(g, '\t\t\t\t\t\t\t\t\t<taxonset id="tip.%s" spec="TaxonSet">\n', est_tip_time{s});
                            fprintf(g, '\t\t\t\t\t\t\t\t\t\t<taxon idref="%s" spec="Taxon"/>\n', est_tip_time{s});
                            fprintf(g, '\t\t\t\t\t\t\t\t\t</taxonset>\n');
                            tmp = strsplit(est_tip_time{s}, '|');
                            tmp_ = strsplit(tmp{4}, '-');
                            if ~isempty(strfind(tmp{4}, 'XX-XX'))
                                fprintf(g, '\t\t\t\t\t\t\t\t\t\t<Uniform id="Unform.%s" name="distr" lower="%s" upper="%d"/>\n', est_tip_time{s}, tmp_{1}, str2double(tmp_{1})+1);
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
                            fprintf(g, '\t\t\t\t\t\t\t\t\t\t<Uniform id="Unform.%s" name="distr" lower="%d" upper="%d"/>\n', est_tip_time{s}, lower, upper);
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
                        fprintf(g, '\t\t\t\t<operator spec="MultiTipDatesRandomWalker" useGaussian="true" windowSize="0.1" weight="0.01">\n');
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

                    fprintf(g, '\t\t\t\t<operator id="CoalescentConstantTreeScaler.t:%s" spec="ScaleOperator" scaleFactor="0.5" tree="@%s.tree" weight="3.0"/>\n', use_segs{s}, use_segs{s});
                    fprintf(g, '\t\t\t\t<operator id="CoalescentConstantTreeRootScaler.t:%s" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@%s.tree" weight="3.0"/>\n', use_segs{s}, use_segs{s});
                    fprintf(g, '\t\t\t\t<operator id="CoalescentConstantUniformOperator.t:%s" spec="Uniform" tree="@%s.tree" weight="30.0"/>\n', use_segs{s}, use_segs{s});
                    fprintf(g, '\t\t\t\t<operator id="CoalescentConstantSubtreeSlide.t:%s" spec="SubtreeSlide" tree="@%s.tree" weight="15.0"/>\n', use_segs{s}, use_segs{s});
                    fprintf(g, '\t\t\t\t<operator id="CoalescentConstantNarrow.t:%s" spec="Exchange" tree="@%s.tree" weight="15.0"/>\n', use_segs{s}, use_segs{s});
                    fprintf(g, '\t\t\t\t<operator id="CoalescentConstantWide.t:%s" spec="Exchange" isNarrow="false" tree="@%s.tree" weight="3.0"/>\n', use_segs{s}, use_segs{s});
                    fprintf(g, '\t\t\t\t<operator id="CoalescentConstantWilsonBalding.t:%s" spec="WilsonBalding" tree="@%s.tree" weight="3.0"/>\n', use_segs{s}, use_segs{s});


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
            elseif contains(line, 'insert_seg_tree')
                for i = 1 : length(use_segs)
                     fprintf(g, '\t\t\t\t\t<up idref="%s.tree"/>\n', use_segs{i});                 
                end
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
            else
                fprintf(g, line);
            end
        end
        fclose(f);fclose(g);
    end 
end
end