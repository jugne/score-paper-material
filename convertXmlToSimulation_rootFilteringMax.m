function [] = convertXmlToSimulation_rootFilteringMax(folder, type_names)
% converts an xml to run a reassorty analyses to perform simulations using
% the same tip dates and the mean estimated reassortment rate and effective
% population size
cd(folder)
for run =1:10
%     cd(folder)
    log_files = dir(['current/h5n1_reject_30DaysDiff_08_16/run_' num2str(run) '/output/root_filter_max/h5n1_combined.rootFilter.log']);
    for i = 1 : length(log_files)
        t = importdata(['current/h5n1_reject_30DaysDiff_08_16/run_' num2str(run) '/output/root_filter_max/' log_files(i).name]);
        % look for the reassormtment and pop size headers
        rea_ind = cell(0,0); pop_ind = cell(0,0); mig_ind = cell(0,0);
        for j = 1 : length(t.textdata)
            for k = 1 : length(type_names)
                if strcmp(t.textdata{j}, ['reassortmentRate.' type_names{k}])
                    rea_ind{k} = j;
                elseif strcmp(t.textdata{j}, ['popSize.t.' type_names{k}])
                    pop_ind{k} = j;
                end
                for kk=1 : length(type_names)
                    if (kk~=k && strcmp(t.textdata{j}, ['b_migrationRate.' type_names{k} '_to_' type_names{kk}]))
                        mig_ind{k} = j;
                    end
                end
            end
        end
            
%         if rea_ind~=-cell(0,0)
            mean_rea = cell(0,0);
            mean_coal = cell(0,0);
            mean_mig = cell(0,0);
                
            for k = 1 : length(type_names)
                mean_rea{k} = mean(t.data(:,rea_ind{k}));
                mean_coal{k} = 1/mean(t.data(:,pop_ind{k}));
                mean_mig{k} = mean(t.data(:,mig_ind{k}));
            end

            % get the taxa and dates from the xmls corresponding to the log files
            f = fopen(['current/h5n1_reject_30DaysDiff_08_16/run_' num2str(run) '/xmls_coupled/' strrep(log_files(i).name, '_combined.rootFilter.log', '_rep0.xml')]);
            taxa_and_dates = cell(0,0);
            types = cell(0,0);
            while ~feof(f)
                line = fgets(f);
%                 if contains(line, '<taxonset id="taxonSet" spec="TaxonSet">')
%                     line = fgets(f);
%                     while ~contains(line, '</taxonset>')
%                         taxa{end+1} = line;
%                         line = fgets(f);
%                     end           
                if contains(line, '<trait spec="TraitSet" traitname="date-forward" id="traitSet">')
                    line = fgets(f);
                    while ~contains(line, '</trait>')
                        taxa_and_dates{end+1,1} = line;
                        line = fgets(f);
                    end
                elseif contains(line, '<typeTrait spec="TraitSet" traitname="type" id="typeTrait">')
                    line = fgets(f);
                    while ~contains(line, '<taxa idref="taxonSet"/>')
                        types{end+1,1} = line;
                        line = fgets(f);
                    end
                end
            end
            fclose(f);
            % build a simulation file with that info
            f = fopen('../simulation_template2.xml');
            system('mkdir simulation/posterior_rates/root_filter_max');
            system(sprintf('mkdir simulation/posterior_rates/root_filter_max/run_%s' , num2str(run)));
            g = fopen(['simulation/posterior_rates/root_filter_max/run_' num2str(run) '/' strrep(log_files(i).name, '.log', '.sim.xml')], 'w');
            while ~feof(f)
                line = fgets(f);
                if contains(line, 'insert_coalescent_rate')
                    fprintf(g, strrep(line, 'insert_coalescent_rate', num2str([mean_coal{:}])));
                elseif contains(line, 'insert_reassortment_rate')
                    fprintf(g, strrep(line, 'insert_reassortment_rate', num2str([mean_rea{:}])));
                elseif contains(line, 'insert_migration_rate')
                    fprintf(g, strrep(line, 'insert_migration_rate', num2str([mean_mig{:}])));
                elseif contains(line, 'insert_types') 
                    for j = 1 : length(types)
                        fprintf(g, '%s', types{j});
                    end
                elseif contains(line, 'insert_taxa_and_dates')
                    for j = 1 : length(taxa_and_dates)
                        fprintf(g, '%s', taxa_and_dates{j});
                    end
                else
                    fprintf(g, line);
                end        
            end
            fclose('all');
        end
    end
end