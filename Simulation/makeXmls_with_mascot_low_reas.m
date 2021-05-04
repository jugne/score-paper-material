% makes sim_compare xmls for the coalescent with reassortment in structured
% populations
clear

rng(1);

% define the number of repetitions
nr_reps = 100;

% rebuild the xml dirs
system('rm -r score mascot sim_compare');
system('mkdir score mascot sim_compare');

% define the number of samples
nr_samples = 100;

% define the sampling interval
sampling_interval = 20;

% file that keeps track of the ne and reassortment rates
h = fopen('rates_compare.csv', 'w');fprintf(h, 'run,Ne_1,Ne_2,reassortment_1,reassortment_2,migration_1, migration_2 \n');

% define params of the lognormal distribution of the Ne
mean_ne = 2;
sigma_ne = 0.5;
mu_ne = log(mean_ne) - sigma_ne^2/2;

% define params of the lognormal distribution of the reassortment rates
% high
% mean_rea = 0.5;
% low
mean_rea = 0.1;
sigma_rea = 0.5;
mu_rea = log(mean_rea) - sigma_rea^2/2;

% define params of the lognormal distribution of the migration rates
mean_mig = 0.2;
sigma_mig = 0.5;
mu_mig = log(mean_mig) - sigma_mig^2/2;

% define which evolutionary rates to use
use_rates = {'high','high','high','high';
            'high','high','low','low';
            'low','low','low','low'};
evol_name = {'high', 'mixed', 'low'};


% make nr reps number of xmls
for i = 1 : nr_reps
    f_sim = fopen('sim_template.xml');
    
    % sample the Ne and reassortment rates
    Ne_1 = lognrnd(mu_ne,sigma_ne);
    Ne_2 = lognrnd(mu_ne,sigma_ne);
    reassortment_1 = lognrnd(mu_rea, sigma_rea);
    reassortment_2 = lognrnd(mu_rea, sigma_rea);
    migration_1 = lognrnd(mu_mig, sigma_mig);
    migration_2 = lognrnd(mu_mig, sigma_mig);
%     migration_1 = exprnd(mean_mig);
%     migration_2 = exprnd(mean_mig);
    
    fprintf(h, '%d,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n', i, Ne_1, Ne_2, reassortment_1, reassortment_2, migration_1, migration_2);
    
    % open the sim_compare xml
    g = fopen(sprintf('sim_compare/sim_%d.xml', i), 'w');
    
    while ~feof(f_sim)
        line = fgets(f_sim);
        if contains(line, 'insert_taxa')
            for j = 1 : nr_samples
                fprintf(g,'\t\t\t\t<taxon spec="Taxon" id="t%d"/>\n', j);
            end             
        elseif contains(line, 'insert_sampling_times')
            time = zeros(nr_samples,1);
            for j = 1 : nr_samples
                time(j) = rand*sampling_interval;
                if j==nr_samples
                    fprintf(g,'\t\t\t\tt%d=%f\n', j, time(j));
                else
                    fprintf(g,'\t\t\t\tt%d=%f,\n', j, time(j));
                end
            end 
        elseif contains(line, 'insert_types')
            type = zeros(nr_samples,1);
            for j=1:nr_samples
                type(j) = randi([0,1]);
                if j==nr_samples
                    fprintf(g,'\t\t\t\tt%d=%i\n', j, type(j));
                else
                    fprintf(g,'\t\t\t\tt%d=%i,\n', j, type(j));
                end
            end   
        elseif contains(line, 'insert_coalescence')
            fprintf(g, strrep(line, 'insert_coalescence', num2str(1/Ne_1)+" "+num2str(1/Ne_2)));
        elseif contains(line, 'insert_reassortment')
            fprintf(g, strrep(line, 'insert_reassortment', num2str(reassortment_1)+" "+num2str(reassortment_2)));
        elseif contains(line, 'insert_migration')
            fprintf(g, strrep(line, 'insert_migration', num2str(migration_1)+" "+num2str(migration_2)));
        else
            fprintf(g, '%s', line);
        end
    end
    fclose(f_sim); fclose(g);
    

    for uncertainty = 1 : length(evol_name)
        % make 3 replicates
%         for r = 1 : 3
            % build the inf_compare xml
            f_inf = fopen('inf_template.xml');
            tmp_ne_1 = exprnd(1);
            tmp_ne_2 = exprnd(1);

            % open the inf_compare xml
            g = fopen(sprintf('score/inf_%s_%d.xml', evol_name{uncertainty}, i), 'w');
            
            % keep track of the segment count for the nexus file name
            segmentcount = 1;segmentcount2=1;
            
            while ~feof(f_inf)
                line = fgets(f_inf);
                if contains(line, 'insert_taxa')
                    for j = 1 : nr_samples
                        fprintf(g,'\t\t\t\t<taxon spec="Taxon" id="t%d"/>\n', j);
                    end             
                elseif contains(line, 'insert_sampling_times')
                    for j = 1 : nr_samples
                        if j==nr_samples
                            fprintf(g,'\t\t\t\tt%d=%f\n', j, time(j));
                        else
                            fprintf(g,'\t\t\t\tt%d=%f,\n', j, time(j));
                        end
                    end
                elseif contains(line, 'insert_types')
                    for j=1:nr_samples
                        if j==nr_samples
                            fprintf(g,'\t\t\t\tt%d=%i\n', j, type(j));
                        else
                            fprintf(g,'\t\t\t\tt%d=%i,\n', j, type(j));
                        end
                    end
                elseif contains(line, 'initial_Ne')
                    fprintf(g, strrep(line, 'initial_Ne', num2str(tmp_ne_1)+" "+num2str(tmp_ne_2)));
                elseif contains(line, 'initial_coalescence')
                    fprintf(g, strrep(line, 'initial_coalescence', num2str(1/tmp_ne_1)+" "+num2str(1/tmp_ne_2)));
                elseif contains(line, 'initial_reassortment')
                    fprintf(g, strrep(line, 'initial_reassortment', reassortment_1+" "+reassortment_2));
%                     fprintf(g, strrep(line, 'initial_reassortment', num2str(exprnd(1))+" "+num2str(exprnd(1))));
                elseif contains(line, 'initial_migration')
                    fprintf(g, strrep(line, 'initial_migration', num2str(exprnd(1))+" "+num2str(exprnd(1))));
                elseif contains(line, 'insert_sim_file_name')
                    line = strrep(line, 'insert_sim_file_name', sprintf('../sim_compare/sim_%d', i) );
                    line = strrep(line, 'inset_evol_rates',  use_rates{uncertainty, segmentcount});
                    fprintf(g, line);
                    segmentcount = segmentcount + 1;
                elseif contains(line, 'insert_clock_rate')
                    line = strrep(line, 'insert_clock_rate',  use_rates{uncertainty, segmentcount2});
                    fprintf(g, line);
                    segmentcount2 = segmentcount2 + 1;

                else
                    fprintf(g, '%s', line);
                end
            end

            fclose(f_inf); fclose(g);
%         end
    end
    

    
    for uncertainty = 1 : length(evol_name)
        if contains(evol_name{uncertainty}, 'mixed')
            rate = {'low', 'high'};
        else
            rate = {'high'};
        end
        for r = 1 : length(rate)
        % make 3 replicates
%         for r = 1 : 3
            % build the inf_compare xml
            f_inf = fopen('inf_template_mascot.xml');
            tmp_ne_1 = exprnd(1);
            tmp_ne_2 = exprnd(1);

            % open the inf_compare xml
            if length(rate) > 1
                g = fopen(sprintf('mascot/inf_%s_%s_%d.xml', evol_name{uncertainty}, rate{r}, i), 'w');
            else
                g = fopen(sprintf('mascot/inf_%s_%d.xml', evol_name{uncertainty}, i), 'w');
            end
            
            % keep track of the segment count for the nexus file name
            segmentcount = 1;segmentcount2=1;
            
            while ~feof(f_inf)
                line = fgets(f_inf);
                if contains(line, 'insert_taxa')
                    for j = 1 : nr_samples
                        fprintf(g,'\t\t\t\t<taxon spec="Taxon" id="t%d"/>\n', j);
                    end             
                elseif contains(line, 'insert_sampling_times')
                    for j = 1 : nr_samples
                        if j==nr_samples
                            fprintf(g,'\t\t\t\tt%d=%f\n', j, time(j));
                        else
                            fprintf(g,'\t\t\t\tt%d=%f,\n', j, time(j));
                        end
                    end
                elseif contains(line, 'insert_types')
                    for j=1:nr_samples
                        if j==nr_samples
                            fprintf(g,'\t\t\t\tt%d=%i\n', j, type(j));
                        else
                            fprintf(g,'\t\t\t\tt%d=%i,\n', j, type(j));
                        end
                    end
                elseif contains(line, 'initial_Ne')
                    fprintf(g, strrep(line, 'initial_Ne', num2str(tmp_ne_1)+" "+num2str(tmp_ne_2)));
                elseif contains(line, 'initial_coalescence')
                    fprintf(g, strrep(line, 'initial_coalescence', num2str(1/tmp_ne_1)+" "+num2str(1/tmp_ne_2)));
                elseif contains(line, 'initial_migration')
                    fprintf(g, strrep(line, 'initial_migration', num2str(exprnd(1))+" "+num2str(exprnd(1))));
                elseif contains(line, 'insert_sim_file_name')
                    line = strrep(line, 'insert_sim_file_name', sprintf('../sim_compare/sim_%d', i) );
                    line = strrep(line, 'inset_evol_rates',  use_rates{uncertainty, r});
                    fprintf(g, line);
%                     segmentcount = segmentcount + 1;
                elseif contains(line, 'insert_clock_rate')
                    line = strrep(line, 'insert_clock_rate',  use_rates{uncertainty, r});
                    fprintf(g, line);
%                     segmentcount2 = segmentcount2 + 1;

                else
                    fprintf(g, '%s', line);
                end
            end

            fclose(f_inf); fclose(g);
%         end
        end
    end

end
fclose(h);
