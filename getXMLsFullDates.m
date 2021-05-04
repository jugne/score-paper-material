clear
currdir = cd;

virus = 'h5n1';
directory = 'h5n1_reject_30DaysDiff_08_16';
template = 'template_rejection.xml';
fasta_file = 'fasta_fullDates_30DaysDiff_08_16_used';
n_runs = 10;
nrsequences = 200;
% from = 2005;
% to = 2015;
from = 2008;
to = 2016;
temperature = 0.01;
seeds = [975, 2785, 6324, 9134, 1270, 8147, 9058, 5469, 9575, 1576];


% getXMLallSegmentsFullDates(virus, directory, fasta_file, template, n_runs, nrsequences, from, to, temperature, seeds)

directory = 'h5n1_mascot_30DaysDiff_08_16';
template = 'template_mascot.xml';
fasta_file = 'fasta_fullDates_30DaysDiff_08_16_used';

getXMLallSegmentsFullDatesMascot(virus, directory, fasta_file, template, n_runs, nrsequences, from, to, temperature, seeds)
