% Single_Cell_JIVE.m
%  Author: Meilei

%% Data Processing
load Grun_reads_cutoff.mat

datablock{1} = Expression;
datablock{2} = Spikein; 

%% JIVE Integration
% plot scree plot 
JIVEPreVisualQF(datablock);

% try rank 
rank{1} = [4 10 20];
rank{2} = [2 3];
JIVEPreVisualQF(datablock,rank);

% select rank as vecr = [2, 3]
paramstruct = struct('ioutput', [1, 1, 1, 1, 1, 1, 1, 1, 1]);
vecr = [20, 3];
outstruct = JIVEMainQF(datablock, vecr, paramstruct);
rjoint = 1;
rIndExpression = 19;
rIndSpikein = 2;

%% visualize matrices
joint_Expression = outstruct.joint{1};
indiv_Expression = outstruct.individual{1};

joint_Spikein = outstruct.joint{2};
indiv_Spikein = outstruct.individual{2};

noise_Expression = Expression - joint_Expression - indiv_Expression;
noise_Spikein = Spikein - joint_Spikein - indiv_Spikein;

JIVEdecompVisualQF(joint_Expression, indiv_Expression)

HeatmapVisualQF(Expression,'Expression Matrix')
HeatmapVisualQF(joint_Expression,'Expression Joint')
HeatmapVisualQF(indiv_Expression,'Expression Individual')
HeatmapVisualQF(noise_Expression,'Expression Noise')

HeatmapVisualQF(Spikein,'Spikein Matrix')
HeatmapVisualQF(joint_Spikein,'Spikein Joint')
HeatmapVisualQF(indiv_Spikein,'Spikein Individual')
HeatmapVisualQF(noise_Spikein,'Spikein Noise')
%% save the result
save grun_reads_cutoff_reasult1.mat joint_Expression joint_Spikein indiv_Expression indiv_Spikein rjoint rIndExpression rIndSpikein;
