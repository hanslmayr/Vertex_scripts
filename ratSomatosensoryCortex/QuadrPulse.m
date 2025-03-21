%% This is a tutorial script from the Vertex toolbox that has been modified by Simon Hanslmayr
% 
% I have applied a number of changes to setUpBipolarElectrodeStimulation,
% so that you can input the desired stimulation intensity and you would need 
% that modified file for it to run;
%
% I also one stimulation pulse (not 4) which starts at 500 ms and induces a
% constant current for 1000 ms;
%
% This script only computes the data and stores it in the specified folders
% (see below); To visualise the results use Comp_Stim_Intens.m

% To run full simulation will require a high performance computing node.
% This has been run on 12 cores with 200GB RAM.
% Smaller scale models can be run by reducing the neuron density in
% loadRatTissueandNeuronParams.
% 
%% You May have to add the path 
addpath(genpath('I:\scripts\Matlab\tbx\Vertex_git-master'))

% Setting the model Tissue, Neuron, and Connection parameters.
% These are based on the Neocortical Collaborative Portal and are stored in
% connectionslayers23to6.mat, rat_no_neurons.mat, and ratlayerthickness.mat

loadRatTissueandNeuronParams;
TissueParams.neuronDensity = 103730;
loadRatConnectionParams;
%% Load the stimulating electrode field and set on and off times
% Default times are for single pulse at 1500 ms.
mA=80;
[TP]=setUpBipolarElectrodeStimulation(mA);
%% For paired pulse stimulation set the on and off times. 
% inherit parameters from setUpBipolarElectrodeStimulation
TissueParams.StimulationOn = TP.StimulationOn; 
TissueParams.StimulationOff = TP.StimulationOff;
TissueParams.StimulationField=TP.StimulationField;
TissueParams.scale=TP.scale;
% e.g. for 100 ms interval:
% TissueParams.StimulationOn = [500 600 ]; 
% TissueParams.StimulationOff = [ 500.5 600.5];
% e.g. for 150 ms interval:
TissueParams.StimulationOn = [500]; 
TissueParams.StimulationOff = TissueParams.StimulationOn + 1000;
% e.g. for 200 ms interval:
% TissueParams.StimulationOn = [500 700 ]; 
% TissueParams.StimulationOff = [ 500.5 700.5];
% e.g. for 250 ms interval:
% TissueParams.StimulationOn = [500 750 ]; 
% TissueParams.StimulationOff = [ 500.5 750.5];

loadSimulationSettings

%% Set the random seed, this should go from 1001 to 1005 to reproduce results in the paper where they are based on repeated simulations. 
SimulationSettings.randomSeed = 1001;
out_fn=['I:\WT_resub\Vertex_simulations\F1000_Data\pairedpulse\cont1sec_' num2str(mA) 'microAmps_' num2str(SimulationSettings.randomSeed)];
%% Set the location to save results.
RecordingSettings.saveDir = out_fn;
%% For single or paired pulse
runRatSimulationQP
% For theta burst stimulation
% runRatSimulationTBS

