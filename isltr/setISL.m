restoredefaultpath;

% Set path
p = path;
path(p,[pwd '/lib/']);
p = path;
path(p,[pwd '/lib/model']);
p = path;
path(p,[pwd '/lib/model/llh']);
p = path;
path(p,[pwd '/lib/model/susc']);
p = path;
path(p,[pwd '/lib/model/rt']);
p = path;
path(p,[pwd '/lib/chart']);
p = path;
path(p,[pwd '/lib/optim']);
p = path;
path(p,[pwd '/lib/sys']);

p = path;
path(p,[pwd '/main/likelihood_estimation']);
p = path;
path(p,[pwd '/main/calculate_R0']);
p = path;
path(p,[pwd '/main/extract_antibody_titres']);
p = path;
path(p,[pwd '/main/mcmc']);
p = path;
path(p,[pwd '/main/mcmc/mode2']);
p = path;
path(p,[pwd '/main/mcmc/mode4']);
p = path;
path(p,[pwd '/main/plot']);
p = path;
path(p,[pwd '/main/mplot']);
p = path;
path(p,[pwd '/main/mplot/back']);
p = path;
path(p,[pwd '/main/produce_samples']);

% Declare global variables
global proj dat Antibody
proj = 'dat/20141206/hk_ph1n1/'; %exlude vaccinated individuals; paired sera
load([proj 'h1n1_titres.mat']);


