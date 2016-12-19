% Plot Posterior SAMPLES
% original file: mcmc_posterior_plot
hFig = figure;

vars = PosteriorSamples.Properties.VariableNames;
totalvars = length(vars);
for i=1:length(vars) 
    if totalvars<=8
        subplot(totalvars,1,i);
        y = PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i));
        plot(table2array(y));
        str = char(vars(i));
        if strfind(str,'immune_alpha')
           str = 'immune_\alpha'; 
        end
        if strfind(str,'immune_beta')
           str = 'immune_\beta';
        end
            xlabel(str);
    end

    if totalvars>8
        subplot(5,2,i);
        y = PosteriorSamples(:,PosteriorSamples.Properties.VariableNames(i));
        plot(table2array(y));
        xlabel(char(vars(i)));
    end
end

set(hFig, 'Position', [10 10 800 120*totalvars])