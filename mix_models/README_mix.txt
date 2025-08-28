Rock glacier dielectric mixing model analysis code

Tyler Meng

Plots results of field measurements applied to various dielectric mixing models. 
%%%%%%%%%%%%%%%%%%%%
Functions and Dependencies:

mix_rule.m = calculates the fraction of inclusions for a measured bulk dielectric constant given end member dielectric constants for the background and inclusion material, along with a model identifier which may only take the values of 0, 2, or 3, depending on the desired version of the unified mixing model described by (Sihvola, 2008). The plot_mix.m script is dependent upon this function for plotting the ice fraction results from the cmp analyses.

%%%%%%%%%%%%%%%%%%%%%%%
Scripts:

plot_mix.m: calculates the ice fraction results of four mixing model variationsand plots these results with the measured rock glacier wave speeds. Three mixing models are variations of the unified mixing model (Sihvola, 2008) using the mix_rule.m function, and we also calculate the ice fraction using the CRIM model (Knight and others, 1990). Generates Figure S7.
