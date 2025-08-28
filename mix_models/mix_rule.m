function f = mix_rule(epsh,epsr, epsf, nu)
%MIX_RULE calculate volume fraction of inclusions based on mixing rule
%   using bulk rock glacier dielectric constant 
%INPUTS
%epsh = dielectric constant of half-space (background)
%epsr = dielectric constant of rock (inclusions)
%epsf = measured effective dielectric constant
%nu = mixing model identifier
%   0 = Maxwell-Garnet 
%   2 = Polder van-Santen
%   3 = Coherent potential approx. 


%OUTPUT
% f = volume fraction of inclusions (rock in rock glaciers)

%%%Reference: Shivola, 2008: EM Mixing Formulas & Applications, section
%%%9.3.1

A = (epsf-epsh)./(epsf+2*epsh+nu.*(epsf-epsh));
B = (epsr-epsh)./(epsr+2*epsh+nu.*(epsf-epsh));
f = A./B;

end

