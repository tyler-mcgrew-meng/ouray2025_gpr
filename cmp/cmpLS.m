function [v,t0,sig,cov,d_hat1] = cmpLS(x,t)
%cmpLS Estimate apparent NMO velocity and zero offset 
%two way travel time from CMP using least squares
%x and t must be vectors of equal length 

%v = best fit wave speed 
%t0 = best fit zero-offset two way travel time
%sig = RMS error in units of time squared
%cov = covariance
%d_hat1 = best fit data (for plotting)

%x = antenna separation 
%t = two way travel time at separation x

d2 = t.^2; %square data to linearize equation

G = [ones(length(t),1), x.^2]; %design matrix
m = inv(G'*G)*G'*d2;

%return best fit velocity and zero-offset travel time
t0 = sqrt(m(1)); 
v = sqrt(1/m(2));

%predict best fit data for vector x
d_hat1 = G*m;

%calculate residuals and asses error statistics
res1 = d2 - d_hat1; 
%variance (s^4)
Var = res1'*res1/(length(t)-length(m));
%standard deviation (units is s^2)
sig = sqrt(Var);
%covariance matrix
cov = Var*inv(G'*G);

end

