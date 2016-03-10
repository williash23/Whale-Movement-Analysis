% load output files from WinBUGS as saved for the CODA 
% S-Plus diagnostic package. 
% Each MCMC chain is in a separate file showing the
% iteration number and value

load out-1.txt;
load out-2.txt;
load out-3.txt;
load out-4.txt;

load indkey; 	% this file contains a description of which 
% lines of the outup file correspond to 
% which variable - this is the CODA .ind file.

load elkdata	% an ascii file with observed steps and turning
% angles

nreps = 5000; 	% number of replicates for the posterior predictive
% check
elk = elkdata;
n = length(elk);	% size of movement path
sim = [];	% empty array to hold the values from the MCMC samples

% read MCMC samples
for k = 1:4
   simi = [];
   if k == 1
      a = out-1;
   elseif k==2
      a = out-2;
   elseif k==3
      a = out-3;
   else
      a = out-4;
   end
         
   for i = 1:length(key)
      simi = [simi a(key(i,1):key(i,2),2)];
   end
     
   sim = [sim; simi];
end

s = size(sim); 	% size of the MCMC samples (all chains)
% samples are in rows and variables in columns

% create some variables to hold results
sqd = ones(nreps,1).*NaN;
L = sqd;
LW = L;
LWC = L;
AC = [];
X = [];
Y = [];

hh = waitbar(0,'Please wait...');
for j = 1:nreps
   
   waitbar(j/nreps,hh)
   
   i = ceil(rand*s(1)); % choose a MCMC chain at random
   
   camp = find(sim(i,1:n)==1); 	% find observations classified as 
% "encamped"
   expl = find(sim(i,1:n)==2); 	% find observations classified as 
% "exploratory"
   
   % set some values to zero
   sqdev = 0;
   sqdeve = 0;
   lWc = 0;
   lWe = 0;
   lWCc = 0;
   lWCe = 0;
   simdatal = zeros(n,1);
   simdatat = zeros(n,1);
   
   if ~isempty(camp)
      
      % likelihoods (wcauchylike and weiblike return negative log 
% likelihoods)
      lWCc = 2 .* wcauchylike([sim(i,n+7) sim(i,n+3)],elk(camp,2));
      lWc = 2 .* WEIBLIKE([sim(i,n+5) sim(i,n+1)],elk(camp,1));
      
      % simulate values for step and turs using parameters from the 
% MCMC chain
      lpred = weibrnd(sim(i,n+5),sim(i,n+1),length(camp),1);
      tpred = wcauchy(sim(i,n+7),sim(i,n+3),length(camp),1);
      % squared deviations
      sqdev = (elk(camp,1)-lpred).^2+(elk(camp,2)-tpred).^2;
      simdatal(camp') = lpred;
      simdatat(camp') = tpred;
            
   end
   
   % do the same for exploratory state
   if ~isempty(expl)
      
      lWCe = 2 .* wcauchylike([sim(i,n+8) sim(i,n+4)],elk(expl,2));
      lWe = 2 .* WEIBLIKE([sim(i,n+6) sim(i,n+2)],elk(expl,1));
      lprede = weibrnd(sim(i,n+6),sim(i,n+2),length(expl),1);
      tprede = wcauchy(sim(i,n+8),sim(i,n+4),length(expl),1);
      sqdeve = (elk(expl,1)-lprede).^2+(elk(expl,2)-tprede).^2;
      simdatal(expl') = lprede;
      simdatat(expl') = tprede;
      
   end
   
   % build simulated movement paths
   x = zeros(n,1);
   y = x;
   dir = rand*2*pi;
   x(2) = cos(dir).*simdatal(1);
   y(2) = sin(dir).*simdatal(1);
   
   for k = 2:n
      x(k+1) = x(k) + cos(simdatat(k-1) + dir) .* simdatal(k);
      y(k+1) = y(k) + sin(simdatat(k-1) + dir) .* simdatal(k);
      dir = dir + simdatat(k-1);
   end
   
   % calculate and save the autocorrelation function
   AC = [AC; acf(simdatal)];
   X = [X x];
   Y = [Y y];
   
   % total squared deviations and likelihoods
   sqd(j) = sum(sum(sqdev)) + sum(sum(sqdeve));
   LWC(j) = sum(sum(lWCc)) + sum(sum(lWCe));
   LW(j) = sum(sum(lWc)) + sum(sum(lWe));
   L(j) = LW(j)+LWC(j);
   
end
close(hh)

% calculate Deviance for tetha hat
indi = median(sim(:,1:n));

camp = find(indi == 1);
expl = find(indi == 2);

lWc = 2 .* WEIBLIKE([mean(sim(:,n+5)) mean(sim(:,n+1))],elk(camp,1));
lWe = 2 .* WEIBLIKE([mean(sim(:,n+6)) mean(sim(:,n+2))],elk(expl,1));
lWCc = 2 .* wcauchylike([meandirection(sim(:,n+7)) mean(sim(:,n+3))],elk(camp,2));
lWCe = 2 .* wcauchylike([meandirection(sim(:,n+8)) mean(sim(:,n+4))],elk(expl,2));


Dtetha = lWc + lWe + lWCc + lWCe;

% calculate expected Deviance
Dbar = mean(L);

DIC = Dtetha + 2 * (Dbar - Dtetha);

% display Deviance results
[Dbar Dtetha Dbar-Dtetha DIC]

% plot the acf 
figure

x = 0:1:n-1;
x = x';
aca = acf(elk(:,1));
aca(1) = NaN;
ac = sort(AC);
ha=plot(x,aca,'.-k'); set(ha,'MarkerSize',20,'LineWidth',2);

hold on
ac(:,1) = NaN;
hacl = plot(x,ac(5000-125,:),'k'); set(hacl,'LineWidth',1)
hacu = plot(x,ac(125,:),'k'); set(hacu,'LineWidth',1)

AC(:,1) = NaN;
h = plot(AC','.k'); set(h,'MarkerSize',3);
AXIS([0 60 -.2 0.8]);

percentilAC = [ac(125,:) ac(5000-125,:)];

% save results
save DICelk L LWC LW lWc lWe lWCc lWCe DIC Dbar Dtetha percentilAC   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB function to simulate pseudo random numbers with Wrapped Cauchy distruibution

function [t] = wcauchy(mu,p,M,N)

%  [t] = wcauchy(mu,p,M,N)
%  pseudo-random number generation of the wrapped cauchy distribution with mean m and 
%  mean resultant lenght p. 
%  wcauchy(mu,p) returns a single value
%  wcauchy(mu,p,M,N) returns a M by N array 
%  The circular dispersion is
%  (1-p^2)/(2p^2)
%  circular variance v = 1-p
%  from Fisher(1993) Statistical analysis of circular data

if nargin == 2
   u = rand;
   V = cos(2*pi*u);
   c = 2*p/(1+p^2);
   
   t = sign(rand - .5) * acos((V+c)/(1+c.*V))+ mu;
   t = mod(t,2*pi);
   
   elseif nargin == 4
   
   u = rand(M,N);
   V = cos(2.*pi.*u);
   c = 2 .* p ./ (1 + p.^2);
   
   t =  sign(rand(M,N) - 0.5) .* acos((V+c)./(1+c.*V))+ mu;
   t = mod(t,2*pi);
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB code for negative log likelihood of Wrapped Cauchy

function logL = wcauchylike(params,data)
%   logL = wcauchylike(params,data)
%   log likelihood for wrapped Cauchy distribution

if nargin < 2, 
    error('Requires at least two input arguments'); 
end

[n, m] = size(data);

if nargout == 2 & max(m,n) == 1
    error('To compute the 2nd output, the 2nd input must have at least two elements.');
end

if n == 1
   data = data';
   n = m;
end

rho = params(2);
mu = params(1);

rho = rho(ones(n,1),:);

mu = mu(ones(n,1),:);

x = (1/(2*pi)) .* (1 - rho.^2)./(1+rho.^2 - 2.*rho.*cos(data-mu)) + eps;

logL = -sum(log(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m = meandirection(X)

%	m = meandirection(X)
%	from Fisher(1993), page 31
%	ave 8/28/00

C = sum(cos(X));
S = sum(sin(X));
%R2 = C^2 + S^2;
if C < 0
   m = atan(S/C) + pi;
elseif S > 0 
   m = atan(S/C);
else
   m = atan(S/C) + 2*pi;
end

