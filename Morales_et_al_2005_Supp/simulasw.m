% Juan M. Morales, Daniel T. Haydon, Jacqui Frair, Kent E. Holsinger and John M. Fryxell
% 2003. Extracting More out of Relocation Data: Building Movement Models as Mixtures of  Random Walks. 
% Ecology VOL: pp-pp.
%
% simulation of two random walks with switching probabilities

% define parameters
nsteps = 202; % number of steps to run
nind = 10;    % number of individuals to run

% population level parameters
% parameters for step length
af = 0.03;           % mean scale parameter for fast movement
as = 2.7;            % mean scale parameter for slow movement
bf = 2;              % mean shape parameter for fast movement
bs = 1;              % mean shape parameter for slow movement
sigma_af = 0.003;    % SD for distribution of scale parameters in fast movement
sigma_as = 0.27;     % SD for distribution of scale parameters in slow movement
sigma_bf = 0.2;      % SD for distribution of shape parameters in fast movement
sigma_bs = 0.1;      % SD for distribution of shape parameters in slow movement

% parameters for turning angles
muf = 0;             % mean turn for fast movement
mus = pi;            % mean turn for slow movement
sigma_muf = 0.05;    % SD for distribution of mean turns in fast movement 
sigma_mus = 0.1;     % SD for distribution of mean turns in slow movement 
arhof = 45;          % parameter a for the Beta distribution of mean cosines of turning angles for fast movement
arhos = 15;          % parameter a for the Beta distribution of mean cosines of turning angles for slow movement
brhof = 5;           % parameter b for the Beta distribution of mean cosines of turning angles for fast movement
brhos = 15;          % parameter b for the Beta distribution of mean cosines of turning angles for slow movement

% parameters for switching rates
afs = 15;            % parameter a for the Beta distribution of switching rates from slow to fast movement
asf = 5;             % parameter a for the Beta distribution of switching rates from slow to fast movement
bfs = 7.5;           % parameter b for the Beta distribution of switching rates from fast to slow movement
bsf = 15;            % parameter b for the Beta distribution of switching rates from fast to slow movement

% empty arrays for results (not very efficient Matlab code)
res = [];
XX = [];
YY = [];
indis = [];

%a = p;
%b = q; % fast to slow

% generate parameter values for all the simulated individuals
% step length parameters
aaf = normrnd(af,sigma_af,nind,1);
aas = normrnd(as,sigma_as,nind,1);
bbf = normrnd(bf,sigma_bf,nind,1);
bbs = normrnd(bs,sigma_bs,nind,1);

% turning angle parameters
mmuf = normrnd(muf,sigma_muf,nind,1);
mmus = normrnd(mus,sigma_mus,nind,1);
rrhof = betarnd(arhof,brhof,nind,1);
rrhos = betarnd(brhos,brhos,nind,1);

% switching parameters
aa = betarnd(asf,bsf,nind,1);
bb = betarnd(afs,bfs,nind,1);

% save space for results 
steps = zeros(nsteps,nind);
turns = steps;
indis = steps;

% iterate over all individuals
for j = 1:nind
    
    a = aa(j); % slow to fast
    b = bb(j); % fast to slow
    
    eqS = b/(a+b); % proportion of individuals in slow state at eq
    eqF = a/(a+b);
    
    if rand <= eqF
        indis(1,j) = 1; % fast movement
    else
        indis(1,j) = 2;
    end
    %indi = 1; %rand(nsteps+1,1)<0.8;
    clue = zeros(nsteps,1);
    x = 0;
    y = 0;
    %dir = rand;
    X = [];
    Y = [];
    t = zeros(nsteps,1);
    dir = zeros(nsteps+1,1);
    dir(1) = rand * 2*pi;
    dsq = zeros(nsteps,1);
    
    for i = 1:nsteps-1
        
        if indis(i,j) ~= 1
            
            % move according to movement mode
            %largo = normrnd(0, fl); %exprnd(fl); % %
            largo = weibrnd(aaf(j),bbf(j));
            t(i) = wcauchy(mmuf(j),rrhof(j));
            dir(i+1) = dir(i) + t(i);
            ct = cos(dir(i+1));
            st = sin(dir(i+1));
            turns(i,j) = t(i);
            steps(i,j) = largo;
            
            % find out whether there is a change in mode
            if rand < b
                indis(i+1,j) = 1;
            else
                indis(i+1,j) = 2;
            end
            
            
        else
            %largo = normrnd(0, sl); %exprnd(sl); % 
            largo = weibrnd(aas(j),bbs(j));
            t(i) = wcauchy(mmus(j),rrhos(j));
            dir(i+1) = dir(i) + t(i);
            ct = cos(dir(i+1));
            st = sin(dir(i+1));
            turns(i,j) = t(i);
            steps(i,j) = largo;
            
            if rand < a
                indis(i+1,j) =  2;
            else
                indis(i+1,j) = 1;
            end
            
        end
        
        x = x + ct * largo;
        y = y + st * largo;
        X = [X;x];
        Y = [Y;y];
        
    end
    dsq = (X-X(1)).^2+(Y-Y(1)).^2;
    res = [res, dsq];
    XX = [XX; X'];
    YY = [YY; Y'];
end
