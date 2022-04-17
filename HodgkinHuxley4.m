function HodgkinHuxley4
% clear
count=1;
for Ne=160% excitatory neurons
    Ni=40;% inhib neurons
%% INITIALIZATION
DT = 0.025;     % time step (ms)
TMAX = 4000;      % simulation time (ms)
t = 0:DT:TMAX;

VREST = -65;      % resting potential

GNA = 100.0;	% NA channel conductivity (mS/cm^2)
GK = 80.0;	    % K channel conductivity (mS/cm^2)
GLEAK = 0.1;	% leakage "channel" conductivity (mS/cm^2)
GAHP = 0.3;	%

ENA = 50.0;	% Na reversal potential (mV)
EK = -100.0;   	% K reversal potential (mV)
ELEAK = -67.0;	% leakage current reversal potential (mV)

C_M = 1.0;	    % membrane capacitance (uF/cm^2)
% Ne=80;      
% Ni=20;    
N=Ne+Ni;
v = zeros(N,1);
m = zeros(N,1);
n = zeros(N,1);
h = zeros(N,1);
w = zeros(N,1);

%synaptic conductance
gei=0.05;
gee=0.01;
gie=5;
gii=10;

se=0; %total excit
si=0; %total inhib
Te=2.4;
Ti=12;

ae=20;ai=1;% /ms
sge = zeros(Ne,1);% synap gate
sgi = zeros(Ni,1);
filt1(1:length(t),1)=0;
v(:) = VREST;                                           % initial resting potential v
m(:) = alpha_m(VREST)/(alpha_m(VREST)+beta_m(VREST));	% initial gating variable m
h(:) = alpha_h(VREST)/(alpha_h(VREST)+beta_h(VREST));	% initial gating variable h
n(:) = alpha_n(VREST)/(alpha_n(VREST)+beta_n(VREST));	% initial hating variable n
w(:) = w_i(VREST);
counter=1;
noise(1:N,40002)=0;
vv = [ones(Ne,1);zeros(Ni,1)];%Dirac delta function to compute variable w for excitatory neurons
firings=[];% fired neurons: for rastor plot
flist=[]; % temporary variable, to save fired neurons
I=0; % synap. current
%% main simulation loop STARTS
for i=1:Ne
    noise(i,:)=0.5*wiener(40002)';
end;
for i=Ne+1:N
    noise(i,:)=0.0*wiener(40002)';
end;
for i=0:DT:TMAX
    x1=.6;x2=2; %excit stimulation
    random_vector1 = (x2-x1)*rand(Ne,1) + x1;
    x1=1.;x2=1.1;    %inhibitory stimulation
    random_vector2 = (x2-x1)*rand(Ni,1) + x1;
    Iapp=[random_vector1;random_vector2];
%         noise(1:Ne,1)=rand(Ne,1);
%         noise(Ne+1:N,1)=0.5*rand(Ni,1);
    counter=counter+1;
    if (mod(counter,40000)==0)
        for i=1:Ne
            noise(i,:)=0.5*wiener(40002)';
        end;
        for i=Ne+1:N
            noise(i,:)=0*wiener(40002)';
        end;
        counter=1;
    end;
    %  noise=0;
    gNa = GNA * m.^3 .* h; % Na channel conductance
    gK  = GK  .* n.^4;          % K channel conductance
    gAHP = GAHP * w;
    % note: the leakage conductance is constant
    mdot = alpha_m(v).*(1-m) - beta_m(v).*m; % differential equations governing the gating variables
    hdot = alpha_h(v).*(1-h) - beta_h(v).*h;
    ndot = alpha_n(v).*(1-n) - beta_n(v).*n;
    wdot = (w_i(v)-w)./w_tau(v);
    
    vdot = (noise(:,counter) + Iapp - I - gNa.*(v-ENA) - gK.*(v-EK) - GLEAK.*(v-ELEAK)-vv.*gAHP.*(v-EK))/C_M; % differential equation governing the  potential
    
    m = m + mdot*DT;		% Euler integration
    h = h + hdot*DT;
    n = n + ndot*DT;
    v = v + vdot*DT;
    w = w + wdot*DT;
    
    % compute synaptic current
    sge = sge + 0.025*( ( (ae*(1+tanh(v(1:Ne)/4))).*(1-sge) )-(sge/Te));
    sgi = sgi + 0.025*( ( (ai*(1+tanh(v(Ne+1:N)/4))).*(1-sgi) )-(sgi/Ti));
    se = sum(sge)/(Ne);%total excitatory
    si = sum(sgi)/Ni;%total inhibitory
    I(1:Ne,1)=(gee*se*v(1:Ne))+(gie*si*(v(1:Ne)+80));%total synaptic current for post-synaptic neuron i
    I(Ne+1:N,1)=(gei*se*v(Ne+1:N))+(gii*si*(v(Ne+1:N)+80));

    %0,0 - 0,1 - (1,1) - 1,0
    % check for fired neurons
    fired=find(v>30);
    dif = setdiff(fired,flist);% the spike is calculated once, see the below comment
    firings=[firings; i+0*dif,dif];
     flist=[flist; fired];% fired neurons are added to this variable
     flist=flist.*(v(flist)>0);% once it's action potential exceeds 30 mV
     flist=flist(find(flist>0));% and removed from the list when becomes < 0
end
%% main simulation loop ENDS

  
% 
  Fs = 1000;                    % Sampling frequency

 T = 1/Fs;                     % Sample time
 L = 1000;                     % Length of signal
 t = (0:L-1)*T;                % Time vector
 NFFT = 2^nextpow2(L); % Next power of 2 from length of y
 f = Fs/2*linspace(0,1,NFFT/2+1);
 times=zeros(40000,1);
 firings1=firings(find(firings(:,1)>3000),:);
  firings1(:,1)=firings1(:,1)-3000;
firings1(:,1)=firings1(:,1)*100000;
firings1(:,1)=floor( firings1(:,1));
% firings1(find(firings1(:,2)>160))=0;
counter=1;
for t=2500:2500:100000000
    times(counter)=size(find(firings1(:,1)==t),1);
    counter=counter+1;
end;
times1=downsample(times,40);
Y = fft(times1,NFFT)/L;
YY=2*abs(Y(1:NFFT/2+1));

end;
figure(1)
plot(firings1(:,1),firings1(:,2),'.')
% 'test'
figure(2)
plot(f(12:30),YY(12:30))

% Y1 = fft(filt1(40002:end),NFFT)/L;
% YY1=2*abs(Y1(1:NFFT/2+1));
% figure(3)
% plot(f(18:55),YY1(18:55))
% toc
%% ========================================================================
%% EMPIRICAL SUBFUNCTIONS
% they contribute to the gaiting control of the different channels

%% alpha_h
function rate = alpha_h(v)
rate = 0.128*exp(-(v+50.0)/18.0);

%% alpha_m
function rate = alpha_m(v)
rate = (0.32*(v+54.0)) ./ (1-exp(-(v+54.0)/4.0));

%% alpha_n
function rate = alpha_n(v)
rate = (0.032*(v+52.0)) ./ (1-exp(-(v+52.0)/5.0));

%% beta_h
function rate = beta_h(v)
rate =  4.0 ./ (exp(-(v+27.0)/5.0) + 1.0);

%% beta_m
function rate =  beta_m(v)
rate =  (0.28*(v+27.0)) ./ (exp((v+27.0)/5.0) - 1.0);

%% beta_n
function rate = beta_n(v)
rate = 0.5*exp(-(57+v)/40.0);

function rate = w_i(v)
rate = (1./(1+exp(-(v+35.0)/10.0)));

function rate = w_tau(v)
rate = (400./(3.3*(exp((v+35.0)/20.0))+exp(-(v+35.0)/20.0)));

