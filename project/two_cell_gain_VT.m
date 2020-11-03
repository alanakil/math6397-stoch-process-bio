%% Analysis of correlation modulation in Balanced networks. 


clear

seed = 9;

% Number of neurons in ffwd layer
N=2;
Nx=4000;

% Ffwd connection probs
Px=[.1; .1];

% Correlation between the spike trains in the ffwd layer
c=0.;
% Timescale of correlation
taujitter=5;

Jxm1=[50]/sqrt(Nx);
Jxm2=[-50]/sqrt(Nx);

% Time (in ms) for sim
T=2000000;

% Time discretization
dt=.1; %ms

% FFwd spike train rate (in kHz)
rx1=10/1000;
rx2=10/1000;

% Number of time bins
Nt=round(T/dt);
time=dt:dt:T;

% Synaptic timescales
taux=10;
taue=8;
taui=4;

% Maximum number of spikes for all neurons
% in simulation. Make it 50Hz across all neurons
% If there are more spikes, the simulation will terminate
maxns=ceil(10*T); % was 0.05.

% Number of time bins to average over when recording
nBinsRecord=10;
dtRecord=nBinsRecord*dt;
timeRecord=dtRecord:dtRecord:T;
Ntrec=numel(timeRecord);

% Integer division function
IntDivide=@(n,k)(floor((n-1)/k)+1);

num_points_gain = 8;
threshold_vector = linspace(-65,-55,num_points_gain);
rx_vector = linspace(10/1000,20/1000,num_points_gain);

num_points_input = 1;
jestim_vector = linspace(-0.2,0.1,num_points_input);
mean_rate_1 = zeros(num_points_input,num_points_gain);
mean_rate_2 = zeros(num_points_input,num_points_gain);
Corr_1_2 = zeros(num_points_input,num_points_gain);

Cm=1;
gL=1/15;
EL=-72;
Vth=-50;
Vre=-75;
DeltaT=1;

winsize=250; 
T1=T/2; % Burn-in period
T2=T;   % Compute covariances until end of simulation

dtRate=1000; % ms


parfor k = 1:num_points_gain
    % Neuron parameters
 
    VT=threshold_vector(k);
    
    rx1 = rx_vector(1);
    rx2 = rx_vector(1);
    
    rng(seed);
    %%% Make (correlated) Poisson spike times for ffwd layer
    %%% See Appendix 1 of Akil et al 2020 for more details
    tic
    if(c<1e-5) % If uncorrelated
        nspikeX1=poissrnd(Nx*rx1*T);
        st=rand(nspikeX1,1)*T;
        sx1=zeros(2,numel(st));
        sx1(1,:)=sort(st);
        sx1(2,:)=randi(Nx,1,numel(st)); % neuron indices
    else % If correlated
        rm=rx1/c; % Firing rate of mother process
        nstm=poissrnd(rm*T); % Number of mother spikes
        stm=rand(nstm,1)*T; % spike times of mother process    
        maxnsx=round(T*rx1*Nx*1.2); % Max num spikes
        sx1=zeros(2,maxnsx);
        ns=0;
        for j=1:Nx  % For each ffwd spike train
            ns0=binornd(nstm,c); % Number of spikes for this spike train
            st=randsample(stm,ns0); % Sample spike times randomly
            st=st+taujitter*randn(size(st)); % jitter spike times
            st=st(st>0 & st<T); % Get rid of out-of-bounds times
            ns0=numel(st); % Re-compute spike count
            sx1(1,ns+1:ns+ns0)=st; % Set the spike times and indices        
            sx1(2,ns+1:ns+ns0)=j;
            ns=ns+ns0;
        end

        % Get rid of padded zeros
        sx1 = sx1(:,sx1(1,:)>0);

        % Sort by spike time
        [~,I] = sort(sx1(1,:));
        sx1 = sx1(:,I);
        nspikeX1=size(sx1,2);
    end
    tGenx=toc;
    disp(sprintf('\nTime to generate ffwd spikes: %.2f sec',tGenx))

    %%% Make (correlated) Poisson spike times for ffwd layer
    %%% See Appendix 1 of Akil et al 2020 for more details
    tic
    if(c<1e-5) % If uncorrelated
        nspikeX2=poissrnd(Nx*rx2*T);
        st=rand(nspikeX2,1)*T;
        sx2=zeros(2,numel(st));
        sx2(1,:)=sort(st);
        sx2(2,:)=randi(Nx,1,numel(st)); % neuron indices
    else % If correlated
        rm=rx2/c; % Firing rate of mother process
        nstm=poissrnd(rm*T); % Number of mother spikes
        stm=rand(nstm,1)*T; % spike times of mother process    
        maxnsx=round(T*rx2*Nx*1.2); % Max num spikes
        sx2=zeros(2,maxnsx);
        ns=0;
        for j=1:Nx  % For each ffwd spike train
            ns0=binornd(nstm,c); % Number of spikes for this spike train
            st=randsample(stm,ns0); % Sample spike times randomly
            st=st+taujitter*randn(size(st)); % jitter spike times
            st=st(st>0 & st<T); % Get rid of out-of-bounds times
            ns0=numel(st); % Re-compute spike count
            sx2(1,ns+1:ns+ns0)=st; % Set the spike times and indices        
            sx2(2,ns+1:ns+ns0)=j;
            ns=ns+ns0;
        end

        % Get rid of padded zeros
        sx2 = sx2(:,sx2(1,:)>0);

        % Sort by spike time
        [~,I] = sort(sx2(1,:));
        sx2 = sx2(:,I);
        nspikeX2=size(sx2,2);
    end
    tGenx=toc;
    disp(sprintf('\nTime to generate ffwd spikes: %.2f sec',tGenx))

    
    % Preallocate memory
    Ix1=zeros(N,1);
    Ix2=zeros(N,1);

    % Set initial voltage
    % Random initial voltages
    V0=rand(N,1)*(VT-Vre)+Vre;
    V=V0;

    iFspike1=1;
    iFspike2=1;
    s=zeros(2,maxns);
    nspike=0;
    TooManySpikes=0;

    J_NeNx = [ones(1,Nx/4),zeros(1,Nx/4),zeros(1,Nx/4),ones(1,Nx/4);...
        ones(1,Nx/4),ones(1,Nx/4),zeros(1,Nx/4),zeros(1,Nx/4)];    %binornd(1,Px(1),N,Nx2);
    % Generate full connectivity matrices
    tic
    Jx1=[J_NeNx.*Jxm1(1)];

    % Generate full connectivity matrices
    Jx2=[J_NeNx.*Jxm2(1)];

    tGen=toc;
    disp(sprintf('\nTime to generate connections: %.2f sec',tGen))

    % Extra stimulus: Istim is a time-dependent stimulus
    % it is delivered to all neurons with weights given by JIstim.
    % Specifically, the stimulus to neuron j at time index i is:
    % Istim(i)*JIstim(j)
    Istim=zeros(size(time)); 
    Istim(time>0)=0;
    jestim=1; 
    Jstim=sqrt(Nx)*[jestim]; 

    tic
    for i=1:numel(time)
        % Propogate E ffwd spikes
        while(sx1(1,iFspike1)<=time(i) && iFspike1<nspikeX1)
            jpre=sx1(2,iFspike1);
            Ix1=Ix1+Jx1(:,jpre)/taux;
            iFspike1=iFspike1+1;
        end

        % Propogate I ffwd spikes
        while(sx2(1,iFspike2)<=time(i) && iFspike2<nspikeX2)
            jpre=sx2(2,iFspike2);
            Ix2=Ix2+Jx2(:,jpre)/taux;
            iFspike2=iFspike2+1;
        end

        % Euler update to V
        V=V+(dt/Cm)*(Istim(i)*Jstim+Ix1+Ix2+gL*(EL-V)+gL*DeltaT*exp((V-VT)/DeltaT));

        % Find which neurons spiked
        Ispike=find(V>=Vth);    

        % If there are spikes
        if(~isempty(Ispike))

            % Store spike times and neuron indices
            if(nspike+numel(Ispike)<=maxns)
                s(1,nspike+1:nspike+numel(Ispike))=time(i);
                s(2,nspike+1:nspike+numel(Ispike))=Ispike;
            else
                TooManySpikes=1;
                break;
            end

            % Update cumulative number of spikes
            nspike=nspike+numel(Ispike);
        end            

        % Euler update to synaptic currents
        Ix1=Ix1-dt*Ix1/taux;
        Ix2=Ix2-dt*Ix2/taux;

        % This makes plots of V(t) look better.
        % All action potentials reach Vth exactly. 
        V(Ispike)=Vth;

        % Store recorded variables
        ii=IntDivide(i,nBinsRecord); 
    %     VRec(:,ii)=VRec(:,ii)+V(Irecord);

        % Reset membrane potential
        V(Ispike)=Vre;

        if mod(i*dt,T/5) == 0 % print every x iterations.
            fprintf('At time %d...\n',i*dt/T);
        end
    end
    % Get rid of padding in s
    s=s(:,1:nspike); 
    tSim=toc;
    disp(sprintf('\nTime for simulation: %.2f min',tSim/60))
    Tburn=T/2;
    r1Sim=hist(s(2,s(1,:)>Tburn & s(2,:)==1),1)/(T-Tburn);
    r2Sim=hist(s(2,s(1,:)>Tburn & s(2,:)==2),2)/(T-Tburn);

    % Time-dependent mean rates
    eRateT=hist(s(1,s(2,:)==1),dtRate:dtRate:T)/(dtRate);
    iRateT=hist(s(1,s(2,:)==2),dtRate:dtRate:T)/(dtRate);

    % Mean rate over E and I pops
    mean_rate_1(1,k)=mean(eRateT);
    mean_rate_2(1,k)=mean(iRateT);

    %% Compute spike count covariances and correlations
    % Compute spike count covariances over windows of size
    % winsize starting at time T1 and ending at time T2.
    tic
    C=SpikeCountCov(s,N,T1,T2,winsize);
    toc

    % Compute spike count correlations
    % This takes a while, so make it optional
    ComputeSpikeCountCorrs=1;
    if(ComputeSpikeCountCorrs)

        % Get correlation matrix from cov matrix
        tic
        R=corrcov(C);
        toc

        Corr_1_2(1,k) = R(1,2);

    end
end

%% Plot Corrs vs gain

figure; hold on
plot(threshold_vector, Corr_1_2 , 'linewidth',3)
xlabel('V_T')
ylabel('Corr(n1,n2)')


figure; hold on
plot(threshold_vector,1000*mean_rate_1,'linewidth',3,'color','r')
plot(threshold_vector,1000*mean_rate_2,'linewidth',3,'color','b')
xlabel('V_T')
ylabel('Rate (Hz)')


save('./gain_mod_twocell_VT.mat', 'c','Corr_1_2','DeltaT','dtRate',...
    'dt','EL','gL','jestim_vector','Jxm1','Jxm2','mean_rate_1',...
    'mean_rate_2','num_points_input','num_points_gain',...
    'Nx','Px','rx_vector','T','T1','T2','threshold_vector',...
    'Vre','Vth','winsize');

