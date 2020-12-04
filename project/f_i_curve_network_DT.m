%% Analysis of correlation modulation in Balanced networks. 

clear

seed = 39;

% Number of neurons in each population
N = 4000;
Ne=0.8*N;
Ni=0.2*N;

% Number of neurons in ffwd layer
Nx=0.2*N;

% Recurrent net connection probabilities
P=[0.1 0.1; 0.1 0.1];

% Ffwd connection probs
Px=[.1; .1];

% Correlation between the spike trains in the ffwd layer
c=0.1;
% Timescale of correlation
taujitter=5;

% Mean connection strengths between each cell type pair
Jm=[25 -150; 112.5 -250]/sqrt(N);
Jxm=[180; 135]/sqrt(N);

% Time (in ms) for sim
T=200000;

% Time discretization
dt=.1; %ms

% Proportion of neurons in each population.
qe=Ne/N;
qi=Ni/N;
qf=Nx/N;

% Number of time bins
Nt=round(T/dt);
time=dt:dt:T;

% Build mean field matrices
Q=[qe qi; qe qi];
Qf=[qf; qf];
W=P.*(Jm*sqrt(N)).*Q;
Wx=Px.*(Jxm*sqrt(N)).*Qf;

% Synaptic timescales
taux=10;
taue=8;
taui=4;

% Neuron parameters
Cm=1;
gL=1/15;
EL=-72;
Vth=-50;
Vre=-75;
DeltaT=1;
VT=-55;

% Maximum number of spikes for all neurons
% in simulation. Make it 50Hz across all neurons
% If there are more spikes, the simulation will terminate
maxns=ceil(.05*N*T); % was 0.05.

% Number of time bins to average over when recording
nBinsRecord=10;
dtRecord=nBinsRecord*dt;
timeRecord=dtRecord:dtRecord:T;
Ntrec=numel(timeRecord);

% Integer division function
IntDivide=@(n,k)(floor((n-1)/k)+1);

num_points_input = 8;
jestim_vector = linspace(-0.05,0.4,num_points_input);
rx_vector = linspace(5/1000,50/1000,num_points_input);

num_points_gain = 3;
threshold_vector = linspace(-65,-55,num_points_gain);
delta_vector = linspace(0.1,10,num_points_gain);


mean_e_rate = zeros(num_points_input,num_points_gain);
mean_i_rate = zeros(num_points_input,num_points_gain);
mRee = zeros(num_points_input,num_points_gain);
mRei = zeros(num_points_input,num_points_gain);
mRii = zeros(num_points_input,num_points_gain);

% Get mean spike count covariances over each sub-pop
[II,JJ]=meshgrid(1:N,1:N);

winsize=250; 
T1=T/2; % Burn-in period
T2=T;   % Compute covariances until end of simulation

Istim=zeros(size(time)); 
Istim(time>0)=0;

Tburn=T/2;

VT=-55;


for thr_indx = 1:num_points_gain
    
    DeltaT=delta_vector(thr_indx);
    
    parfor k=1:num_points_input

        rng(seed);

        % FFwd spike train rate (in kHz)
        rx=rx_vector(k);

        %%% Make (correlated) Poisson spike times for ffwd layer
        %%% See Appendix 1 of Akil et al 2020 for more details
        tic
        if(c<1e-5) % If uncorrelated
            nspikeX=poissrnd(Nx*rx*T);
            st=rand(nspikeX,1)*T;
            sx=zeros(2,numel(st));
            sx(1,:)=sort(st);
            sx(2,:)=randi(Nx,1,numel(st)); % neuron indices
        else % If correlated
            rm=rx/c; % Firing rate of mother process
            nstm=poissrnd(rm*T); % Number of mother spikes
            stm=rand(nstm,1)*T; % spike times of mother process    
            maxnsx=round(T*rx*Nx*1.2); % Max num spikes
            sx=zeros(2,maxnsx);
            ns=0;
            for j=1:Nx  % For each ffwd spike train
                ns0=binornd(nstm,c); % Number of spikes for this spike train
                st=randsample(stm,ns0); % Sample spike times randomly
                st=st+taujitter*randn(size(st)); % jitter spike times
                st=st(st>0 & st<T); % Get rid of out-of-bounds times
                ns0=numel(st); % Re-compute spike count
                sx(1,ns+1:ns+ns0)=st; % Set the spike times and indices        
                sx(2,ns+1:ns+ns0)=j;
                ns=ns+ns0;
            end

            % Get rid of padded zeros
            sx = sx(:,sx(1,:)>0);

            % Sort by spike time
            [~,I] = sort(sx(1,:));
            sx = sx(:,I);
            nspikeX=size(sx,2);
        end
        tGenx=toc;
        disp(sprintf('\nTime to generate ffwd spikes: %.2f sec',tGenx))


        % Preallocate memory
        Ie=zeros(N,1);
        Ii=zeros(N,1);
        Ix=zeros(N,1);
        % VRec=zeros(numrecord,Ntrec);

        % Set initial voltage
        % Random initial voltages
        V0=rand(N,1)*(VT-Vre)+Vre;
        V=V0;

        iFspike=1;
        s=zeros(2,maxns);
        nspike=0;
        TooManySpikes=0;

        J_NeNx = binornd(1,P(1,1),Ne,Nx);
        J_NiNx = binornd(1,P(1,1),Ni,Nx);
        % Generate full connectivity matrices
        tic
        J=[Jm(1,1)*binornd(1,P(1,1),Ne,Ne) Jm(1,2)*binornd(1,P(1,1),Ne,Ni); ...
           Jm(2,1)*binornd(1,P(1,1),Ni,Ne) Jm(2,2)*binornd(1,P(1,1),Ni,Ni)];
        Jx=[J_NeNx.*Jxm(1); J_NiNx.*Jxm(2)];

        tGen=toc;
        disp(sprintf('\nTime to generate connections: %.2f sec',tGen))

        % Extra stimulus: Istim is a time-dependent stimulus
        % it is delivered to all neurons with weights given by JIstim.
        % Specifically, the stimulus to neuron j at time index i is:
        % Istim(i)*JIstim(j)
        jestim=jestim_vector(k); 
        jistim=jestim_vector(k);
        Jstim=sqrt(N)*[jestim*ones(Ne,1); jistim*ones(Ni,1)]; 

        tic
        for i=1:numel(time)
            % Propogate ffwd spikes
            while(sx(1,iFspike)<=time(i) && iFspike<nspikeX)
                jpre=sx(2,iFspike);
                Ix=Ix+Jx(:,jpre)/taux;
                iFspike=iFspike+1;
            end

            % Euler update to V
            V=V+(dt/Cm)*(Istim(i)*Jstim+Ie+Ii+Ix+gL*(EL-V)+gL*DeltaT*exp((V-VT)/DeltaT));

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

                % Update synaptic currents
                Ie=Ie+sum(J(:,Ispike(Ispike<=Ne)),2)/taue;    
                Ii=Ii+sum(J(:,Ispike(Ispike>Ne)),2)/taui;               

                % Update cumulative number of spikes
                nspike=nspike+numel(Ispike);
            end            

            % Euler update to synaptic currents
            Ie=Ie-dt*Ie/taue;
            Ii=Ii-dt*Ii/taui;
            Ix=Ix-dt*Ix/taux;

            % This makes plots of V(t) look better.
            % All action potentials reach Vth exactly. 
            % This has no real effect on the network sims
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

        reSim=hist(s(2,s(1,:)>Tburn & s(2,:)<=Ne),1:Ne)/(T-Tburn);
        riSim=hist(s(2,s(1,:)>Tburn & s(2,:)>Ne)-Ne,1:Ni)/(T-Tburn);

        % Mean rate over E and I pops
        mean_e_rate(k,thr_indx)=mean(reSim);
        mean_i_rate(k,thr_indx)=mean(riSim);

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

            mRee(k,thr_indx)=nanmean(R(II<=Ne & JJ<=II & isfinite(R)));
            mRei(k,thr_indx)=nanmean(R(II<=Ne & JJ>Ne & isfinite(R)));
            mRii(k,thr_indx)=nanmean(R(II>Ne & JJ>II & isfinite(R)));
        end
    end
end
%% Plot mean f-I curve
figure; hold on
for thr_indx = 1:num_points_gain
%     plot(jestim_vector,1000*mean_e_rate,'linewidth',3,'color','blue')
%     plot(jestim_vector,1000*mean_i_rate,'linewidth',3,'color','red')
    plot(rx_vector,1000*(0.8*mean_e_rate(:,thr_indx)+0.2*mean_i_rate(:,thr_indx)),'linewidth',3)
end
xlabel('Input')
ylabel('Frequency')

%% Plot mean corrs
figure; hold on
for thr_indx = 1:num_points_gain
%     plot(jestim_vector,mRee,'color','blue','linewidth',3)
%     plot(jestim_vector,mRei,'color','green','linewidth',3)
%     plot(jestim_vector,mRii,'color','red','linewidth',3)
    plot(rx_vector,(mRee(:,thr_indx)+2*mRei(:,thr_indx)+mRii(:,thr_indx))/4,'linewidth',3)
end
xlabel('Input')
ylabel('Correlations')


save('./f_i_curve_network_DT.mat', 'c','Cm','DeltaT','delta_vector',...
    'dt','EL','gL','jestim_vector','Jm','Jxm','mean_e_rate',...
    'mean_i_rate','mRee','mRei',...
    'mRii','N','Ne','Ni','num_points_input','num_points_gain',...
    'Nx','P','rx_vector','T','T1','T2',...
    'Vre','VT','Vth','W','winsize','Wx');