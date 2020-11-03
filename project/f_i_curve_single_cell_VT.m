%% Analysis of correlation modulation in Balanced networks. 


clear

seed = 39;

% Number of neurons in ffwd layer
Nx=100;

% Ffwd connection probs
Px=[.1; .1];

% Correlation between the spike trains in the ffwd layer
c=0;
% Timescale of correlation
taujitter=5;

Jxm=[180]/sqrt(Nx);

% Time (in ms) for sim
T=500000;

% Time discretization
dt=.1; %ms

% FFwd spike train rate (in kHz)
rx=10/1000;

% Number of time bins
Nt=round(T/dt);
time=dt:dt:T;

% Synaptic timescales
taux=10;
taue=8;
taui=4;

%%% Make (correlated) Poisson spike times for ffwd layer
%%% See Appendix 1 of Akil et al 2020 for more details
tic
if(c<1e-5) % If uncorrelated
    nspikeX=poissrnd(Nx*rx*T);
    st=rand(nspikeX,1)*T;
    sx=zeros(2,numel(st));
    sx(1,:)=sort(st);
    sx(2,:)=randi(Nx,1,numel(st)); % neuron indices
    clear st;
else % If correlated
    rm=rx/c; % Firing rate of mother process
    nstm=poissrnd(rm*T); % Number of mother spikes
    stm=rand(nstm,1)*T; % spike times of mother process    
    maxnsx=T*rx*Nx*1.2; % Max num spikes
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

num_points_gain = 2;
threshold_vector = linspace(-65,-55,num_points_gain);

num_points_input = 16;
jestim_vector = linspace(-0.4,.1,num_points_input);
mean_rate = zeros(num_points_input,num_points_gain);

for k = 1:num_points_gain
    % Neuron parameters
    Cm=1;
    gL=1/15;
    EL=-72;
    Vth=-50;
    Vre=-75;
    DeltaT=1;
    VT=threshold_vector(k);
    parfor j=1:num_points_input

        rng(seed);
        % Preallocate memory
        Ix=zeros(1,1);

        % Set initial voltage
        % Random initial voltages
        V0=rand(1,1)*(VT-Vre)+Vre;
        V=V0;

        iFspike=1;
        s=zeros(2,maxns);
        nspike=0;
        TooManySpikes=0;

        J_NeNx = binornd(1,0.1,1,Nx);
        % Generate full connectivity matrices
        tic
        Jx=[J_NeNx.*Jxm(1)];

        tGen=toc;
        disp(sprintf('\nTime to generate connections: %.2f sec',tGen))

        % Extra stimulus: Istim is a time-dependent stimulus
        % it is delivered to all neurons with weights given by JIstim.
        % Specifically, the stimulus to neuron j at time index i is:
        % Istim(i)*JIstim(j)
        Istim=zeros(size(time)); 
        Istim(time>0)=1;
        jestim=jestim_vector(j); 
        Jstim=sqrt(Nx)*[jestim]; 

        tic
        for i=1:numel(time)
            % Propogate ffwd spikes
            while(sx(1,iFspike)<=time(i) && iFspike<nspikeX)
                jpre=sx(2,iFspike);
                Ix=Ix+Jx(:,jpre)/taux;
                iFspike=iFspike+1;
            end

            % Euler update to V
            V=V+(dt/Cm)*(Istim(i)*Jstim+Ix+gL*(EL-V)+gL*DeltaT*exp((V-VT)/DeltaT));

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

            if mod(i*dt,T/10) == 0 % print every x iterations.
                fprintf('At time %d...\n',i*dt/T);
            end
        end
        % Get rid of padding in s
        s=s(:,1:nspike); 
        tSim=toc;
        disp(sprintf('\nTime for simulation: %.2f min',tSim/60))
        Tburn=T/2;
        reSim=hist(s(2,s(1,:)>Tburn & s(2,:)<=1),1:1)/(T-Tburn);

        % Mean rate over E and I pops
        mean_rate(j,k)=mean(reSim);
        
        
    end
end

%% Plot mean f-I curve

figure; hold on
for k = 1:num_points_gain
    plot(jestim_vector,1000*mean_rate(:,k) , 'linewidth',3)
end


save('./f_i_curve_singlecell_VT.mat', 'c','Cm','DeltaT',...
    'dt','EL','gL','jestim_vector','Jxm','mean_rate',...
    'num_points_input','num_points_gain',...
    'Nx','Px','T','threshold_vector',...
    'Vre','VT','Vth');

