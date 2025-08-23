%% Parallel Spiking network modeling and dynamic routing - Four Networks Multiple Trials

% This script simulates the biophysical interplay of four spiking networks
% Networks 1-2 are reciprocally connected
% Network 1 has feedforward connection to Network 4
% Network 2 has feedforward connection to Network 3
% Modified for parallel processing with 100 trials

clear 
close all

% Set up parallel pool with 3 workers
if isempty(gcp('nocreate'))
    parpool(3);
end

% Number of trials
num_trials = 100;

% Pre-allocate output arrays
LFP_all = cell(num_trials, 4);          % LFP for each network
spike_density_all = cell(num_trials, 4); % Spike density for each network  
burst_density_all = cell(num_trials, 4); % Burst density for each network

% Run parallel simulations
parfor trial_idx = 1:num_trials
    fprintf('Running trial %d/%d\n', trial_idx, num_trials);
    
    % Run single trial and get outputs
    [LFP_trial, spike_density_trial, burst_density_trial] = run_single_trial();
    
    % Store results
    for net_idx = 1:4
        LFP_all{trial_idx, net_idx} = LFP_trial{net_idx};
        spike_density_all{trial_idx, net_idx} = spike_density_trial{net_idx};
        burst_density_all{trial_idx, net_idx} = burst_density_trial{net_idx};
    end
end

% Save results
save('simulation_results_Invalid_03.mat', 'LFP_all', 'spike_density_all', 'burst_density_all', '-v7.3');
fprintf('Results saved to simulation_results.mat\n');

% Function to run a single trial
function [LFP_out, spike_density_out, burst_density_out] = run_single_trial()
    
    Net.NumNets = 4; % number of networks
    Net.NumCells = 1000; % number of neurons per network
    Ntot = Net.NumNets*Net.NumCells;
    % number of excitatory and inhibitory neurons
    Net.EI.prop{1,1} = [.80 .20];  % ratio for network 1
    Net.EI.prop{2,1} = [.80 .20];  % ratio for network 2
    Net.EI.prop{3,1} = [.80 .20];  % ratio for network 3
    Net.EI.prop{4,1} = [.80 .20];  % ratio for network 4

    Net.EItypes.prop{1,1} = [.7 .2 .1];
    Net.EItypes.prop{1,2} = [.25 .25 .25 .25]; % prop for PV, CCK, CB, CR in Network 1

    Net.EItypes.prop{2,1} = [.7 .2 .1];
    Net.EItypes.prop{2,2} = [.25 .25 .25 .25]; % prop for PV, CCK, CB, CR in Network 2

    Net.EItypes.prop{3,1} = [.7 .2 .1];
    Net.EItypes.prop{3,2} = [.25 .25 .25 .25]; % prop for PV, CCK, CB, CR in Network 3

    Net.EItypes.prop{4,1} = [.7 .2 .1];
    Net.EItypes.prop{4,2} = [.25 .25 .25 .25]; % prop for PV, CCK, CB, CR in Network 4

    Fs = 2E3;
    tau = 1E3/Fs;
    trialtime = 4000;
    W11 = 1;
    W22 = 1;
    W33 = 1;
    W44 = 1;
    Nmua = 125;

    % proportion of internetwork projection type
    inter_E = .80;
    inter_I = .20;

    % NEW: Percentage of neurons receiving inter-network inputs
    inter_receive_ratio = 1; % 30% of neurons receive inter-network inputs

    % Inter-network synaptic delay (in ms)
    inter_delay_ms = 10; % 10ms delay
    inter_delay_samples = round(inter_delay_ms / tau); % Convert to samples

    WIA = .1; % inter-network connection ratio
    W12 = WIA; % Network 1 to 2
    W21 = WIA; % Network 2 to 1
    W14 = WIA; % Network 1 to 4 (feedforward)
    W23 = WIA; % Network 2 to 3 (feedforward)
    W34 = WIA; % Network 3 to 4
    W43 = WIA; % Network 4 to 3

    w12 = 1;
    w21 = 1;
    w14 = .75;
    w23 = .75;
    w41 = .1;
    w32 = .1;
    w34 = 1;
    w43 = 1;

    load('WF_mat.mat')
    % Ensure WF_mat matches total neuron count
    if size(WF_mat,1) < Ntot
        % Repeat WF_mat if needed
        WF_mat = repmat(WF_mat, ceil(Ntot/size(WF_mat,1)), 1);
    end
    WF_mat = WF_mat(1:Ntot,901-20:5:901+20);

    %% Making a cubic structure for each network with recording sites
    Net_Label = {'Network1','Network2','Network3','Network4'};

    NNs = Net.NumNets; % number of networks
    N = Net.NumCells; % number of neurons per network
    Ntot = NNs*N; % total number of neurons in the simulation

    Np = nthroot(N,3);
    Sp = round(nthroot(N,3)/2);
    Ep = (Np-1); % end neuron point
    Nd = -Ep:2:Ep;
    Ed = 2*(-Sp+1:2:Sp-1);
    SpN = length(Ed);

    [nx,~,~] = meshgrid(Nd);
    [ex,~,~] = meshgrid(Ed);

    % current nodes
    nx = repmat(Nd,Np.^2,1);
    nx = reshape(nx,numel(nx),1);
    ny = meshgrid(Nd);
    ny = repmat(reshape(ny,numel(ny),1),Np,1);
    nz = repmat(Nd',Np.^2,1);
    Nc = [nx,ny,nz];

    % recording nodes (4 disks for 4 networks)
    ex = 0;
    ey = 0;
    ez = 15;
    Ec = [ex,ey,ez];
    NE = length(ex);

    % pair wise distance of R nodes to C nodes
    dE = zeros(NE,N);
    dNs = [];
    for i=1:NE
        dE(i,:) = abs(Ec(i,3)-Nc(:,3));
    end
    for i=1:N
        dNs(i,:) = sqrt(sum((Nc(i,:)-Nc).^2,2));
    end
    gC = exp(-1*(dNs/std(dNs,[],'all')).^2)+diag(-1*ones(1,N));

    gC(isinf(gC)) = 0;
    gC = 1*gC./max(gC(gC>0),[],'all');

    [rn,~,r_ind] = unique(dE);
    r_ind = reshape(r_ind,NE,N);

    cntr1 = 0;
    cntr2 = 0;

    Level1_name = []; 
    Level2_name = [];
    Level3_name = [];

    Net.Level1.Net_ID = [1:NNs]';
    Net.Level1.Net_Size = repmat(N,NNs,1);
    Net.Level1.Net_Ind = [N:N:Ntot]';
    Level1_name = fields(Net.Level1);
    Net.Level1.Net_Spec = table(Net.Level1.Net_ID,Net.Level1.Net_Size,Net.Level1.Net_Ind,'VariableNames',Level1_name);
    Net.Level1.Net_Mat = [Net.Level1.Net_ID,Net.Level1.Net_Size,Net.Level1.Net_Ind];
    Net.Level1.Net_Mat = [zeros(1,size(Net.Level1.Net_Mat,2));Net.Level1.Net_Mat];

    for i=1:NNs % network #
        for j=1:2 % E/I type
            cntr1 = cntr1+1;
            Net.Level2.Net_ID(cntr1,1) = i;
            Net.Level2.Net_Size(cntr1,1) = Net.Level1.Net_Size(i);
            Net.Level2.Net_Ind(cntr1,1) = Net.Level1.Net_Ind(i);
            Net.Level2.EI_ID(cntr1,1) = j;
            Net.Level2.EI_Size(cntr1,1) = Net.EI.prop{i,1}(j) * Net.Level1.Net_Size(i);
            if cntr1>1 
                Net.Level2.EI_Ind(cntr1,1) = Net.Level2.EI_Ind(cntr1-1,1) + Net.Level2.EI_Size(cntr1,1);
            else
                Net.Level2.EI_Ind(cntr1,1) = Net.Level2.EI_Size(cntr1,1);
            end
                      
            for k=1:length(Net.EItypes.prop{i,j}) % E/I subtypes
                cntr2 = cntr2+1;
                Net.Level3.Net_ID(cntr2,1) = i;
                Net.Level3.Net_Size(cntr2,1) = Net.Level1.Net_Size(i);
                Net.Level3.Net_Ind(cntr2,1) = Net.Level1.Net_Ind(i);
                Net.Level3.EI_ID(cntr2,1) = j;
                Net.Level3.EI_Size(cntr2,1) = Net.Level2.EI_Size(cntr1,1);
                Net.Level3.EI_Ind(cntr2,1) = Net.Level2.EI_Ind(cntr1,1);
                Net.Level3.EItypes_ID(cntr2,1) = k;
                
                % Calculate subtype size, adjusting last one for rounding
                if k < length(Net.EItypes.prop{i,j})
                    Net.Level3.EItypes_Size(cntr2,1) = round(Net.EItypes.prop{i,j}(k) * Net.Level2.EI_Size(cntr1,1));
                else
                    % Last subtype gets remainder to ensure exact total
                    prev_sum = 0;
                    if cntr2 > 1
                        % Find start of current E/I type
                        start_idx = find(Net.Level3.EI_ID == j & Net.Level3.Net_ID == i, 1);
                        if ~isempty(start_idx) && start_idx < cntr2
                            prev_sum = sum(Net.Level3.EItypes_Size(start_idx:cntr2-1));
                        end
                    end
                    Net.Level3.EItypes_Size(cntr2,1) = Net.Level2.EI_Size(cntr1,1) - prev_sum;
                end
                
                if cntr2>1 
                    Net.Level3.EItypes_Ind(cntr2,1) = Net.Level3.EItypes_Ind(cntr2-1,1) + Net.Level3.EItypes_Size(cntr2,1);
                else
                    Net.Level3.EItypes_Ind(cntr2,1) = Net.Level3.EItypes_Size(cntr2,1);
                end
            end
        end
        Net.rand_Ind{i} = randperm(N); % randomly assigning each neuron to an index of the geometric structure
        [~,Net.rand_IndS2C{i}] = sort(Net.rand_Ind{i}); % index of the structure element assigned to the Neuron ID
    end

    Level2_name = fields(Net.Level2);
    Net.Level2.EI_Spec = table(Net.Level2.Net_ID,Net.Level2.Net_Size,Net.Level2.Net_Ind,Net.Level2.EI_ID,...
       Net.Level2.EI_Size,Net.Level2.EI_Ind,'VariableNames',Level2_name);
    Net.Level2.Net_Mat = [Net.Level2.Net_ID,Net.Level2.Net_Size,Net.Level2.Net_Ind,...
       Net.Level2.EI_ID,Net.Level2.EI_Size,Net.Level2.EI_Ind];
    Net.Level2.Net_Mat = [zeros(1,size(Net.Level2.Net_Mat,2));Net.Level2.Net_Mat];
              
    Level3_name = fields(Net.Level3);
    Net.Level3.EI_Spec = table(Net.Level3.Net_ID,Net.Level3.Net_Size,Net.Level3.Net_Ind,Net.Level3.EI_ID,...
                          Net.Level3.EI_Size,Net.Level3.EI_Ind,Net.Level3.EItypes_ID,Net.Level3.EItypes_Size,Net.Level3.EItypes_Ind,...
                          'VariableNames',Level3_name);
    Net.Level3.Net_Mat = [Net.Level3.Net_ID,Net.Level3.Net_Size,Net.Level3.Net_Ind,Net.Level3.EI_ID,...
                          Net.Level3.EI_Size,Net.Level3.EI_Ind,Net.Level3.EItypes_ID,Net.Level3.EItypes_Size,Net.Level3.EItypes_Ind];
    Net.Level3.Net_Mat = [zeros(1,size(Net.Level3.Net_Mat,2)); Net.Level3.Net_Mat];

    Wcon = nan(Ntot);
    l_ind = [0;Net.Level1.Net_Ind];
    for i=1:NNs
        Wcon(l_ind(i)+1:l_ind(i+1),l_ind(i)+1:l_ind(i+1)) = gC(Net.rand_Ind{i},Net.rand_Ind{i});
    end
    GL = Wcon;
    GL(GL==0) = 1;
    GL(isnan(GL)) = 0;
    Wcon(isnan(Wcon)) = 0.1; % baseline inter-network connection probability

    % Calculate network indices
    net_idx = zeros(NNs+1,1);
    for i=1:NNs+1
        net_idx(i) = (i-1)*N;
    end

    %% Create inter-network connectivity masks
    % Initialize masks for neurons that can receive inter-network inputs
    inter_receive_mask = zeros(Ntot, 1);
    for i = 1:NNs
        % Select random neurons in each network to receive inter-network inputs
        net_neurons = net_idx(i)+1:net_idx(i+1);
        n_inter_receivers = round(N * inter_receive_ratio);
        selected_neurons = randperm(N, n_inter_receivers) + net_idx(i);
        inter_receive_mask(selected_neurons) = 1;
    end

    % Create a mask matrix for inter-network connections
    inter_conn_mask = zeros(Ntot, Ntot);
    for src = 1:NNs
        for tgt = 1:NNs
            if src ~= tgt
                src_neurons = net_idx(src)+1:net_idx(src+1);
                tgt_neurons = net_idx(tgt)+1:net_idx(tgt+1);
                % Only neurons marked to receive inter-network inputs can have connections
                tgt_receivers = tgt_neurons(inter_receive_mask(tgt_neurons) == 1);
                inter_conn_mask(src_neurons, tgt_receivers) = 1;
            end
        end
    end

    %% main connectivity 
    WEE = [W11 w12 0   w14;   % Net1->1,2,4
           w21 W22 w23 0;     % Net2->1,2,3
           0   w32   W33 w34;   % Net3->3,4
           w41   0   w43 W44];  % Net4->3,4

    WEI = [W11 w12 0   w14;   % E->I connections
           w21 W22 w23 0;
           0 0   W33 w34;
           0 0   w43 W44];

    WIE = [W11 w12 0   0;   % E->I connections
           w21 W22 0   0;
           0   0   W33 w34;
           0   0   w43 W44];

    WII = [W11 w12 0   0;   % E->I connections
           w21 W22 0   0;
           0   0   W33 w34;
           0   0   w43 W44];

    Net = set_initial_network_iEEG(WEE,WEI,WIE,WII,Net);
    Net.C = Net.C/2;

    l_ind = [0;Net.Level3.EItypes_Ind];
    l_size = Net.Level3.EItypes_Size; 
    l_ID = Net.Level3.EI_ID;  
    [I,J] = size(Net.C);

    % connectivity mat
    for i = 1:I
        for j = 1:J
            PSCSi = ((-1).^(l_ID(i)-1));
            PSCSj = ((-1).^(l_ID(i)-1));
            Net.Cmat(l_ind(i)+1:l_ind(i+1),l_ind(j)+1:l_ind(j+1)) = Net.C(i,j).*(.9+.1*rand(l_size(i),l_size(j))).*PSCSi;         
        end
    end

    %% parameter settings 
    SPK_time = [];             % spike-times
    tspan = 0:tau:trialtime;
    Nsample = length(tspan);
    SPK_time = zeros(Ntot,Nsample);
    SPK_timeZ = SPK_time;
    SPK_Ones = ones(Ntot,Nsample);
    ext_Input = zeros(Ntot,Nsample);

    for ipar=1:9
        for j=1:J
            if ipar<9
                Net.ParMat(l_ind(j)+1:l_ind(j+1),ipar) = Net.par(ipar,j).* ones(l_size(j),1)+Net.par(ipar+9,j).*(.5-rand(l_size(j),1));
            else
                ext_Input(l_ind(j)+1:l_ind(j+1),:) = Net.par(ipar,j).* ones(l_size(j),Nsample)+.2*Net.par(ipar+9,j).*(.5-rand(l_size(j),Nsample)); 
            end
        end
    end

    a = Net.ParMat(:,1);
    b = Net.ParMat(:,2);
    c = Net.ParMat(:,3);
    d = Net.ParMat(:,4);
    tr = Net.ParMat(:,5);
    td = Net.ParMat(:,6);
    Syn_amp = Net.ParMat(:,7);
    v = Net.ParMat(:,8);

    cn = dsp.ColoredNoise(2,'SamplesPerFrame',Nsample,...
        'NumChannels',Ntot);
     
    Bnoise = cn()';
    Bnoise = Bnoise./max(abs(Bnoise));
    Bnoise = Bnoise(:,1:Nsample);

    Syn_cn = dsp.ColoredNoise(2,'SamplesPerFrame',Nsample,...
        'NumChannels',Ntot);
    NG_syn = Syn_cn()';
    NG_syn = NG_syn./max(abs(NG_syn));
    NG_syn = (NG_syn/50);
      
    %% Simulation 
    C = Net.Cmat.*Wcon;
    PC = rand(size(C))<abs(C);

    % Apply inter-network connectivity mask
    C = C .* (1 - inter_conn_mask) + C .* inter_conn_mask;

    % Define network ranges dynamically
    IA = cell(NNs,1);
    for i=1:NNs
        IA{i} = net_idx(i)+1:net_idx(i+1);
    end
    IA1 = IA{1}; IA2 = IA{2}; IA3 = IA{3}; IA4 = IA{4};

    C = sign(C).*PC;

    %% External Input
    t1s = linspace(-2*pi,2*pi,Nsample/50);
    t2s = linspace(-pi,pi,Fs/2);
    x1 = square(t2s,70)+2;
    x2 = -square(t2s,70)+2;
    GWin = gausswin(500./tau,4)';
    GWin = [GWin];
    Input_pulse1 = [];
    Input_pulse2 = [];
    Input_pulse3 = [];
    Input_pulse4 = [];

    % desynchronized input ratio
    i1 = [99, 99];
    i2 = [99, 99];
    i3 = [99, 99];
    i4 = [99, 99];

    InputL = length(i1);

    for i=1:InputL
        in1 = InputFun_02(Ntot/4, (Nsample-1)/InputL, Fs, 1, i1(i));
        in2 = InputFun_02(Ntot/4, (Nsample-1)/InputL, Fs, 1, i2(i));
        in3 = InputFun_02(Ntot/4, (Nsample-1)/InputL, Fs, 1, i3(i));
        in4 = InputFun_02(Ntot/4, (Nsample-1)/InputL, Fs, 1, i4(i));
        if i==InputL
            in1 = InputFun_02(Ntot/4, (Nsample-1)/InputL+1, Fs, 1, i1(i));
            in2 = InputFun_02(Ntot/4, (Nsample-1)/InputL+1, Fs, 1, i2(i));
            in3 = InputFun_02(Ntot/4, (Nsample-1)/InputL+1, Fs, 1, i3(i));
            in4 = InputFun_02(Ntot/4, (Nsample-1)/InputL+1, Fs, 1, i4(i));
        end
        Input_pulse1 = [Input_pulse1,in1];
        Input_pulse2 = [Input_pulse2,in2];
        Input_pulse3 = [Input_pulse3,in3];
        Input_pulse4 = [Input_pulse4,in4];
    end

    coherence = 0.01;
    freqShift = 0;
    CosMat1 = Input_Cosine(Ntot/4,Fs,trialtime/1000,.01,rand(1)/100,freqShift);
    CosMat2 = Input_Cosine(Ntot/4,Fs,trialtime/1000,.01,rand(1)/100,freqShift);
    CosMat3 = Input_Cosine(Ntot/4,Fs,trialtime/1000,.01,rand(1)/100,freqShift);
    CosMat4 = Input_Cosine(Ntot/4,Fs,trialtime/1000,.01,rand(1)/100,freqShift);

    s = attention_pulse(Fs,round(trialtime/1000));

    CosMat1 = round((CosMat1+1)./2);
    CosMat2 = round((CosMat2+1)./2);
    CosMat3 = round((CosMat3+1)./2);
    CosMat4 = round((CosMat4+1)./2);

    Input_p1 = ones(size(CosMat1));
    Input_p2 = ones(size(CosMat2));
    Input_p3 = ones(size(CosMat3));
    Input_p4 = ones(size(CosMat4));

    Zs1 = randperm(N,N/2+round(N*rand(1)/2));
    Zs2 = randperm(N,N/2+round(N*rand(1)/2));
    Zs3 = randperm(N,N/2+round(N*rand(1)/2));
    Zs4 = randperm(N,N/2+round(N*rand(1)/2));

    Input_p1(Zs1,:) = 0;
    Input_p2(Zs2,:) = 0;
    Input_p3(Zs3,:) = 0;
    Input_p4(Zs4,:) = 0;

    Input_pulseC = [CosMat1;CosMat2;CosMat3;CosMat4];
    Input_pulse = [Input_pulse1;Input_pulse2;Input_pulse3;Input_pulse4];
    Input_pulse = conv2(Input_pulse,GWin,'same');
    Input_pulse(Input_pulse<.5) = 0;
    Input_pulse(Input_pulse>=.5) = 1;
    Input_pulse = Input_pulse./max([Input_pulse,ones(Ntot,1)],[],'all');

    Input_pulse(IA{1},:) = Input_pulseC(IA{1},:)+1*(s(3,:)).*Input_pulseC(IA{1},:)+2*(s(4,:)).*Input_pulseC(IA{1},:);
    Input_pulse(IA{2},:) = Input_pulseC(IA{2},:)+1*(s(1,:)).*Input_pulseC(IA{1},:)+2*(s(4,:)).*Input_pulseC(IA{1},:);
    Input_pulse(IA{3},:) = Input_pulseC(IA{3},:)+1*(s(2,:)).*Input_pulseC(IA{3},:)+2*(s(4,:)).*Input_pulseC(IA{1},:);
    Input_pulse(IA{4},:) = Input_pulseC(IA{4},:)+1*(s(2,:)).*Input_pulseC(IA{4},:)+2*(s(4,:)).*Input_pulseC(IA{1},:);

    %% Run simulation
    vv = [];
    VVrmv = [];
    Nspiked = [];
    u = b.*v;
    AttenL = 25*Fs/1000;
    Atten = [1-hamming(2*AttenL+1)]';

    % NEW: Initialize delay buffers for inter-network connections
    % Create separate buffers for local and delayed inputs
    G_syn_local = zeros(Ntot, 1);
    G_syn_delayed = zeros(Ntot, inter_delay_samples);
    Gout = zeros(Ntot, 1);
    delay_idx = 1; % Current position in circular buffer

    % initialization of synaptic parameters
    T = 1E6*ones(Ntot,1);
    v_syn = -(-(1./tr).*exp((-T-tau)./tr)+(1./td).*exp((-T-tau)./td));
    G_syn = zeros(Ntot,1);
    tm = log(tr./td)./(1./td -1./tr);
    Mm = -exp(-tm./tr)+exp(-tm./td);

    % initialization function of Synaptic voltage gradient
    init_V = @(Tr,Td,Ts) -(-(1./Tr).*exp((-Ts+tau)./Tr)+(1./Td).*exp((-Ts+tau)./Td));

    % Separate connectivity matrices for local and inter-network connections
    C_local = C .* (1 - inter_conn_mask);
    C_inter = C .* inter_conn_mask;

    for t=1:Nsample
        % Get delayed synaptic input from the circular buffer
        delayed_input = G_syn_delayed(:, delay_idx);
        
        % Combine local and delayed inputs
        Gout_combined = G_syn_local + delayed_input;
        
        I = (ext_Input(:,t)-1).*(Input_pulse(:,t)) + (Bnoise(:,t));
        
        Nspiked = find(v > 30);
        SPK_time(Nspiked,t) = 1;

        % Calculate synaptic input with separated local and delayed components
        I = I + C_local'*((G_syn_local + NG_syn(:,t)).*Syn_amp)*2 + ...
                C_inter'*((delayed_input + NG_syn(:,t)).*Syn_amp)*2;
        
        v = v + tau*(0.04*v.^2+5*v+140-u+I);
        u = u + tau*a.*(b.*v-u);
        v(Nspiked) = c(Nspiked);
        u(Nspiked) = u(Nspiked) + d(Nspiked);

        % Update synaptic variables
        T(Nspiked) = tr(Nspiked).*(td(Nspiked)./(tr(Nspiked)+td(Nspiked))).*(G_syn_local(Nspiked)/(.63));
        v_syn(Nspiked) = init_V(tr(Nspiked),td(Nspiked),T(Nspiked));
        v_syn = v_syn+tau*(-v_syn.*((1./td)+(1./tr)) - ((1./(td.*tr)).*G_syn));
        G_syn = G_syn+v_syn.*tau;
        G_syn_local = G_syn./Mm;
        
        % Store current synaptic output in delay buffer for inter-network connections
        G_syn_delayed(:, delay_idx) = G_syn_local;
        
        % Update circular buffer index
        delay_idx = mod(delay_idx, inter_delay_samples) + 1;
        
        GV(:,t) = Gout_combined;
        
        VVrmv = [VVrmv,v];
        vv = [vv,v];
        
        mInd = min([AttenL,t-1]);
        MInd = min([AttenL,Nsample-t]);    
        SPK_timeZ(Nspiked,t-mInd:t+MInd) = 1;
        SPK_Ones(Nspiked,t-mInd:t+MInd) = SPK_Ones(Nspiked,t-mInd:t+MInd).*Atten(AttenL+1-mInd:AttenL+1+MInd);
    end

    %% Calculate LFP
    l_ind = [0;Net.Level3.EItypes_Ind];  
    l_EItype_ID = Net.Level3.EItypes_ID;
    l_EI_ID = Net.Level3.EI_ID;
    l_Net_ID = Net.Level3.Net_ID;  
    l_Net_Ind = [0;Net.Level1.Net_Ind];

    r_ind1 = r_ind(:,Net.rand_IndS2C{1});
    r_ind2 = r_ind(:,Net.rand_IndS2C{2}); 
    r_ind3 = r_ind(:,Net.rand_IndS2C{3});
    r_ind4 = r_ind(:,Net.rand_IndS2C{4});

    Vout0 = -(VVrmv+70);
    Vout = Vout0;
    A = [];
    for i=1:Ntot
        A(i,:) = conv(SPK_time(i,:),WF_mat(i,:),'same');
    end
    sd = std(Vout0,[],2);
    B = 100.*A;

    PinkN = dsp.ColoredNoise(2,'SamplesPerFrame',Nsample,...
        'NumChannels',Ntot);
    CC = PinkN()';
    cn = dsp.ColoredNoise(2,'SamplesPerFrame',Nsample,...
        'NumChannels',Ntot);

    CC = CC(1:Ntot,1:Nsample);
    CC = CC+cn()';
    normCC = (CC./range(CC,2));
    sdN = nanstd(normCC,[],2);

    Spkwin = find(SPK_timeZ);
    Vout(Spkwin) = Vout0(Spkwin).*SPK_Ones(Spkwin);
    Vout = Vout0;

    % Ohmic LFP signals on the recording sites
    Vr1 = (10+rn(r_ind1)./min(rn(r_ind1),[],2));
    Vr2 = (10+rn(r_ind2)./min(rn(r_ind2),[],2));
    Vr3 = (10+rn(r_ind3)./min(rn(r_ind3),[],2));
    Vr4 = (10+rn(r_ind4)./min(rn(r_ind4),[],2));

    idx1 = Vr1<18;
    idx2 = Vr2<18;
    idx3 = Vr3<18;
    idx4 = Vr4<18;

    % Ohmic Monopole
    mVEs1 = ((1./Vr1).^1)* (Vout(IA{1},:)+(sign(nansum(C(IA{1},:),2)).*GV(IA{1},:)));
    mVEs2 = ((1./Vr2).^1)* (Vout(IA{2},:)+(sign(nansum(C(IA{2},:),2)).*GV(IA{2},:)));
    mVEs3 = ((1./Vr3).^1)* (Vout(IA{3},:)+(sign(nansum(C(IA{3},:),2)).*GV(IA{3},:)));
    mVEs4 = ((1./Vr4).^1)* (Vout(IA{4},:)+(sign(nansum(C(IA{4},:),2)).*GV(IA{4},:)));

    % Ohmic Dipole
    dVEs1 = ((1./Vr1(idx1)).^2)* Vout(IA{1}(idx1),:) + ((1./Vr1(idx1)).^1)*(sign(nansum(C(IA{1}(idx1),:),2)).*GV(IA{1}(idx1),:));
    dVEs2 = ((1./Vr2(idx2)).^2)* Vout(IA{2}(idx2),:) + ((1./Vr2(idx2)).^1)*(sign(nansum(C(IA{2}(idx2),:),2)).*GV(IA{2}(idx2),:));
    dVEs3 = ((1./Vr3(idx3)).^2)* Vout(IA{3}(idx3),:) + ((1./Vr3(idx3)).^1)*(sign(nansum(C(IA{3}(idx3),:),2)).*GV(IA{3}(idx3),:));
    dVEs4 = ((1./Vr4(idx4)).^2)* Vout(IA{4}(idx4),:) + ((1./Vr4(idx4)).^1)*(sign(nansum(C(IA{4}(idx4),:),2)).*GV(IA{4}(idx4),:));

    %% Adaptive Non-Ohmic Filter
    fL          = 10;        % [Hz] upper bound of the "low-frequency" band
    alpha_fast  = 2e-4;       % steep roll-off when low-freq power is weak
    alpha_slow  = 3e-3;       % gentle roll-off when low-freq power is strong

    signals_in  = {dVEs1, dVEs2, dVEs3, dVEs4};
    signals_out = cell(size(signals_in));

    for k = 1:numel(signals_in)

        x = signals_in{k}(:).';                 % ensure row-vector
        padLength  = round(numel(x) / 4);       % mirror padding length

        % mirror-pad, FFT
        x_padded   = [fliplr(x(1:padLength)), x, fliplr(x(end-padLength+1:end))];
        X          = fft(x_padded);

        % frequency axis
        n_padded   = numel(x_padded);
        f          = (0:n_padded-1) * (Fs / n_padded);       % 0 … Fs(1-1/n)

        % low-band / total power ratio
        mag2       = abs(X).^2;
        idxLow     = f <= fL | f >= Fs - fL;                 % two-sided FFT
        rho        = sum(mag2(idxLow)) / sum(mag2);          % 0 … 1

        % build adaptive attenuation curve
        attenuation = (1-rho) .* exp(-alpha_fast .* f) + ...
                        rho  .* exp(-alpha_slow .* f);

        % filter & inverse FFT
        X_filt     = X .* attenuation;
        x_filt_pad = real(ifft(X_filt));

        % remove padding
        x_filt     = x_filt_pad(padLength+1 : end-padLength);

        % optional start-point stabilisation
        if numel(x_filt) >= 100
            x_filt(1:100) = x_filt(100);
        end

        % collect
        signals_out{k} = x_filt;
    end

    [filtered_dVEs1, filtered_dVEs2, filtered_dVEs3, filtered_dVEs4] = ...
        deal(signals_out{:});

    %% Detect bursts
    freq_band = [65 115]; 
    Num_cycle = 1.5;
    separating_cycle = 1;

    out1 = detect_LFP_burst(filtered_dVEs1, freq_band, Num_cycle, separating_cycle, Fs, 'DensityWindow', 0.25,'BurstDurationCycles', 1.5);
    out2 = detect_LFP_burst(filtered_dVEs2, freq_band, Num_cycle, separating_cycle, Fs, 'DensityWindow', 0.25,'BurstDurationCycles', 1.5);
    out3 = detect_LFP_burst(filtered_dVEs3, freq_band, Num_cycle, separating_cycle, Fs, 'DensityWindow', 0.25,'BurstDurationCycles', 1.5);
    out4 = detect_LFP_burst(filtered_dVEs4, freq_band, Num_cycle, separating_cycle, Fs, 'DensityWindow', 0.25,'BurstDurationCycles', 1.5);

    %% Prepare outputs
    LFP_out = {filtered_dVEs1, filtered_dVEs2, filtered_dVEs3, filtered_dVEs4};
    
    spike_density_out = cell(4,1);
    for i = 1:4
        spike_density_out{i} = sum(SPK_time(IA{i},:));
    end
    
    burst_density_out = {out1.burst_density, out2.burst_density, out3.burst_density, out4.burst_density};
end