function s = attention_pulse(fs,T) 

%  Input pulse to the network
% fs             % sample rate [Hz]
% T              % total duration [s]

t  = (0 : 1/fs : T).';   

% Event descriptions
ev(1) = struct('t0',2,'tp',2.2,'te',4,'A',1);   % s1
ev(2) = struct('t0',3,'tp',3.5,'te',4,'A',1);   % s2
ev(3) = struct('t0',rand(1)/5,'tp',(rand(1)+1)/2,'te',rand(1)+10,'A',1);   % s2

double_exp = @(e,tt) ...
    arrayfun(@(x) transient_sample(e,x), tt);

function y = transient_sample(e,tt)
    if tt < e.t0           % before onset
        y = 0;
    elseif tt < e.tp       % rising phase
        tau_r = (e.tp - e.t0)/4;               
        rise  = 1 - exp(-(tt - e.t0)/tau_r);
        scale = 1 / (1 - exp(-(e.tp - e.t0)/tau_r));  % normalise to 1 at tp
        y = e.A * rise * scale;
    else                   % decay phase
        tau_d = (e.te - e.tp)/log(20);           % ensures 5 % at te
        y = e.A * exp(-(tt - e.tp)/tau_d);
    end
end

s1 = double_exp(ev(1), t)/3;
s2 = double_exp(ev(2), t)/2.5;
s4 = double_exp(ev(3), t)/2;
s3 = s1 + s2;

s = [s1';s2';s3';s4'];
end