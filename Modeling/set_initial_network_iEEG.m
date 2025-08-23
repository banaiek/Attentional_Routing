function Net = set_initial_network_iEEG(WEE,WEI,WIE,WII,Net)

%   Net = SET_INITIAL_NETWORK_ECOG_V05(WEE, WEI, WIE, WII, Net)
%
%       WEE(rs,rt) multiplies *E→E* synapses
%       WEI(rs,rt) multiplies *E→I* synapses
%       WIE(rs,rt) multiplies *I→E* synapses
%       WII(rs,rt) multiplies *I→I* synapses


if nargin < 5 || isempty(Net)
    Net = struct();
end

if ~(isequal(size(WEE),size(WEI)) && isequal(size(WEE),size(WIE)) && isequal(size(WEE),size(WII)))
    error('All weight matrices must have identical dimensions.');
end
nNet   = size(WEE,1);
NcType = 7;            % 3 excitatory classes + 4 inhibitory classes

% intrinsic parameters
base_par = [...
   0.02 0.02 0.02  0.10 0.05 0.02 0.02 ;  % a0
   0.20 0.22 0.20  0.21 0.21 0.21 0.21 ;  % b0
  -65  -65  -65   -65  -65  -65  -65 ;    % c0
     8    4    2     2    2    2    2 ;   % d0
     1    1    1     5   5    5    5 ;   % synRise
   2.4  2.4  2.4   6.4 16.4 16.4 24.4 ;   % synDecay
     1    1    1     1    1    1    1 ;   % synAmp
  -65  -65  -65   -65  -65  -65  -65 ;    % v0
     3    3    3   3  3  3  3 ;   % extInputMean
  0.001 0.001 0.001 0.001 0.001 0.001 0.001 ;
   0.01 0.01 0.01 0.01 0.01 0.01 0.01 ;
   0.01 0.01 0.01 0.01 0.01 0.01 0.01 ;
   0.01 0.01 0.01 0.01 0.01 0.01 0.01 ;
   0.01 0.01 0.01 0.01 0.01 0.01 0.01 ;
   0.01 0.01 0.01 0.01 0.01 0.01 0.01 ;
   0.01 0.01 0.01 0.01 0.01 0.01 0.01 ;
   0.01 0.01 0.01 0.01 0.01 0.01 0.01 ;
     1    1    1   0.1  0.1  0.1  0.1 ];

Net.par = repmat(base_par,1,nNet);



Cintra_base = [ ...
 .5  .5   .5   0.5 0.5 0.5 0.5 ;  % E rows → E/I cols
 .5  .5   .5   0.5 0.5 0.5 0.5 ;
 .5  .5   .5   0.5 0.5 0.5 0.5 ;
  0.6 0.6 0.6 0.5 0.5 0.5 0.5 ;
  0.6 0.6 0.6 0.5 0.5 0.5 0.5 ;
  0.6 0.6 0.6 0.6 0.6 0.5 0.5 ;
  0.05 0.05 0.05 0.05 0.05 0.6 0.05 ];

% Feed‑forward base (src < tgt)
Cfwd_base = [ ...
 .5  .5   .5   0.5 0.5 0.5 0.5 ;  % E rows → E/I cols
 .5  .5   .5   0.5 0.5 0.5 0.5 ;
 .5  .5   .5   0.5 0.5 0.5 0.5 ;
  0.6 0.6 0.6 0.5 0.5 0.5 0.5 ;
  0.6 0.6 0.6 0.5 0.5 0.5 0.5 ;
  0.6 0.6 0.6 0.6 0.6 0.5 0.5 ;
  0.05 0.05 0.05 0.05 0.05 0.6 0.05 ];

% Feedback base (src > tgt)
Cbwd_base = [ ...
 .5  .5   .5   0.5 0.5 0.5 0.5 ;  % E rows → E/I cols
 .5  .5   .5   0.5 0.5 0.5 0.5 ;
 .5  .5   .5   0.5 0.5 0.5 0.5 ;
  0.6 0.6 0.6 0.5 0.5 0.5 0.5 ;
  0.6 0.6 0.6 0.5 0.5 0.5 0.5 ;
  0.6 0.6 0.6 0.6 0.6 0.5 0.5 ;
  0.05 0.05 0.05 0.05 0.05 0.6 0.05 ];

%  assemble the full connectivity matrix 
Net.C = zeros(NcType*nNet);

for src = 1:nNet
    idx_s = (src-1)*NcType + (1:NcType);
    for tgt = 1:nNet
        idx_t = (tgt-1)*NcType + (1:NcType);

        % select base block
        if src == tgt
            base_blk = Cintra_base;
        elseif src < tgt
            base_blk = Cfwd_base;
        else
            base_blk = Cbwd_base;
        end

        % weight multiplier matrix
        mult = ones(NcType);
        mult(1:3,1:3) = WEE(src,tgt);  % E→E
        mult(1:3,4:7) = WEI(src,tgt);  % E→I
        mult(4:7,1:3) = WIE(src,tgt);  % I→E
        mult(4:7,4:7) = WII(src,tgt);  % I→I

        Net.C(idx_s,idx_t) = base_blk .* mult;
    end
end

Net.nNet = nNet;

end
