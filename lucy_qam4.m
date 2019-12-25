% A skeleton BER script for a wireless link simulation
clear all;clc; close all

%% part a

% For the final version of this project, you must use these 3
% parameter. You will likely want to set numIter to 1 while you debug your
% link, and then increase it to get an average BER.
numIter = 20000;  % The number of iterations of the simulation
nSym = 1000;    % The number of symbols per packet
SNR_Vec = 0:2:16;
lenSNR = length(SNR_Vec);
m_ary = [2, 4, 16];        % The M-ary number, 2 corresponds to binary modulation

%chan = 1;          % No channel
chan = [1, 0.2, 0.4];
%%chan = [0.227 0.460 0.688 0.460 0.227]';   % Not so invertible, severe ISI

tic;

%we use 300 training bits in this case
trainlen = 300;
displayStr = ["BER-2 with ISI","BER-4 No ISI", "BER-16 No ISI"];

for it=1:length(m_ary)
    M = m_ary(it);
    berVec = zeros(numIter, lenSNR);
    parfor ii = 1:numIter

        msg = randi([0, M-1], nSym*(log2(M)), 1);     % Generate random bits
        % New bits must be generated at every
        % iteration

        % If you increase the M-ary number, as you most likely will, you'll need to
        % convert the bits to integers. See the BIN2DE function
        % For binary, our MSG signal is simply the bits
    
        % We reshape bits so that there are a proper number of bits per row,
        % Then we convert each row to decimal and move on.
        %msg = reshape(bits,[nSym, log2(m_ary)]);
        %msg = bi2de(msg,'left-msb');
        %msg = bits;
        bits = de2bi(msg, 'left-msb').'; %transpose here 
        bits = bits(:);

        for jj = 1:lenSNR % one iteration of the simulation at each SNR Value
            tx = qammod(msg,M);  % BPSK modulate the signal
            
            %if m_ary == 4:
            
            if M == 2
                if isequal(chan,1)
                    txChan = tx;
                    txNoisy = txChan;
                else
                    txChan = filter(chan,1,tx);  % Apply the channel.
                    txNoisy = awgn(txChan,SNR_Vec(jj)); % Add AWGN
                    
                    %equalizer
                    %lineq = comm.LinearEqualizer('Algorithm','LMS', 'NumTaps',6,'StepSize',0.01);
                    %p = lineq(txNoisy, tx(1:trainlen));
                    eq1 = lineareq(6, lms(0.01));
                    txNoisy = equalize(eq1,txNoisy,tx(1:trainlen)); % Equalize.
                    %txNoisy = filter(eq1.weights, 1, txNoisy);
                    reset(eq1);
                end
            else
                txNoisy = awgn(tx + (eps*1j), SNR_Vec(jj) + 10*log10(log2(M)),'measured');
                %channel = comm.AWGNChannel('NoiseMethod', ...
                %    'Signal to noise ratio (SNR)', 'SNR', SNR_Vec(jj));
                %txNoisy = channel(tx);
            end
            rx = qamdemod(txNoisy,M); %,'OutputType', 'integer'); % Demodulate
            %rxMSG = de2bi(rx, [], 2);
            
            % Again, if M was a larger number, I'd need to convert my symbols
            % back to bits here - convert each row to its binary sequence
            % the transpose and the rx(:) is housekeeping - conceptually we are
            % taking each row, appending it after the previous row, but we do
            % this transposed since we are working with columns
            rxTmp = de2bi(rx, 'left-msb').'; %transpose here 
            rxMSG = rxTmp(:);
            
            % Compute and store the BER for this iteration
            % We're interested in the BER, which is the 2nd output of BITERR
            [~, berVec(ii,jj)] = biterr(bits(trainlen+1:end), rxMSG(trainlen+1:end));  

        end  % End SNR iteration
    end      % End numIter iteration

    % Compute and plot the mean BER
    ber = mean(berVec,1);
    
    figure(it);
    semilogy(SNR_Vec, ber, 'DisplayName', displayStr(it))
    hold on;
    
    if M == 2
        berTheory2 = berawgn(SNR_Vec,'psk', 2,'nondiff');
        semilogy(SNR_Vec,berTheory2,'DisplayName', 'Theoretical BER for M=2')
        legend('Location', 'southwest')
    elseif M == 4
        berTheory4 = berawgn(SNR_Vec,'qam', 4,'nondiff');
        semilogy(SNR_Vec,berTheory4,'DisplayName', 'Theoretical BER for M=4')
        legend('Location', 'southwest')
    elseif M == 16
        berTheory16 = berawgn(SNR_Vec,'qam', 16,'nondiff');
        semilogy(SNR_Vec,berTheory16, 'DisplayName', 'Theoretical BER for M=16');
        legend('Location', 'southwest')
    end
         
end

fprintf('Part A: ');
toc

%% part b
clc;clear;close all
tic;

numIterations = 10000;  % The number of iterations of the simulation
numSym = 1000;
n = 30;
%30 represents 2 msg (14bits) and its parity bits (16bits) being used as 
%training bits
%n specifics how many pairs of msg and parity strings will be used as
%training bits
%so numTraining is total number of training bits being used (bits => binary)
numTraining = 30*n;

SNR_Vec = 0:2:16;
SNRlen = length(SNR_Vec);

chan = [1, 0.2, 0.4];
tic;
M = 4;
codeWordLen = 15;
msgLen = 7;
%286 codewords transmitted which is 286*7=2002 bits or 1001 symbols
numWords = ceil(numSym*log2(M)/msgLen);

%number of msg bits that have to be used for training 840
trainingBits = log2(M) * (numTraining/codeWordLen) * msgLen; 

BERvec2 = zeros(numIterations, SNRlen);

enc = comm.BCHEncoder(codeWordLen, msgLen);
dec = comm.BCHDecoder(codeWordLen, msgLen);


for ii=1:numIterations
    %total bits we can generate is numWords*msgLen=2002 bits
    msg = randi([0, 1], msgLen*numWords, 1);
    numBits = size(msg,1);
    %encode with parity bits inserted after every 7 bits aka every msg
    msg_enc = step(enc, msg);
    a = msg_enc;
    
    %pair up the bits to convert to decimal if necessary (M>2)
    msg_enc = reshape(msg_enc,log2(M),length(msg_enc)/log2(M)).';
    b = msg_enc;
    
    %convert to decimal
    msg_enc = bi2de(msg_enc);
    
    parfor jj=1:SNRlen
        %default gray code ordering and integer inputs 
        tx = qammod(msg_enc, M, 'UnitAveragePower', true);

            txChan = filter(chan,1,tx);  % Apply the channel.
            txNoisy = awgn(txChan,SNR_Vec(jj)); % add noise
            
            %make the eq
            %Some previous attempts
                %lineq = comm.LinearEqualizer('Algorithm','LMS', 'NumTaps',6,'StepSize',0.01); % doesn't work on matlab 2018a
                %eq1 = lineareq(6, lms(0.001));
                %txNoisy = filter(eq1.weights, 1, txNoisy);
            eq1 = dfe(12,6, lms(0.01)); 
            eq1.SigConst = qammod(0:M-1, M, 'UnitAveragePower', true);
            eq1.ResetBeforeFiltering = 0;
            
            txNoisy = equalize(eq1,txNoisy,tx(1:numTraining));
             
            reset(eq1); % clean up - to be removed?
        
        rx = qamdemod(txNoisy, M, 'UnitAveragePower', true);
        rxTmp = (de2bi(rx)).'; %transpose here  % 'left-msb'
        rxMsg = rxTmp(:);
        
        dec_msg = step(dec, rxMsg);
        
        [~, BERvec2(ii,jj)] = biterr(msg(trainingBits+1:end), dec_msg(trainingBits+1:end));  
    end
end

ber2 = mean(BERvec2,1);

% grab last non zero value
ber2(find(ber2,1,'last'));

bit_rate = (numBits - trainingBits)/numBits

figure;
semilogy(SNR_Vec, ber2, 'DisplayName', "BER-4 with ISI")
hold on;
berTheory2 = berawgn(SNR_Vec,'qam', 4,'nondiff');
semilogy(SNR_Vec,berTheory2,'DisplayName', 'Theoretical BER for M=4')
legend('Location', 'southwest')
fprintf('Part B: ');
toc
berTheory2
