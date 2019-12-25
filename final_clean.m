% Matt Kribs, Richard Lee, Lucy Xu
% Comm Theory Final Project
% Description: A BER script for a wireless link simulation

%% part a
clear all;clc;close all

%This section has two parts to it. The first part is matching AWGN
%performance for M-ary QAM schemes for M = 4 and 16. The second part to
%this section is to obtain a BER of 10e-4 using BPSK (M=2) after the message
%bits have been modulated, channeled, added to awgn, equalized, and demodulated.

numIter = 10000;              % number of iterations of the simulation
nSym = 1000;                % number of symbols per packet
SNR_Vec = 0:2:16;           % SNR value for the awgn function
lenSNR = length(SNR_Vec);   
m_ary = [2, 4, 16];         % M-ary number, 2 corresponds to binary modulation
%chan = 1;                  % No channel
chan = [1, 0.2, 0.4];       % moderate ISI

% Not so invertible, severe ISI, we did not use this channel for any part
% of this project
%chan = [0.227 0.460 0.688 0.460 0.227]';   

tic;
trainlen = 300;             % number of training bits 
displayStr = ["BER-2 with ISI","BER-4 No ISI", "BER-16 No ISI"]; % plot labels

for it=1:length(m_ary)
    M = m_ary(it);
    
    berVec = zeros(numIter, lenSNR);    % preallocate BER matrix
    
    parfor ii = 1:numIter
        % we generate a random decimal message everytime
        msg = randi([0, M-1], nSym*(log2(M)), 1);  
 
        % we turn the msg into bits here
        % the transpose is so we can columnize
        % we use bits to for the biterr func, but do not transmit this 
        bits = de2bi(msg, 'left-msb').'; 
        bits = bits(:);
        
        for jj = 1:lenSNR        % one iteration of the simulation at each SNR Value
            tx = qammod(msg,M);  % modulate the msg
               
            % for 2-ary we filter, add awgn, and equalize
            if M == 2
                
                % if there is no channel then the received msg is the same
                % as the transmitted msg
                if isequal(chan,1) 
                    txChan = tx;
                    txNoisy = txChan;
                else
                    txChan = filter(chan,1,tx);         % apply the channel
                    txNoisy = awgn(txChan,SNR_Vec(jj)); % add AWGN
                               
                    % Originally we tried the linear equalizer but could
                    % not get it down enough, so we switched to dfe
                   
                    % lineq = comm.LinearEqualizer('Algorithm','LMS', 'NumTaps',6,'StepSize',0.01); 
                    % doesn't work on matlab 2018a
                    
                    % we also started off with simply feedforward taps but
                    % it did not work well enough so we used feedback taps
                    % as well
                    % eq1 = lineareq(6, lms(0.01));
                    
                    eq1 = dfe(12,6, lms(0.01)); 
                    eq1.SigConst = qammod(0:M-1, M, 'UnitAveragePower', true);
                    eq1.ResetBeforeFiltering = 1;
            
                    txNoisy = equalize(eq1,txNoisy,tx(1:trainlen)); % equalize
                end
                
            % if we are 4-ary or 16-ary we only add awgn
            else 
                txNoisy = awgn(tx, SNR_Vec(jj) + 10*log10(log2(M)),'measured');
            end
            
            rx = qamdemod(txNoisy,M); % demodulate the signal
            
            % Convert the received message to bits
            % the transpose and the rx(:) is housekeeping - conceptually we are
            % taking each row, appending it after the previous row, but we do
            % this transposed since we are working with columns
            rxTmp = de2bi(rx, 'left-msb').'; 
            rxBits = rxTmp(:);
            
            % Compute and store the BER for this iteration
            % We're interested in the BER, which is the 2nd output of BITERR
            [~, berVec(ii,jj)] = biterr(bits(trainlen+1:end), rxBits(trainlen+1:end));  
            
        end  % End SNR iteration
    end      % End numIter iteration
    
    % Compute and plot the mean BER
    ber = mean(berVec,1);
    
    % print the BER at SNR = 12 for 2-ary
    if M == 2
        BER_at_12 = ber(7)
    end
    
    % here we plot the figures
    figure(it);
    semilogy(SNR_Vec, ber, 'DisplayName', displayStr(it))
    hold on;
    
    % plot different theoreticals depending on M
    if M == 2
        berTheory2 = berawgn(SNR_Vec,'psk', 2,'nondiff');
        semilogy(SNR_Vec,berTheory2,'DisplayName', 'Theoretical BER for M=2')
    elseif M == 4
        berTheory4 = berawgn(SNR_Vec,'qam', 4,'nondiff');
        semilogy(SNR_Vec,berTheory4,'DisplayName', 'Theoretical BER for M=4')
    elseif M == 16
        berTheory16 = berawgn(SNR_Vec,'qam', 16,'nondiff');
        semilogy(SNR_Vec,berTheory16, 'DisplayName', 'Theoretical BER for M=16');
    end
    
    legend('Location', 'southwest')
    title('Part A');
    xlabel('SNR in dB');
    ylabel('Bit Error Rate');
    
end
fprintf('Part A: ');
toc

%% part b

% In this part, an encoding scheme will be added to the system to improve
% the BER to 10e-6

%% BPSK with 1/2 convolutional encoding
% clear all;clc;close all

% we initially tried 1/2 convolutional encoding since its less tedious to
% implement than block codes but this proved to not be enough so please
% find our trials with BCH

% numIter = 0;
% nSym = 500;    
% SNR_Vec = 0:2:16;
% lenSNR = length(SNR_Vec);
% trainlen = 300;      
% M_values = 2;
% chan = [1, 0.2, 0.4];
% 
% displayStr = ["BER-2 with ISI","BER-4 No ISI", "BER-16 No ISI"];
% 
% tic;
% for it=1:length(M_values)
%     m_ary = M_values(it);
%     berVec = zeros(numIter, lenSNR);
%     for ii = 1:numIter
% 
%         msg = randi([0, m_ary-1], nSym*(log2(m_ary)), 1); 
%         numBits = size(msg,1);
%         % New bits must be generated at every iteration
% 
%         %encoding 1/2 convolutional
%         K = 3;                %constraint length is how long a bit can affect the encoder output
%         g_1 = '1+x+x^2';                
%         g_2 = '1+x^2';        
%         trellis = poly2trellis(K,{g_1,g_2});    %trellis is created 
%         msg_enc = convenc(msg,trellis);         %input bits are encoded
%         
%         bits = de2bi(msg_enc, 'left-msb').'; %transpose here 
%         bits = bits(:);
%         
%         for jj = 1:lenSNR               % one iteration of the simulation at each SNR Value
%             tx = qammod(msg_enc,m_ary); % BPSK modulate the signal
%             
%             m_ary == 2
%                 txChan = filter(chan,1,tx);         % apply the channel
%                 txNoisy = awgn(txChan,SNR_Vec(jj)); % add AWGN
%                     
%                 sigConst = qammod(0:m_ary-1, m_ary, 'UnitAveragePower', true);
%                 eq1 = dfe(12,6,lms(0.01)); 
%                 eq1.SigConst = sigConst; 
%                 eq1.ResetBeforeFiltering = 0;
%                    
%                 txNoisy = equalize(eq1,txNoisy,tx(1:trainlen)); % equalize
%                 reset(eq1);
%                     
%                 rx = qamdemod(txNoisy,m_ary); %,'OutputType', 'integer'); % Demodulate
%                 rxTmp = de2bi(rx, 'left-msb').'; 
%                 
%             tblen = (K - 1)*5;   %positive integer scalar that specifies the traceback depth
%                                  %for rate 1/2, a typical value for tblen is about
%                                  %five times the constraint length - 1
%                                       
%             rxMSG = (vitdec(rxTmp,trellis,tblen,'trunc','hard')).';
%             %trunc: encoder is assumed to have started at the all zero state
%             %hard: code contains binary input values                         
%             
%             % Compute and store the BER for this iteration
%             % We're interested in the BER, which is the 2nd output of BITERR
% 
%             [~, berVec(ii,jj)] = biterr(msg(trainlen+1:end), rxMSG(trainlen+1:end)); 
% 
%         end  % End SNR iteration
%     end      % End numIter iteration
% end

% % Compute and plot the mean BER
%ber = mean(berVec,1);
    
% grab last non zero value and report on screen
%lowest_non_zero_BER = ber(find(ber,1,'last'))

% % calculate bit rate
% bit_rate = (numBits - trainlen)/numBits

% % plot the figure
% figure('Name', 'Part B');
% semilogy(SNR_Vec, ber, 'DisplayName', "BER-2 with ISI")
% hold on;
% berTheory2 = berawgn(SNR_Vec,'psk', 2,'nondiff');
% semilogy(SNR_Vec,berTheory2,'DisplayName', 'Theoretical BER for M=2')
% legend('Location', 'southwest')
% title('Reducing Error to 10e-6 for BPSK');
% xlabel('SNR in dB');
% ylabel('Bit Error Rate');
% fprintf('Part B: ');
% toc

%% BPSK with BCH 15-7 encoding
% % clc;clear;close all
% 
% % We thought that we had to stick to BPSK (2-ary) so we have this section
% % where we used BCH 15-7 encoding on 2-ary but after we realized we could
% % have used QAM, we tried that out in the next section
% 
% numIterations = 1;  
% numSymbols = 1000;
% 
% % We are allowed roughly 1000 symbols: we use this to find roughly how many words
% % we can generate for BCH 15-7
% % we use a ceiling function below to calculate precisely how many codewords
% % we can send
%  
% % fter encoding, the number of total bits (including parity) that we use for training
% numTraining = 150;
% 
% SNR_Vec = 0:2:16; 
% SNRlen = length(SNR_Vec);
% 
% % same channel, Moderate ISI
% chan = [1, 0.2, 0.4];
% 
% tic;
% 
% M = 2;
% codeWordLen = 15;
% msgLen = 7;
% 
% % the ceiling is to round up the number of words
% 
% % later, we use num words * msgLen to figure out how many bits we can
% % generate, knowing that encoding will add 8 parity bits to each symbol
% numWords = ceil(numSymbols/codeWordLen);
% 
% % number of training bits that we had to take from the original message
% trainingBits = (numTraining/codeWordLen) * msgLen;
% 
% % make a 0 vector
% BERvec2 = zeros(numIterations, SNRlen);
% 
% % We use the comm BCH encoder and decoder objects, make them once and reuse
% enc = comm.BCHEncoder(codeWordLen, msgLen);
% dec = comm.BCHDecoder(codeWordLen, msgLen);
% 
% for ii=1:numIterations
%     % make a msg that is number of msg bits long, such that after encoding, 
%     % there are approximately 100- symbols transmitted 
%     msg = randi([0, M-1], msgLen * numWords, 1);
%     numBits = size(msg,1);
%     
%     % BCH encode it. 469 bits generated above, 1005 symbols transmitted
%     msg_enc = step(enc, msg);
%     
%     for jj=1:SNRlen
%         
%         tx = qammod(msg_enc, M);            % modulate the signal
%         txChan = filter(chan,1,tx);         % apply the channel
%         txNoisy = awgn(txChan,SNR_Vec(jj)); % add AWGN
%         
%         eq1 = dfe(12,6, lms(0.01));         % create equalizer object
%         eq1.SigConst = qammod(0:M-1, M, 'UnitAveragePower', true);
%         eq1.ResetBeforeFiltering = 1;
%        
%         txNoisy = equalize(eq1,txNoisy,tx(1:numTraining));  % equalize
%         rx = qamdemod(txNoisy, M);                          % demodulate
%         dec_msg = step(dec, rx);                            % decode
%         
%         % calculate BER
%         [~, BERvec2(ii,jj)] = biterr(msg(trainingBits+1:end), dec_msg(trainingBits+1:end));  
%     end
% end
% 
% % average BER for each SNR
% ber2 = mean(BERvec2,1);
% 
% % grab last non zero value and report on screen
% lowest_non_zero_BER = ber2(find(ber2,1,'last'))
% 
% % calculate bit rate
% bit_rate = (numBits - trainingBits)/numBits
% 
% % plot the figure
% figure('Name', 'Part B');
% semilogy(SNR_Vec, ber2, 'DisplayName', "BER-2 with ISI")
% hold on;
% berTheory2 = berawgn(SNR_Vec,'psk', 2,'nondiff');
% semilogy(SNR_Vec,berTheory2,'DisplayName', 'Theoretical BER for M=2')
% legend('Location', 'southwest')
% title('Reducing Error to 10e-6 for BPSK');
% xlabel('SNR in dB');
% ylabel('Bit Error Rate');
% fprintf('Part B: ');
% toc

%% QAM with BCH 15-7
% clc;clear;close all

% QAM with 4-ary should be able to imporve the bit rate because the number 
% of bits per symbol is doubled 

% same parameters from part a
numIterations = 100;
numSym = 1000;
SNR_Vec = 0:2:16;
SNRlen = length(SNR_Vec);
chan = [1, 0.2, 0.4];

% number of codewords*2 that is used for training
n = 60;    
% number of symbols that is used for training aka number of msg bits*2 
% that is used for training
numTraining = 15*n;     

M = 4;
codeWordLen = 15;
msgLen = 7;

% 286 codewords transmitted for the given parameters which is 286*7 = 2002 
% bits aka 1001 symbols
numWords = ceil(numSym*log2(M)/msgLen);

% number of bits from 'msg' that have to be used for training
trainingBits = log2(M)*(numTraining/codeWordLen) * msgLen; 

tic;
% preallocate BER matrix like previously
BERvec2 = zeros(numIterations, SNRlen);     

% create encoder and decoder objects for the given BCH parameters
enc = comm.BCHEncoder(codeWordLen, msgLen); 
dec = comm.BCHDecoder(codeWordLen, msgLen);

parfor ii=1:numIterations
    
    % total bits we can generate is numWords*msgLen = 2002 bits = 1001 symbols
    msg = randi([0, 1], msgLen*numWords, 1);
    numBits = size(msg,1);
    
    % encode msg with parity bits inserted after every 7 bits aka every
    % message
    msg_enc = step(enc, msg);
    a = msg_enc;                % checkpoint
    
    % pair up the bits to convert to decimal
    msg_enc = reshape(msg_enc,log2(M),length(msg_enc)/log2(M)).';
    b = msg_enc;                % checkpoint
    
    % convert to decimal
    msg_enc = bi2de(msg_enc);
    
    for jj=1:SNRlen
        % qammod default parameters is gray code ordering and integer inputs 
        tx = qammod(msg_enc, M, 'UnitAveragePower', true);

            txChan = filter(chan,1,tx);         % apply the channel
            txNoisy = awgn(txChan,SNR_Vec(jj),'measured'); % add awgn

            eq1 = dfe(12,6, lms(0.01));         % create equalizer objecct
            eq1.SigConst = qammod(0:M-1, M, 'UnitAveragePower', true);
            eq1.ResetBeforeFiltering = 0;
            
            txNoisy = equalize(eq1,txNoisy,tx(1:numTraining));  % equalize
            reset(eq1);
        
        rx = qamdemod(txNoisy, M, 'UnitAveragePower', true);    % demodulate
        
        % convert back to binary so that recieved msg can be decoded
        rxTmp = (de2bi(rx)).';
        rxMsg = rxTmp(:);           % columnize
        
        dec_msg = step(dec, rxMsg); % decode
        
        % calculate BER
        [~, BERvec2(ii,jj)] = biterr(msg(trainingBits+1:end), dec_msg(trainingBits+1:end));  
    end
end

ber2 = mean(BERvec2,1);

% grab last non zero value of BER vector 
lowest_non_zero_BER = ber2(find(ber2,1,'last'))

bit_rate = (numBits - trainingBits)/numBits

figure;
semilogy(SNR_Vec, ber2, 'DisplayName', "BER-4 with ISI")
title('Part B')
hold on;
berTheory2 = berawgn(SNR_Vec,'qam', 4,'nondiff');
semilogy(SNR_Vec,berTheory2,'DisplayName', 'Theoretical BER for M=4')
legend('Location', 'southwest')
fprintf('Part B: ');
toc
