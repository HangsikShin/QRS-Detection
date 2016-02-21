%%
% [loc, time] = QRSdetection(wdata,Fs,wtime,ftype)
%    
% This code is for QRS complex detection of Electrocardiogram
%
% loc and time is indices and times of QRS location, respectively.   
% wdata is ECG waveform data, Fs is the sampling frequency, wtime is the
% time information of the wdata and ftype is the filter type (MIT-BIH or
% AHA)
%
% This source code is for the following article : Jinkwon Kim and Hangsik
% Shin, "Simple and Robust Realtime QRS Detection Algorithm based on
% Spatiotemporal Characteristic of the QRS Complex", Plos One, 2016
%
% article link : TBD
%
% Jinkwon Kim and Hangsik Shin (hangsik.shin@gmail.com)
% @ Healthcare Solution Laboratory, Chonnam National University
%
% 2015.12.31      
%
function [loc, time] = QRSdetection(wdata,Fs,wtime,ftype)

if(ftype==0)  % MIT-BIH
    load('filterCoeff');
else % AHA
   load('filterCoeffAHA');
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Parameters
%
%   wsize1 : MAF size for Energy Level Detection
%   wsize2 : MAF size for Energy Variation Detection
%   refractory_time : Refrqctory Period
%   thEL0 : threshold for Energy Level Detection
%   stabLevel : Stabilization Reference Voltage
%   r_a : application rate for weight adjustment
%   r_b : application rate for weight adjustment 
%   r_nr : application rate of noise level (for signal threshold)
%   r_s : application rate of signal level (for noise threshold)
%   r_d : decay rate for adaptive threshold
%   r_n : application rate of noise level (for signal threshold)
%   Weight : weight for adjustment signal level
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wsize1 = 0.15;  % size of 1st MAF
wsize2 = 0.2;   % size of 2nd MAF
refractory_time = 0.15;   
thEL0 = 0.1;  % threshold for Energy Level Detection
stabLevel = 0.5;
r_a = 0.1;
r_b = 0.05;
r_nr = 1.75;
r_s = 0.001;
r_d = 0.05;
r_n = 0.03;
Weight = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Window Configuration
%
% 	winsizeEL : window size for Energy Level (EL)
% 	winsizeEV : window size for Energy Variation (EV)
% 	thEL : Threshold for EL (Adaptive threshold)
% 	thEL0 : Initial value for EL threshold
%   refrctoryP : refractory period
%   thEV : Threshold for EV (Hard threshold)
%   thEVmax : Maximum thEV
%   thEVmin : Minimum thEV
%   cDecay : decay constant for thEL
%   nk1 : constant for adjustment of noise threshold (mult)
%   nk2 : constant for adjustment of noise threshold (add)
%   nk3 : constant for adjustment of noise threshold (inc)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

winsizeEL = round(wsize1*Fs);
winsizeEV = round(wsize2*Fs);
diffWinsize = winsizeEV - winsizeEL;
refractoryP = round(refractory_time*Fs);    

thEVlimit = 1*Fs/(0.2*Fs*20);
thEVub = 0.45*Fs/(0.2*Fs*20);
thEVlb = -0.45*Fs/(0.2*Fs*20);

thEVub2 = 20*Fs/(0.2*Fs*20);
thEVlb2 = -20*Fs/(0.2*Fs*20);

ArrayL = 5; 
decayF = 1-1/( (0.40-refractory_time)*Fs ); 

checker2=0; 
loc=[];
isStart=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Preprocessing
%
%   fSig : Signal after BPF
%   sSig : Signal after squaring
%   dSig : Signal after differentiated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fSig=filtfilt(Num,1,wdata);
sSig=sqrt(fSig.^2);
dSig=[0; diff(sSig)]*Fs;
sigLen = length(sSig);

%% Memory Allcation
ELQRS = zeros(sigLen,1);
EVQRS = zeros(sigLen,1);
thEL = ones(sigLen,1) * thEL0;
thEV = zeros(sigLen-1,1);
thN = zeros(sigLen,1);
maxVArry=zeros(ArrayL,1);
maxDifBuf=zeros(ArrayL,1);
minDifBuf=zeros(ArrayL,1);
BUF1=zeros(winsizeEL,1);
BUF2=zeros(winsizeEV,1);

LargeWin = winsizeEV;

maxV=1;
QRScount = 0;
    
for kk=LargeWin:sigLen-1
    % MOVING AVERAGE with WEIGHT
    BUF1(kk,1)=sum( sSig(kk-winsizeEL+1:kk,1) )/ winsizeEL;
    BUF2(kk,1)=sum( dSig(kk-winsizeEV+1:kk,1) )/ winsizeEV;
    
    ELQRS(kk,1) = sum( BUF1(kk-winsizeEL+1:kk,1) )/winsizeEL;  
    EVQRS(kk,1) = sum( BUF2(kk-winsizeEV+1:kk,1) )/winsizeEV;
    
    % STEP 1: ENERGY LEVEL DETECTION
    if (isStart==0) && (ELQRS(kk) >= thEL(kk))
        thEL(kk)=ELQRS(kk);
        maxV=ELQRS(kk);
        maxP=kk;
        isStart = 1;
    end
    if(ELQRS(kk)<thN(kk))
        thN(kk)=ELQRS(kk);
    end
    if(isStart)
        if (ELQRS(kk)>=maxV)
            thEL(kk)=ELQRS(kk);
            maxV=ELQRS(kk);
            maxP=kk;
            Timer=refractoryP;
        else
            Timer=Timer-1;
            thEL(kk)=maxV;
            if(Timer==0)
                isStart=0;
                checker2=1;
                TimerOfPeak=winsizeEV-(refractoryP-winsizeEL);
                maxP_Buf=maxP;
                maxV_Buf=maxV;
            end
        end
    end
    
    % STEP 2: ENERGY VARIATION DETECTION  
    if(checker2==1)
        TimerOfPeak=TimerOfPeak-1;
        if(TimerOfPeak==0)    
            checker2=0;
            if maxP_Buf-winsizeEL<1
                BufStartP2=1;
            else
                BufStartP2=maxP_Buf-winsizeEL;
            end
            if maxP_Buf+2*diffWinsize>sigLen
                BufEndP2=size(M,1);
            else
                BufEndP2=maxP_Buf+2*diffWinsize*2;
            end

            DiffSumCheck1=max(EVQRS(BufStartP2:maxP_Buf+diffWinsize,1));
            DiffSumCheck2=min(EVQRS(maxP_Buf+diffWinsize:BufEndP2,1));
            
            if ( isempty(loc) || (DiffSumCheck1-DiffSumCheck2>thEVlimit)&&(DiffSumCheck1*DiffSumCheck2<0)&&(DiffSumCheck1>thEVub)&&(DiffSumCheck2<thEVlb)&&(DiffSumCheck1<thEVub2)&&(DiffSumCheck2>thEVlb2) )
                QRScount = QRScount + 1;
                                
                loc=[loc; maxP_Buf-winsizeEL+2];        
                % STEP 3: WEIGHT ADJUSTMENT
                maxVArry(mod(QRScount,ArrayL)+1,1)=maxV_Buf;
                maxDifBuf(mod(QRScount,ArrayL)+1,1)=max(EVQRS(BufStartP2:BufEndP2,1));
                minDifBuf(mod(QRScount,ArrayL)+1,1)=min(EVQRS(BufStartP2:BufEndP2,1));
                if stabLevel>mean(maxVArry) 
                    AdujR1=min((stabLevel-median(maxVArry))*r_a, stabLevel*r_b);
                else
                    AdujR1=min((stabLevel-median(maxVArry))*r_a, stabLevel*-r_b);
                end
                Weight=Weight+AdujR1;
            end
        end
    end
    thN(kk+1)=thN(kk)+r_s*ELQRS(kk);
    if(exist('maxV_Buf','var'))
        thEL(kk+1)=thEL(kk)* ( decayF  * (1-r_d*(thN(kk)/maxV_Buf))) + r_n*thN(kk);
    else
        thEL(kk+1)=thEL(kk)* ( decayF);
    end
    if thEL(kk+1)<r_nr*thN(kk)
        thEL(kk+1)=r_nr*thN(kk);
    end
end
time=wtime(loc);
end

% End of the Code