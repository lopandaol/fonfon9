close all; clear all;
Fs=8000; bits = 16 ;
lengthpad = 0.7*Fs;                       % % t=0.7s
%%%%%%%%%%%%%% record
% disp('Dont speak')
% recnoise = audiorecorder(Fs,bits,1);    % 1 channel
% recordblocking(recnoise, 2);            % 2 seconds
% noise = getaudiodata(recnoise);         % noisefloor
% 
% recvoice = audiorecorder(Fs,bits,1);    % 1 channel
% disp('Start speaking.')
% recordblocking(recvoice, 2);            % 2 seconds
% voice = getaudiodata(recvoice);         % noisefloor

noise = audioread('noisetesttwomos1.wav');
voice = audioread('testtwomos1.wav');
%%%%%%%%%%%%%% cut
A = 3*max(noise);
L = length(voice);                              
xcut = zeros(L,1);
for i=1:L
   if (abs(voice(i))>=A)
   	    xcut(i) = voice(i);
   end    
end
PF = find(xcut,1,'first');
PE = find(xcut,1,'last');

y =  voice(PF:PE);                          % cut voice
if length(y)<= lengthpad
        y = y;
else y = y(1:lengthpad);
end
    
%%%%%%%%%%%%%% pre-emphasis
b = [1  -0.95]; a = 1; 
preem = filter(b,a,y);                      % Pre-emphasised signal

%%%%%%%%%%%%%% pad
pad = zeros(lengthpad,1);                    	
lpreem = length(preem);
for i=1:lpreem;
    if (i<=lpreem)
        pad(i)=preem(i);
    else pad(i)=pad(i);                   % y:zero-paded signal
    end
end

% %%%%%%%%%%%%%% framing & hamming window
L = 255+1;   
w = hamming(L);
numfram = round((lengthpad-L)/99)-1 ;
k = zeros(numfram+1,L);
% ham = zeros(numfram+1,L);
k(1,:)= pad(1:256);
ham = k(1,:).*w' ; 

for N=1:numfram                            % Number_of_frames=54; (t=0.7s) %     
    start = 100+(N-1)*99;
    stop = start+255;

    k(N+1,:) = pad(start:stop);
    ham(N+1,:) = k(N+1,:).*w' ;       
end

%%%%%%%%%%%%%% DFT
   
for i=1:numfram+1      
        x = ham(i,:);
        m = length(x);                                    % Window length
        n = pow2(nextpow2(m));                            % Transform length
%         n = 256;
        y = fft(x,n);                                       % dft
        f = (0:n/2-1)*(Fs/n);                               % Frequency range
        sptr(i,:) = (abs(y(1:128)));  
end

%%%%%%%%%%%%%% Mel-Frequency
nfft = 256;
lower = 300;
upper = Fs/2;
hfs = upper;
mN = 28;

melf = zeros(1,mN);
convf = zeros(1,mN);
fi = zeros(1,mN);
filbk = zeros(mN-2,nfft/2);

mf_upper = 1125*log(1+ (upper/700));
mf_lower = 1125*log(1+ (lower/700));

mf_gap = (mf_upper-mf_lower)/(mN - 1);

melf(1) = mf_lower;
convf(1) = lower;

fi(1)=floor(((nfft+1)*convf(1))/Fs);

fVal = Fs*[0:nfft/2-1]/nfft;

for i=2:mN
    melf(i) = melf(i-1) + mf_gap;
    convf(i) = 700*(exp(melf(i)/1125)-1);
    fi(i)=floor(((nfft+1)*convf(i))/Fs);
end

nFb = length(melf);
H = zeros(nFb-2,nfft/2);

for k=2:nFb-1
   for j=1:nfft/2
       if (j < fi(k-1)),
           H(k,j) = 0;
       elseif ((j >= fi(k-1)) && (j <= fi(k))),
           H(k,j) = ( j - fi(k-1) ) / (fi(k) - fi(k-1));
       elseif (( j >= fi(k) ) && ( j <= fi(k+1) )),
           H(k,j) = ( fi(k+1) - j )/( fi(k+1) - fi(k) );
       elseif ( j > fi(k+1) ),
           H(k,j) = 0;
       end
   end
end

for i=2:27
   filbk(i-1,:) = H(i,:);
end

fllfrm = zeros(54,128);

%%%%%%%%%%%%%% filter bank
for jj = 1:54
        for kk = 1:26
            fltbked = sptr(jj,:).*filbk(kk,:);
            sumspect = sum(fltbked);
            cep(jj,kk) = (sumspect)^2;
        end
end

ceplog = log10(cep);

for ii=1:54
        aa = sum(cep(ii,:));
        if (aa == 0),
            cep(ii,:) = 1;
        end
end

ceplog = log10(cep);

%%%%%%%%%%%%%% cepstral coefficient
    for jj = 1:54
        mcep(jj,:) = dct(ceplog(jj,:),26);
        CT_KW(jj,:) = mcep(jj,(2:13));         % บันทึกค่าสปส.ที่ได้สำหรับหนึ่งไฟล์เสียง
    end
    
%%%%%%%%%%%%%% DELTA CEPSTRAL COEFFICIENT
y = CT_KW;
prep_del = zeros(1,12);
prep_del = [prep_del; y;];
bN = size(prep_del);

N = bN(1);
DT_KW= zeros(bN(1)-1,bN(2));

for i=2:N-1

    DT_KW(i-1,:) = (prep_del(i+1,:) - prep_del(i-1,:) ) /2;

end

%%%%%%%%%%%%%% energy
for i=1:54
            
    hh = ham(i,:);
    %energy(i) = sum(hh.^2);
    energy(i) = (hh(1).^2) + (hh(2).^2);
    energy = energy';
    
end

%%%%%%%%%%%%%% delta energy
prep_energy = zeros(1);
prep_energy = [prep_energy; energy;];
bN = size(prep_energy);

N = bN(1);
DT_energy = zeros(bN(1)-1,bN(2));

for i=2:N-1

 DT_energy(i-1,:) = (prep_energy(i+1,:) - prep_energy(i-1,:) ) /2;

end

%%%%%%%%%%%%%% 26 cepstrum
aa = CT_KW; bb=energy; cc=DT_KW; dd=DT_energy;
coeff26= [aa, bb, cc, dd];

%%%%%%%%%%%%%% merge MFCC
for m=1:54
    if (m == 1);
        stt = coeff26(m,:);
    else yy = coeff26(m,:);
        stt =  [stt yy];
    end
end


