clear all;
close all;
clc 
fs=8000;

%initial condition========================================
n=0.3;             % time index
st=1; 
name = 'fon'

for rr=1:2
%     figure
    filename = strcat('cutforsagc_testtwo',name, num2str(rr), '.wav');
    [mt1,Fs]=audioread(filename);

    mt1=mt1(1:n*8000);
    mt1=mt1';
    mt1=mt1-mean(mt1);
    mt1= mt1/(max(abs(mt1)));
    
    lp_cutoff = 500;
    [bl,al] = butter(4, lp_cutoff/(Fs/2), 'low');
    mt = filter(bl, al, mt1);                                     %%% mt is prepared signal
    siglength=length(mt);
    t=linspace(0,n,n*8000);
    [pssig,sq,fdsig,sfsig,signal1]=phaseshift(mt,siglength,st,t);
%     lp_cutoff = 2000;
%     [B,A] = butter(4, lp_cutoff/(Fs/2), 'low');
%     output=filter(B,A,real(pssig));                             %%% output is output of SAGC
    output = real(pssig);

%     %%%%%%% fft 
   
        m = length(output);                                    % Window length
        nn = pow2(nextpow2(m));                            % Transform length
        nn = 1024;
        y = fft(output,nn);                                       % dft
        f = (0:nn/2-1)*(fs/nn);                               % Frequency range
        sptr(rr,:) = (abs(y(1:nn/2)));
        
        maxsagc(rr)= max(sptr(rr,:));
%         
%         [r,c] = find(sptr(rr,:)>100);
%         for w = 1:length(c)
%             po = c(w);
%             fre(rr,w) = f(po);
%         end
%         w=0; c=0;
%      figure
%        subplot(4,1,1)
%      plot(mt1)
% %      % plot(t,mt1)
%      grid on
% % %     axis([0 n min(mt1) max(mt1)])
%      title('SPEECH')
% ylabel('Amplitude');xlabel('Sample')
%      subplot(4,1,2)
%      plot(mt);ylabel('Amplitude');xlabel('Sample')
% % %     % plot(t,mt)
%      grid on
% % %     axis([0 n min(mt) max(mt)])
%      title('PREPERATION PROCESS')
%      subplot(4,1,3)
%      plot(output);ylabel('Amplitude');xlabel('Sample')
% % %      plot(t,output)
%      grid on
%      ylabel('Amplitude')
%      ylim([-1-0.5 1+0.5])
%      title('OUTPUT OF SAGC SCHEME');ylabel('Amplitude');xlabel('Sample')
% %     
%      subplot(4,1,4)
%      stem(f,sptr(rr,:))
%      title('FFT of SAGC');ylabel('Amplitude');xlabel('Frequency')
% % % 
% % 
% 
% 
% 




% %         
% %         %%%%% LPF to envelope 
         lp_cutoff = 200;
         [B,A] = butter(4, lp_cutoff/(Fs/2), 'low');
         enve=filter(B,A,real(sq));                             %%% output is output of SAGC
         nn = 1024;
         y = fft(enve,nn);                                       % dft
         f = (0:nn/2-1)*(fs/nn);                               % Frequency range
         sptr_enve(rr,:) = (abs(y(1:nn/2)));
         sptr_enve(rr,1)=0;
%          sptr_enve(rr,2)=0;
%          sptr_enve(rr,3)=0;
%          sptr_enve(rr,4)=0;
%          sptr_enve(rr,5)=0;
%          sptr_enve(rr,6)=0;
%          xxx = zeros(1,length(sptr_enve));
         
         maxenve(rr)= max(sptr_enve(rr,:));
         
% %         
% %         
% 
figure
        subplot(4,1,1);grid on
        plot(mt1); title('Speech'); hold on
        ylabel('Amplitude');xlabel('Sample')
%         plot(enve,'r')
        subplot(4,1,2)
     plot(real(sq));grid on
     title('Pitch Detection');hold on
     ylabel('Amplitude');xlabel('Sample')
%      plot(enve,'r')
     subplot(4,1,3)
     plot(enve);grid on
     title('Envelope of Pitch')
     ylabel('Amplitude');xlabel('Sample')
    subplot(4,1,4)
     stem(f,sptr_enve(rr,:))
     title('FFT of Envelope')
     ylabel('Amplitude');xlabel('Frequency')

% 
% % stem(f,sptr_enve(rr,:))



end



bb = maxenve';
aa = maxsagc';
merge = [aa bb];
% eval(strcat('newSAGCsptr_',name,'=merge;'))
% save(strcat('newSAGCsptr_',name),strcat('newSAGCsptr_',name))

% eval(strcat('SAGCsptr_',name,'=sptr_enve;'))
% save(strcat('SAGCsptr_',name),strcat('SAGCsptr_',name))

% mix = [sptr sptr_enve];
% eval(strcat('newmaxSAGCandenvo_',name,'=merge;'))
% save(strcat('newmaxSAGCandenvo_',name),strcat('newmaxSAGCandenvo_',name))

% 
% eval(strcat('freq_',name,'=fre;'))
% % save(strcat('freq_',name),strcat('freq_',name))
% 
% [r1,c1,v1] = find(fre);
% freqmo = mode(v1);
% for jj=1:rr
%     if jj==1
%         pj = fre(jj,:);
%     else
%         pj = [ pj fre(jj,:)];
%     end
% end
% ma = mode(pj);