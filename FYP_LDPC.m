u=rand(1,1500000)<0.5;    % create a random input of 2000000 bits
Ec=ones(1,70);   % Assign  symbol energy Ec to be 1.
Es=Ec;
Eb1=ones(1,30);  % Bit energy should be 7/3 of Ec that is 7/3.
Eb=(7/3)*Eb1;
N0=[0:0.1:3-0.1];    % create 30 different PSD of noise
sigma=(N0/2).^0.5;
d=Es.^0.5;    % figure out the value of d based on Es
lenu=length(u);
lenc=lenu*(7/3);
BER=zeros(1,30);

G=[1 1 0 0 1 0 0;0 1 1 0 0 1 0;1 0 1 1 0 0 1];    % The generator matrix
H=[1 0 0 0 1 0 1;0 1 0 0 1 1 0;0 0 1 0 0 1 1;0 0 0 1 0 0 1];  % Find the parity check matrix using G
c2=zeros(1,7);
c1=zeros(1,3500000);
c=zeros(1,3500000);
r=zeros(1,lenc);
rdemod=zeros(1,3500000);
rdec=zeros(1,3500000);
udec=zeros(1,1500000);

% Use message bits and generator matrix to create the codewords
for i=1:1:lenu/3
    c2=[u(3*i-2) u(3*i-1) u(3*i)]*G;
    for j=1:1:7
        c1(7*(i-1)+j)=c2(j);
    end
end
c=mod(c1,2);    

% codewords to vectors mapping 
vec=zeros(1,lenc);
for i=1:1:lenc  
    if c(i)==0
        vec(i)=-d(1);
    elseif c(i)==1
        vec(i)=d(1);
    end 
end




% Measure 30 samples
for p=1:1:30
    %%Z-channel
    noise=zeros(1,lenc);    % create a noise vector
    for k=1:1:lenc
        noise(k)=sigma(p)*rand(1);  % create a 1-dimension noise vector with variance we want for each vector
        if noise(k)<=0.65*sigma(p)
            noise(k)=0;
        elseif noise(k)>0.65*sigma(p)
            noise(k)=-d(1)*sigma(p)-1;
        end
        
    end
    
    for i=1:1:lenc
        if vec(i)==-d(1)
            r(i)=-d(1);
        elseif vec(i)==d(1)
            r(i)=vec(i)+noise(i);
        end  % add signal vector to each noise vector(only ones can be flipped)
    end
% demodulator -- for each received vector, decide which symbol it
% should belong to. In this case, I use the decision rule.
    for l=1:1:lenc
        if r(l)<=0
            rdemod(l)=0;
        elseif r(l)>0
            rdemod(l)=1;
        end   
    end
    I=0;
    Imax=20;
    M=zeros(1,7);
    % Error detection
    for i=1:1:lenc/7
        for j=1:1:7
            M(j)=rdemod(7*(i-1)+j);
        end
        
         L1=mod((M(1)+M(5)+M(7)),2);
         L2=mod((M(2)+M(5)+M(6)),2);
         L3=mod((M(3)+M(6)+M(7)),2);
         L4=mod((M(4)+M(7)),2);
         
        for k=1:1:20
            if (I<Imax)&((L1~=0)|(L2~=0)|L3~=0|L4~=0)
                E11_1=M(5)+M(7);
                E11=mod(E11_1,2);
                E15_1=M(1)+M(7);
                E15=mod(E15_1,2);
                E17_1=M(1)+M(5);
                E17=mod(E17_1,2);
                E22_1=M(5)+M(6);
                E22=mod(E22_1,2);
                E25_1=M(2)+M(6);
                E25=mod(E25_1,2);
                E26_1=M(2)+M(5);
                E26=mod(E26_1,2);
                E33_1=M(6)+M(7);
                E33=mod(E33_1,2);
                E36_1=M(3)+M(7);
                E36=mod(E36_1,2);
                E37_1=M(3)+M(6);
                E37=mod(E37_1,2);
                E44=M(7);
                E47=M(4);
                
                if M(1)~=E11
                    M(1)=mod(M(1)+1,2);
                end
                
                if M(2)~=E22
                    M(2)=mod(M(2)+1,2);
                end
                
                if M(3)~=E33
                    M(3)=mod(M(3)+1,2);
                end
                
                if M(4)~=E44
                    M(4)=mod(M(4)+1,2);
                end
                
                if (M(5)~=E15)&(M(5)~=E25)
                    M(5)=mod(M(5)+1,2);
                end
                
                if (M(6)~=E26)&(M(6)~=E36)
                    M(6)=mod(M(6)+1,2);
                end
                
                if (M(7)~=E17)&(M(7)~=E37)&(M(7)~=E47)
                    M(7)=mod(M(7)+1,2);
                end
                
                L1=mod((M(1)+M(5)+M(7)),2);
                L2=mod((M(2)+M(5)+M(6)),2);
                L3=mod((M(3)+M(6)+M(7)),2);
                L4=mod((M(4)+M(7)),2);
                I=I+1;
            end
        end
        for l=1:1:7
            rdec(7*(i-1)+l)=M(l);
        end
    end



    % figure out BER
    for i=1:1:lenc/7
        u2=[rdec(7*i-2) rdec(7*i-1) rdec(7*i)];
        for j=1:1:3
            udec(3*(i-1)+j)=u2(j);
        end
    end
    
    sumbit=0;
    for count=1:1:lenu
        if u(count)~=udec(count)
            sumbit=sumbit+1;    % count the amount of different bits
        end
    end
    BER(p)=sumbit/lenu;
    % used for x-axis
    EbN0(p)=Eb(p)/N0(p);
end
% Plot the BER curve
figure
semilogy(10*log10(EbN0),BER,'r');
xlabel('Eb/N0(dB)')
ylabel('BER')