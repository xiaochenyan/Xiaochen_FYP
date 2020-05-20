u=rand(1,2000000)<0.5;    % create a random input of 2000000 bits
Ec=ones(1,70);   % Assign  symbol energy Ec to be 1.
Es=Ec;
Eb1=ones(1,40);  % Bit energy should be 7/4 of Ec that is 7/4.
Eb=(7/4)*Eb1;
N0=[0:0.1:4-0.1];    % create 40 different PSD of noise
sigma=(N0/2).^0.5;
d=Es.^0.5;    % figure out the value of d based on Es
lenu=length(u);
lenc=lenu*(7/4);
BER=zeros(1,40);
pber=zeros(1,40);

G=[1 1 0 1 0 0 0;0 1 1 0 1 0 0;1 1 1 0 0 1 0;1 0 1 0 0 0 1];    % The generator mastrix
H=[1 0 0 1 0 1 1;0 1 0 1 1 1 0;0 0 1 0 1 1 1];  % Find the parity check matrix using G
c2=zeros(1,7);
c1=zeros(1,3500000);
c=zeros(1,3500000);
rdemod=zeros(1,3500000);
s2=zeros(1,3);
s1=zeros(1,1500000);
s=zeros(1,1500000);
error2=zeros(1,7);
error=zeros(1,3500000);
udecod2=zeros(1,4);
udecod=zeros(1,2000000);
r=zeros(1,3500000);

% Use message bits and generator matrix to create the codewords
for i=1:1:lenu/4
    c2=[u(4*i-3) u(4*i-2) u(4*i-1) u(4*i)]*G;
    for j=1:1:7
        c1(7*(i-1)+j)=c2(j);
    end
end
c=mod(c1,2);    

% codewords to vectors mapping 
vec=zeros(1,lenc);
for i=1:1:lenc  
    if c(i)==0
%        if input 000 symbol='s1'
        vec(i)=-d(1);
    elseif c(i)==1
%        if input 100 symbol='s2'
        vec(i)=d(1);
    end 
end




% Measure 40 samples
for p=1:1:40
    %%AWGN
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
 

    % Error detection(find the syndrone)
    for i=1:1:lenc/7
        s2=[rdemod(7*i-6) rdemod(7*i-5) rdemod(7*i-4) rdemod(7*i-3) rdemod(7*i-2) rdemod(7*i-1) rdemod(7*i)]*H';
        for j=1:1:3
            s1(3*(i-1)+j)=s2(j);
        end
    end
    s=mod(s1,2);

    % Find the error vector using syndrone
    for m=1:1:length(s)/3
        if [s(3*m-2) s(3*m-1) s(3*m)]==[0 0 0]
            error2=[0 0 0 0 0 0 0];
            for n=1:1:7
                error(7*(m-1)+n)=error2(n);
            end
        elseif [s(3*m-2) s(3*m-1) s(3*m)]==[0 0 1]
            error2=[0 0 1 0 0 0 0];
            for n=1:1:7
                error(7*(m-1)+n)=error2(n);
            end
        elseif [s(3*m-2) s(3*m-1) s(3*m)]==[0 1 0]
            error2=[0 1 0 0 0 0 0];
            for n=1:1:7
                error(7*(m-1)+n)=error2(n);
            end
        elseif [s(3*m-2) s(3*m-1) s(3*m)]==[0 1 1]
            error2=[0 0 0 0 1 0 0];
            for n=1:1:7
                error(7*(m-1)+n)=error2(n);
            end
        elseif [s(3*m-2) s(3*m-1) s(3*m)]==[1 0 0]
            error2=[1 0 0 0 0 0 0];
            for n=1:1:7
                error(7*(m-1)+n)=error2(n);
            end
        elseif [s(3*m-2) s(3*m-1) s(3*m)]==[1 0 1]
            error2=[0 0 0 0 0 0 1];
            for n=1:1:7
                error(7*(m-1)+n)=error2(n);
            end
        elseif [s(3*m-2) s(3*m-1) s(3*m)]==[1 1 0]
            error2=[0 0 0 1 0 0 0];
            for n=1:1:7
                error(7*(m-1)+n)=error2(n);
            end
        elseif [s(3*m-2) s(3*m-1) s(3*m)]==[1 1 1]
            error2=[0 0 0 0 0 1 0];
            for n=1:1:7
                error(7*(m-1)+n)=error2(n);
            end
        end
    end

    % Find the original codewords using mod-2 addition of received vectors and
    % error vectors. 
    cdecod1=rdemod+error;
    cdecod=mod(cdecod1,2);

    % find the message bits stored in codewords(systematic coding)
    for i=1:1:length(cdecod)/7
        udecod2=[cdecod(7*i-3) cdecod(7*i-2) cdecod(7*i-1) cdecod(7*i)];
        for j=1:1:4
            udecod(4*(i-1)+j)=udecod2(j);
        end
    end

    % figure out BER
    sumbit=0;
    for count=1:1:lenu
        if u(count)~=udecod(count)
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
