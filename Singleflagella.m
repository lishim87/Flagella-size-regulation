
clear all
N = 200; %total number of motors
Nb = N; %motors at base
Nf = N-Nb; %number of the filament
k_inj = 1; %injection rate
k_d = 0.03; %disassembly
k_f = 2; % choosing filament
b=0;

%Attributing motors (index, location, position, activity)
%index = 1 to 200
%location = base(0) or filament (1)
%position = position on the filament in range (1,L)
%activity = diffusing (0) and active transport (1)
motor = zeros(N ,4);

%initialize the motor attributes
L = 2; %length of flagella
for i=1:N
    motor(i) = i;
    if i<=Nb
        motor(i,2)= 0;
    else
        motor(i,2)= 1;
        motor(i,3)= randi([1 L]);
        motor(i,4)=1;
    end
    
end

motor;

Steps= 2000000; %total number of steps
T= zeros(1,Steps); % time (mins) vector
T(1)=0; %initial time
%tau=zeros(1,Steps); %time steps
m1 = zeros(1, Steps);
m1(1)=2;
Ntot= 1000;
monomers = Ntot;
Thr = 1; %Threshold for injection

for i= 1:Steps
    %i
    k1 = k_inj;
    if Nb<Thr
        k1=0;
    end
    k2 = k_d;
    if (m1(i)<=1)
        k2=0;
    end
    
    k3 = Nf*k_f;
    if Nf<=0
        k3=0;
    end
    %k1=0;
    k0 = k1+k2+k3;
    
   
    
    C1=rand;
    tau(i)=(1/k0)*log(1/C1);
    T(i+1)= T(i)+tau(i);
    
    C2= rand;
    
    if C2<=k1/k0  %injection
        %disp('inject')
        
        h = find(motor(:,2)== 0); %number in base
        %size(h,1)= number at the base
        for j= 1 : Thr
        pos0 = randi(length(h));
        ch = h(pos0);
        motor(ch,2:4)= 1;
        h = h(h~=ch);
        end
        Nb= Nb-Thr;
        Nf= Nf+Thr;
        m1(i+1)=m1(i);
        
    else
        if (C2 > k1/k0) && (C2 <= (k1+k2)/k0) % decay
            %disp('subtract')
            m1(i+1)= m1(i)- 1;
            monomers = monomers + 1;
            
        else      %choose motor activity
            
            h = find(motor(:,2)==1);
            pos = randi(length(h));
            ch1 = h(pos);
            %disp('translate')
            
            if  (motor(ch1,4)==1)
                %   disp('translate+')
                motor(ch1,3)= motor(ch1,3)+ 1; %increment the position
                if (motor(ch1,3)>= m1(i))
                    m1(i+1)= m1(i)+ 1; %increase the length
                    motor(ch1,4)=0; %unbind 
                    motor(ch1,3)= m1(i);
                    monomers=monomers-1; %decrese free monomers
                else
                    
                    m1(i+1)= m1(i);
                end
                
            else %diffuse
                m1(i+1)=m1(i);
                %disp('here')
                p1= rand;
                if p1< 0.5
                    %disp('here+')
                    if (motor(ch1,3)< m1(i))
                       motor(ch1,3)= motor(ch1,3)+2; 
                    else
                           
                        %motor(ch1,3)= m1(i)- 2;
                        motor(ch1,3)= motor(ch1,3)- 2;
                    end
                    
                    
                else
                    %disp('here-')
                    if motor(ch1,3)<= m1(1)
                        motor(ch1,2:4)=0;
                        Nf= Nf-1;
                        Nb=Nb+1;
                       % disp('base')
                        b=b+1;
                    else
                        motor(ch1,3)= motor(ch1,3)-2; 
                                          
                    end
                end
            end
        end
        
        
    end
end

plot(T,m1)

