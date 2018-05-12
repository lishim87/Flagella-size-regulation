close all
clear all
N = 20; %total number of motors
Nb = N; %motors at base
Nf = N-Nb; %number of the filament
k_inj = 1; %injection rate
k_d = 0; %disassembly
k_f = 2; % choosing filament
b=0;

%Attributing motors (index, location, position, activity, filament id)
%index = 1 to 200
%location = base(0) or filament (1)
%position = position on the filament in range (1,L)
%activity = diffusing (0) and active transport (1)
%filament id = 1 or 2
motor = zeros(N ,5);

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
        motor(i,5)= randi(2);
    end
    
end

%motor

Steps= 50; %total number of steps
T= zeros(1,Steps); % time (mins) vector
T(1)=0; %initial time
%tau=zeros(1,Steps); %time steps
m1= zeros(1, Steps);
m1(1)=2;
m2= zeros(1, Steps);
m2(1)=2;

Ntot= 1000;
monomers = Ntot;
Thr= 1; %Threshold for injection

for i= 1:Steps
    i
    m=zeros(1, Steps);
    k1 = Nb*k_inj; %inject in to 1/2
    if Nb<Thr
        k1=0;
    end
    
    k2= k_d;
    if (m1(i)<=m1(1))
        k2=0;
    end
    
    k3 = k_d;
    if (m2(i)<=m2(1))
        k3=0;
    end
    
    k4 = Nf*k_f;
    if monomers==0
        k4=0;
    end
    
    k0 = k1+k2+k3+k4;
    
    
    
    C1=rand;
    tau(i)=(1/k0)*log(1/C1);
    T(i+1)= T(i)+tau(i);
    
    C2= rand;
    
    if C2<=k1/k0  %injection
        %disp('inject')
        
        h = find(motor(:,2)==0); %number in base
        %size(h,1)= number at the base
        for j=1:Thr
            j;
            pos0 = randi(length(h));
            ch = h(pos0);
            motor(ch,2:4)= 1;
            motor(ch, 5)= randi(2);
            h = h(h~=ch);
        end
        %  motor
        Nb= Nb-Thr;
        %Nf= Nf+Thr;
        Nf1=size(motor(motor(:,5)==1),1);
        Nf2=size(motor(motor(:,5)==2),1);
        Nf= Nf1+Nf2;
        
        m1(i+1)=m1(i);
        m2(i+1)=m2(i);
        
    else
        if (C2 > (k1/k0)) && (C2 <= (k1+k2)/k0) % decay
            disp('subtract1')
            m1(i+1)= m1(i)- 1;
            m2(i+1)= m2(i);
            monomers = monomers + 1;
            b = b-1;
            
        else
            
            if (C2 > (k1+k2)/k0) && (C2 <= (k1+k2+k3)/k0) % decay
                disp('subtract2')
                m2(i+1)= m2(i)- 1;
                m1(i+1)=m1(i);
                monomers = monomers + 1;
                b=b-1;
                
            else 
                if (C2 > (k1+k2+k3)/k0) && (C2 <= (k1+k2+k3+k4)/k0)%choose motor activity
                    h = find(motor(:,2)==1);
                    pos = randi(length(h));
                    ch1 = h(pos);
                    %disp('translate')
                    %motor(ch1,5)
                    if motor(ch1,5)==1
                        m = m1;
                    else
                        if motor(ch1,5)==2
                        m = m2;
                        end
                    end
                    disp('choose filament')
                    motor(ch1,5)
                    m
                    if  (motor(ch1,4)== 1)  %   disp('translate+')
                        motor(ch1,3)= motor(ch1,3)+ 1; %increment the position
                        
                        if (motor(ch1,3)>= m(i))
                            if monomers > 0
                                disp('+length')
                                m(i+1)= m(i)+ 1; %increase the length
                                monomers = monomers - 1; %decrease free monomers
                                b=b+1;
                                motor(ch1,4)=0; %unbind
                                motor(ch1,3)= m(i);
                                if motor(ch1,3)<0
                                    disp('bac')
                                end
                            else
                                disp('no monomers')
                                break;
                            end
                        else
                            disp('+length not stay')
                            m(i+1)= m(i);
                        end
                        
                    else %diffuse
                        m(i+1)=m(i);
                        %disp('here')
                        p1= rand;
                        if p1< 0.5 %forward
                            %disp('here+')
                            if (motor(ch1,3)<= m(i))
                                motor(ch1,3)= motor(ch1,3)+ 1;
                            else
                                motor(ch1,3)= motor(ch1,3)- 1;
                                if motor(ch1,3)<0
                                    disp('for')
                                end
                                
                            end
                        else %backward
                            %disp('here-')
                            if motor(ch1,3)<= (m(1)+1)
                                motor(ch1,2:5)=0;
                                Nf= Nf-1;
                                Nb= Nb+1;
                                % disp('base')
                                %b=b+1;
                            else
                                motor(ch1,3)= motor(ch1,3)-1;
                                if motor(ch1,3)<0
                                    disp('bac')
                                end
                                
                            end
                        end
                        
                    end
                    
                    if motor(ch1,5)==1
                        m1(i+1)= m(i+1);
                        m2(i+1)= m2(i);
                        disp('m1')
                        m1
                        m2
                        
                    else
                        m2(i+1)= m(i+1);
                        m1(i+1)=m1(i);
                        disp('m2')
                        m1
                        m2
                    end
                    

                else
                    disp('boo')
                end
                
            end
        end
    end
   disp('diff')
   m1(i+1)- m1(1)+(m2(i+1)-m2(1))- b 
end

%m1(end)-m1(1)+(m2(end)-m2(1))

plot(T(1:Steps),m1(1:Steps)-m1(1),T(1:Steps),m2(1:Steps)-m2(1))

