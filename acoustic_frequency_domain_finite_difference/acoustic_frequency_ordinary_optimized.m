%implementation of wavefield model of acoustic media with finite-difference  frequency domain approach 

clear all;
% grid parameter
thick = 20;  %PML thickness
dx = 5;
dz = 5;
m = 200+2*thick;
n = 200+2*thick;
mo = m-2*thick; % size of the real model 
% m = 5;
% n = 5;
Dx = dx*dx;
Dz = dz*dz;

dt = 1/1024;  %time sampling

% grid parameters of sparse matrix
M = m*n;
N = m*n;
% material parameter
Vp = 2e3;

den = 2.1e3;
b = 1/den;
K= Vp^2*den*ones(m,n);
% weighting coefficients
a1 = 0.5461;
a2 =1 - a1;
c1 = 0.6248;
c2 = 0.0938;
c3 = (1 - c1 -4*c2)/4;

%source
f=25;
w= 2*pi*f;
W= 30;
F = zeros(N,1);
max_f = 80;  %maximum frequency

%frequency field
P = zeros(mo,mo);

%PML
A = 100;  % dampling factor
complex = sqrt(-1);
ex = ones(2*m-1,2*n-1);
ez = ones(2*m-1,2*n-1);
for i=1:thick*2+1
    ex(:,[i,2*n-i])= 1-complex/w*A*cos(pi/2*(i-1)/(2*thick));
    ez([i,2*m-i],:)= 1-complex/w*A*cos(pi/2*(i-1)/(2*thick));
end

%sparse matrix
S=sparse(M,N,0);

% time max 1s, sampling frequency 1024Hz
max_t = 1024;
% save the frequency field of the points used for ifft
% lines are the frequencies of a single point, rows are different points
% fft and ifft requires time length and frequency length to be the equal
% fft and ifft in Matlab require the frequency series to be
% symmetrical, so before ifft, we need to add the symmetrical part of the
% frequency series
Q=zeros(mo^2,max_t);

for f= 1:1:max_f
    %source
    w= 2*pi*f;
    F(m/2*n+m/2,1) = 2/sqrt(pi)*(f^2)/(W^2)*exp(-(f/W)^2);
    disp('%%%%%%%')
    disp(f);
    disp('%%%%%%%')
    
    for x=1:M
        if x==1 %Line 1
            t = 1;
            S(x,1) = c1*w^2/K(t,1) - a1*1/Dx*b/ex(1+2*(t-1),1+2*(1-1))/ex(1+2*(t-1),1+2*(1-1)+1) -...
                a1*1/Dz*b/ez(1+2*(t-1),1+2*(1-1))/ez(1+2*(t-1)+1,1+2*(1-1));
            S(x,2) = c2*w^2/K(t,1) + a1*1/Dx*b/ex(1+2*(t-1),1+2*(1-1))/ex(1+2*(t-1),1+2*(1-1)+1);
            
            S(x,1+n) = c2*w^2/K(t,1) + a1*1/Dz*b/ez(1+2*(t-1),1+2*(1-1))/ez(1+2*(t-1)+1,1+2*(1-1));
        end
        
        if (x>=2)&&(x<=(n-1))  %from Line 2 to Line n-1
            i = x-1; 
            t = 1;
            S(x,i+0) = c2*w^2/K(t,i+1) + a1*1/Dx*b/ex(1+2*(t-1),1+2*(i+1-1))/ex(1+2*(t-1),1+2*(i+1-1)-1);
            S(x,i+1) = c1*w^2/K(t,i+1) - a1*1/Dx*b/ex(1+2*(t-1),1+2*(i+1-1))/ex(1+2*(t-1),1+2*(i+1-1)+1) -...
                a1*1/Dx*b/ex(1+2*(t-1),1+2*(i+1-1))/ex(1+2*(t-1),1+2*(i+1-1)-1) -...
                a1*1/Dz*b/ez(1+2*(t-1),1+2*(i+1-1))/ez(1+2*(t-1)+1,1+2*(i+1-1));
            S(x,i+2) = c2*w^2/K(t,i+1) + a1*1/Dx*b/ex(1+2*(t-1),1+2*(i+1-1))/ex(1+2*(t-1),1+2*(i+1-1)+1);
            
            
            S(x,i+1+n) = c2*w^2/K(t,i+1) + a1*1/Dz*b/ez(1+2*(t-1),1+2*(i+1-1))/ez(1+2*(t-1)+1,1+2*(i+1-1));
        end
        
        if x==n     %Line n
            i = x-1;
            t = 1;
            S(x,i+0) = c2*w^2/K(t,i+1) + a1*1/Dx*b/ex(1+2*(t-1),1+2*(i+1-1))/ex(1+2*(t-1),1+2*(i+1-1)-1);
            S(x,i+1) = c1*w^2/K(t,i+1) - a1*1/Dx*b/ex(1+2*(t-1),1+2*(i+1-1))/ex(1+2*(t-1),1+2*(i+1-1)-1) -...
                a1*1/Dz*b/ez(1+2*(t-1),1+2*(i+1-1))/ez(1+2*(t-1)+1,1+2*(i+1-1));
            
            
            S(x,i+1+n) = c2*w^2/K(t,i+1) + a1*1/Dz*b/ez(1+2*(t-1),1+2*(i+1-1))/ez(1+2*(t-1)+1,1+2*(i+1-1));
            disp(1);
        end
        
        if (x>=n+1)&&(x<=(m-1)*n)   %from Line n+1 to Line (m-1)*n
            i=ceil(x/n)-1;    
            if (x-i*n)==1 
                t = i+1;
                S(x,1+(i-1)*n) = c2*w^2/K(t,1) + a1*1/Dz*b/ez(1+2*(t-1),1+2*(1-1))/ez(1+2*(t-1)-1,1+2*(1-1));
                
                S(x,1+(i+0)*n) = c1*w^2/K(t,1) - a1*1/Dx*b/ex(1+2*(t-1),1+2*(1-1))/ex(1+2*(t-1),1+2*(1-1)+1) -...
                    a1*1/Dz*b/ez(1+2*(t-1),1+2*(1-1))/ez(1+2*(t-1)+1,1+2*(1-1)) - ...
                    a1*1/Dz*b/ez(1+2*(t-1),1+2*(1-1))/ez(1+2*(t-1)-1,1+2*(1-1));
                S(x,2+(i+0)*n) = c2*w^2/K(t,1) + a1*1/Dx*b/ex(1+2*(t-1),1+2*(1-1))/ex(1+2*(t-1),1+2*(1-1)+1);
                
                S(x,1+(i+1)*n) = c2*w^2/K(t,1) + a1*1/Dz*b/ez(1+2*(t-1),1+2*(1-1))/ez(1+2*(t-1)+1,1+2*(1-1));
            end
            
            if ((x-i*n)>=2)&&((x-i*n)<=(n-1)) 
                j = x-i*n-1;
                t = i+1;
                S(x,j+1+(i-1)*n) = c2*w^2/K(t,j+1) + a1*1/Dz*b/ez(1+2*(t-1),1+2*(j+1-1))/ez(1+2*(t-1)-1,1+2*(j+1-1));
                
                S(x,j+0+(i+0)*n) = c2*w^2/K(t,j+1) + a1*1/Dx*b/ex(1+2*(t-1),1+2*(j+1-1))/ex(1+2*(t-1),1+2*(j+1-1)-1);
                S(x,j+1+(i+0)*n) = c1*w^2/K(t,j+1) - a1*1/Dx*b/ex(1+2*(t-1),1+2*(j+1-1))/ex(1+2*(t-1),1+2*(j+1-1)+1) -...
                    a1*1/Dx*b/ex(1+2*(t-1),1+2*(j+1-1))/ex(1+2*(t-1),1+2*(j+1-1)-1) -...
                    a1*1/Dz*b/ez(1+2*(t-1),1+2*(j+1-1))/ez(1+2*(t-1)+1,1+2*(j+1-1)) - ...
                    a1*1/Dz*b/ez(1+2*(t-1),1+2*(j+1-1))/ez(1+2*(t-1)-1,1+2*(j+1-1));
                S(x,j+2+(i+0)*n) = c2*w^2/K(t,j+1) + a1*1/Dx*b/ex(1+2*(t-1),1+2*(j+1-1))/ex(1+2*(t-1),1+2*(j+1-1)+1);
                
                S(x,j+1+(i+1)*n) = c2*w^2/K(t,j+1) + a1*1/Dz*b/ez(1+2*(t-1),1+2*(j+1-1))/ez(1+2*(t-1)+1,1+2*(j+1-1));
            end
            
            if (x-i*n)==n 
                j = x-i*n-1;
                t = i+1;
                S(x,j+1+(i-1)*n) = c2*w^2/K(t,j+1) + a1*1/Dz*b/ez(1+2*(t-1),1+2*(j+1-1))/ez(1+2*(t-1)-1,1+2*(j+1-1));
                
                S(x,j+0+(i+0)*n) = c2*w^2/K(t,j+1) + a1*1/Dx*b/ex(1+2*(t-1),1+2*(j+1-1))/ex(1+2*(t-1),1+2*(j+1-1)-1);
                S(x,j+1+(i+0)*n) = c1*w^2/K(t,j+1) - a1*1/Dx*b/ex(1+2*(t-1),1+2*(j+1-1))/ex(1+2*(t-1),1+2*(j+1-1)-1) -...
                    a1*1/Dz*b/ez(1+2*(t-1),1+2*(j+1-1))/ez(1+2*(t-1)+1,1+2*(j+1-1)) -...
                    a1*1/Dz*b/ez(1+2*(t-1),1+2*(j+1-1))/ez(1+2*(t-1)-1,1+2*(j+1-1));
                
                S(x,j+1+(i+1)*n) = c2*w^2/K(t,j+1) + a1*1/Dz*b/ez(1+2*(t-1),1+2*(j+1-1))/ez(1+2*(t-1)+1,1+2*(j+1-1));
                
                disp(i+1);
            end
        end
        
        if x==((m-1)*n+1)
            t = m;
            S(x,1+(m-2)*n) = c2*w^2/K(t,1) + a1*1/Dz*b/ez(1+2*(t-1),1+2*(1-1))/ez(1+2*(t-1)-1,1+2*(1-1));
            
            S(x,1+(m-1)*n) = c1*w^2/K(t,1) - a1*1/Dx*b/ex(1+2*(t-1),1+2*(1-1))/ex(1+2*(t-1),1+2*(1-1)+1) -...
                a1*1/Dz*b/ez(1+2*(t-1),1+2*(1-1))/ez(1+2*(t-1)-1,1+2*(1-1));
            S(x,2+(m-1)*n) =  c2*w^2/K(t,1) + a1*1/Dx*b/ex(1+2*(t-1),1+2*(1-1))/ex(1+2*(t-1),1+2*(1-1)+1);
        end
        
        if (x>=((m-1)*n+2))&&(x<=((m-1)*n+n-1))    
            i = x-(m-1)*n-1; 
            t = m;
            S(x,i+1+(m-2)*n) = c2*w^2/K(t,i+1) + a1*1/Dz*b/ez(1+2*(t-1),1+2*(i+1-1))/ez(1+2*(t-1)-1,1+2*(i+1-1));
            
            S(x,i+0+(m-1)*n) = c2*w^2/K(t,i+1) + a1*1/Dx*b/ex(1+2*(t-1),1+2*(i+1-1))/ex(1+2*(t-1),1+2*(i+1-1)-1);
            S(x,i+1+(m-1)*n) = c1*w^2/K(t,i+1) - a1*1/Dx*b/ex(1+2*(t-1),1+2*(i+1-1))/ex(1+2*(t-1),1+2*(i+1-1)+1) -...
                a1*1/Dx*b/ex(1+2*(t-1),1+2*(i+1-1))/ex(1+2*(t-1),1+2*(i+1-1)-1) -...
                a1*1/Dz*b/ez(1+2*(t-1),1+2*(i+1-1))/ez(1+2*(t-1)-1,1+2*(i+1-1));
            S(x,i+2+(m-1)*n) = c2*w^2/K(t,i+1) + a1*1/Dx*b/ex(1+2*(t-1),1+2*(i+1-1))/ex(1+2*(t-1),1+2*(i+1-1)+1);
        end
        
        if x==((m-1)*n+n)    
            i = x-(m-1)*n-1; 
            t = m;
            S(x,i+1+(m-2)*n) = c2*w^2/K(t,i+1) + a1*1/Dz*b/ez(1+2*(t-1),1+2*(i+1-1))/ez(1+2*(t-1)-1,1+2*(i+1-1));
            
            S(x,i+0+(m-1)*n) = c2*w^2/K(t,i+1) + a1*1/Dx*b/ex(1+2*(t-1),1+2*(i+1-1))/ex(1+2*(t-1),1+2*(i+1-1)-1);
            S(x,i+1+(m-1)*n) = c1*w^2/K(t,i+1) - a1*1/Dx*b/ex(1+2*(t-1),1+2*(i+1-1))/ex(1+2*(t-1),1+2*(i+1-1)-1) -...
                a1*1/Dz*b/ez(1+2*(t-1),1+2*(i+1-1))/ez(1+2*(t-1)-1,1+2*(i+1-1));
            
            disp(m);
        end
    end
    
    % %solve the linear system by LU decomposition
    [L,U] = lu(S);
    y = U\(L\F);
    
    for i = 1:mo
        P(i,:) = y(n*(thick+i-1)+thick+1:n*(thick+i-1)+n-thick,1);
    end
    
    % rearrange the frequency field 
    % line is frequency, row is the point
    for i = 1:mo
        for j = 1:mo
             Q((i-1)*mo+j,f+1) = P(i,j); %add the zero frequency
        end
    end
    
    
end

T1=zeros(mo,mo);  %t = 0.1s;
T2=zeros(mo,mo);  %t = 0.2s;
T3=zeros(mo,mo);  %t = 0.4s;

% add the symmetrical part of the frequency used for ifft
for i = 1:max_f
    Q(:,max_t-(i-1)) = conj(Q(:,i+1));
end

freq = zeros(1,max_t);
time = zeros(1,max_t);

% ifft for every point
for i = 1:(mo*mo)
    freq(1,:)= Q(i,:);     %frequency field of a single point
    time(1,:) = ifft(freq);   %time field of a single point
      
    T1(ceil(i/mo),i - (ceil(i/mo)-1)*mo) = time(1,round(0.1*max_t));
    T2(ceil(i/mo),i - (ceil(i/mo)-1)*mo) = time(1,round(0.2*max_t));
    T3(ceil(i/mo),i - (ceil(i/mo)-1)*mo) = time(1,round(0.4*max_t));
end

frequency_slice = reshape(Q(:,30),[mo,mo,]);

figure,hold on
imagesc(real(frequency_slice));figure(gcf);
hold off
axis equal
axis tight

figure,hold on
imagesc(real(T1));figure(gcf);
hold off
axis equal
axis tight

figure,hold on
imagesc(real(T2));figure(gcf);
hold off
axis equal
axis tight

figure,hold on
imagesc(real(T3));figure(gcf);
hold off
axis equal
axis tight


























