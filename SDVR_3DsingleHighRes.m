warning('off','all')
clear all

%% constants
% potential of length l and potential height 100

for a=2:5
    tic

        lx=3.7e-9;%str2double(uinput(1))*10^-9;     %convert string from user into a number

        ly=3.7e-9;%str2double(uinput(2))*10^-9;     %convert string from user into a number

        lz=26.1e-9;%str2double(uinput(3))*10^-9;     %convert string from user into a number

        nx=[0.5 .4 .3 .2 .1 .08 .07 .06 .05 .04 .06 .05 .04 .03 .02 .01].*10^(-9); %str2double(uinput(4))*10^-9;     %convert string from user into a number

        ny=[0.5 .4 .3 .2 .1 .08 .07 .06 .05 .04 .06 .05 .04 .03 .02 .01].*10^(-9); %str2double(uinput(5))*10^-9;     %convert string from user into a number        

        nz=[0.8 .8 .8 .8 .8 1.0 1.0 1.0 1.0 1.0 1.0 .05 .04 .03 .02 .01].*10^(-9); %str2double(uinput(6))*10^-9;     %convert string from user into a number
        
        shape=3; %str2double(uinput(7));        %convert string from user into a number %% 1=ellipsoid 2=rectangular prisim 3=rod

p=10000;
m_h=9.10938356e-31*.45; %mass of hole
m_e=9.10938356e-31*.13; %mass of electron
mass=m_h*m_e/(m_h+m_e); %mass 
%9.10938356e-31*.45; %mass of hole
%9.10938356e-31*.13; %mass of electron
h_bar=1.0545718e-34;

%% define length of x array

x=-lx/2-nx(a):nx(a):lx/2+nx(a);      %did -nx(a) and +nx(a) to go 1 over the shape
y=-ly/2-ny(a):ny(a):ly/2+ny(a);
z=-lz/2-nz(a):nz(a):lz/2+ny(a);
x2=1:length(x).^3;
z2=x2;
%% set up particle in a box potential for plotting purposes only

lenx=length(x);
leny=length(y);
lenz=length(z);
len3=lenx*leny*lenz;
lenlots=len3^2;
len_max=max(max(lenx,leny),lenz);

% k and k_prime are indexes of the H matrix
% xi and yj are the components make up k
% xi_prime and yj_prime are the components of k_prime
count=0;
xi=0;
yj=0;
zj=0;
k=0;
k_prime=0;
T1=0;
T2=0;
T3=0;
V=zeros(len3,1);
aa=zeros(floor(lenlots/5),1);
bb=zeros(floor(lenlots/5),1);
cc=zeros(floor(lenlots/5),1);
H1=0;

switch shape
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
case 1 %ellipsoid
        hbnx1=h_bar*h_bar/(6*mass*nx(a)*nx(a))*pi*pi;
        hbnx2=h_bar*h_bar/(mass*nx(a)*nx(a));
        hbny1=h_bar*h_bar/(6*mass*ny(a)*ny(a))*pi*pi;
        hbny2=h_bar*h_bar/(mass*ny(a)*ny(a));
        hbnz1=h_bar*h_bar/(6*mass*nz(a)*nz(a))*pi*pi;
        hbnz2=h_bar*h_bar/(mass*nz(a)*nz(a));
        
        xi=         zeros(len3,1);
        yj=         zeros(len3,1);
        zj=         zeros(len3,1);
        
for k=1:len3        
      
        xi(k)=         floor((k-1)/(lenz*leny))+1;
        yj(k)=         mod(floor((k-1)/lenz),leny)+1;% mod(k-1,lenx)+1;
        zj(k)=         mod((k-1),lenz)+1;  
    
      if ((x(xi(k))*x(xi(k)))/(lx*lx/4)+(y(yj(k))*y(yj(k)))/(ly*ly/4)+(z(zj(k))*z(zj(k)))/(lz*lz/4))>1 
          V(k)=p;
      end
end
    


    for k=1:len3
    
        %parametrical equations for xi and yj
%         holder= zeros(1,len3);
        
    for k_prime=1:len3
        %parametrical equations for xi_prime and yj_prime
 
   %Kinetic Energy T1 and V
      if yj(k)==yj(k_prime) && zj(k)==zj(k_prime)  
            if xi(k)==xi(k_prime)
                T1=hbnx1+V(k);
            else
                T1=(-1)^(xi(k)-xi(k_prime))*hbnx2/((xi(k)-xi(k_prime))*(xi(k)-xi(k_prime)));
            end
      else
            T1=0; %kronecker delta
        
      end
        %kinetic energy T2
        
      if xi(k)==xi(k_prime) && zj(k)==zj(k_prime)  
            if yj(k)==yj(k_prime)
                   T2=hbny1;
            else
                T2=(-1)^(yj(k)-yj(k_prime))*hbny2/((yj(k)-yj(k_prime))*(yj(k)-yj(k_prime)));
            end
     else
            T2=0; %kronecker delta
     end  
     
     %kinetic energy T3
      
      if xi(k)==xi(k_prime) && yj(k)==yj(k_prime)  
          if zj(k)==zj(k_prime)
                   T3=hbnz1;
          else
               T3=(-1)^(zj(k)-zj(k_prime))*hbnz2/((zj(k)-zj(k_prime))*(zj(k)-zj(k_prime)));
          end
        
      else
            T3=0; 
      end      
      

        %H = T + V   
        H1 = T1 + T2 + T3;
         if H1 ~= 0
            count=count+1;
            aa(count)=k;
            bb(count)=k_prime;
            cc(count)=H1;
         end
      
    end

     
 end    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  case 2 %rectangular Prism

        hbnx1=h_bar*h_bar/(6*mass*nx(a)*nx(a))*pi*pi;
        hbnx2=h_bar*h_bar/(mass*nx(a)*nx(a));
        hbny1=h_bar*h_bar/(6*mass*ny(a)*ny(a))*pi*pi;
        hbny2=h_bar*h_bar/(mass*ny(a)*ny(a));
        hbnz1=h_bar*h_bar/(6*mass*nz(a)*nz(a))*pi*pi;
        hbnz2=h_bar*h_bar/(mass*nz(a)*nz(a));
        
        xi=         zeros(len3);
        yj=         zeros(len3);
        zj=         zeros(len3);
        
for k=1:len3        
      
        xi(k)=         floor((k-1)/(lenz*leny))+1;
        yj(k)=         mod(floor((k-1)/lenz),leny)+1;% mod(k-1,lenx)+1;
        zj(k)=         mod((k-1),lenz)+1;  
    
end
    
for k=1:len3
    
        %parametrical equations for xi and yj
        
    for k_prime=1:len3
        %parametrical equations for xi_prime and yj_prime
 
   %Kinetic Energy T1
      if yj(k)==yj(k_prime) && zj(k)==zj(k_prime)  
            if xi(k)==xi(k_prime)
                T1=hbnx1;
            else
                T1=(-1)^(xi(k)-xi(k_prime))*hbnx2/((xi(k)-xi(k_prime))*(xi(k)-xi(k_prime)));
            end
      else
            T1=0; %kronecker delta
        
      end
        %kinetic energy T2
        
      if xi(k)==xi(k_prime) && zj(k)==zj(k_prime)  
            if yj(k)==yj(k_prime)
                   T2=hbny1;
            else
                T2=(-1)^(yj(k)-yj(k_prime))*hbny2/((yj(k)-yj(k_prime))*(yj(k)-yj(k_prime)));
            end
     else
            T2=0; %kronecker delta
     end  
     
     %kinetic energy T3
      
      if xi(k)==xi(k_prime) && yj(k)==yj(k_prime)  
          if zj(k)==zj(k_prime)
                   T3=hbnz1;
          else
               T3=(-1)^(zj(k)-zj(k_prime))*hbnz2/((zj(k)-zj(k_prime))*(zj(k)-zj(k_prime)));
          end
        
      else
            T3=0; 
      end     
      
%      if (((x(xi)^2)/(lx/2)^2+(y(yj)^2)/(ly/2)^2+(z(zj)^2)/(lz/2)^2)>1 && xi==xi_prime && yj==yj_prime && zj==zj_prime)
%          V=p;
%      else
%          V=0;
%      end

        %H = T + V
%         H1=T1+T2+T3;
%          if H1 ~= 0
%             H(k,k_prime) = H1;
%             count=count+1;
%          end
        H1=T1+T2+T3;
        if (H1) ~= 0
            count=count+1;
            aa(count)=k;
            bb(count)=k_prime;
            cc(count)=H1;
        end
     end
end     
     
        % computations take place here
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
case 3 %rods
        hbnx1=h_bar*h_bar/(6*mass*nx(a)*nx(a))*pi*pi;
        hbnx2=h_bar*h_bar/(mass*nx(a)*nx(a));
        hbny1=h_bar*h_bar/(6*mass*ny(a)*ny(a))*pi*pi;
        hbny2=h_bar*h_bar/(mass*ny(a)*ny(a));
        hbnz1=h_bar*h_bar/(6*mass*nz(a)*nz(a))*pi*pi;
        hbnz2=h_bar*h_bar/(mass*nz(a)*nz(a));
        
        xi=         zeros(len3,1);
        yj=         zeros(len3,1);
        zj=         zeros(len3,1);
        
for k=1:len3        
      
        xi(k)=         floor((k-1)/(lenz*leny))+1;
        yj(k)=         mod(floor((k-1)/lenz),leny)+1;% mod(k-1,lenx)+1;
        zj(k)=         mod((k-1),lenz)+1;  
        
          if ((x(xi(k))*x(xi(k)))/(lx*lx/4)+(y(yj(k))*y(yj(k)))/(ly*ly/4))>1 
          V(k)=p;
      end

end
    


    for k=1:len3
    
        %parametrical equations for xi and yj
%         holder= zeros(1,len3);
        
    for k_prime=1:len3
        %parametrical equations for xi_prime and yj_prime
 
   %Kinetic Energy T1 and V
      if yj(k)==yj(k_prime) && zj(k)==zj(k_prime)  
            if xi(k)==xi(k_prime)
                T1=hbnx1+V(k);
            else
                T1=(-1)^(xi(k)-xi(k_prime))*hbnx2/((xi(k)-xi(k_prime))*(xi(k)-xi(k_prime)));
            end
      else
            T1=0; %kronecker delta
        
      end
        %kinetic energy T2
        
      if xi(k)==xi(k_prime) && zj(k)==zj(k_prime)  
            if yj(k)==yj(k_prime)
                   T2=hbny1;
            else
                T2=(-1)^(yj(k)-yj(k_prime))*hbny2/((yj(k)-yj(k_prime))*(yj(k)-yj(k_prime)));
            end
     else
            T2=0; %kronecker delta
     end  
     
     %kinetic energy T3
      
      if xi(k)==xi(k_prime) && yj(k)==yj(k_prime)  
          if zj(k)==zj(k_prime)
                   T3=hbnz1;
          else
               T3=(-1)^(zj(k)-zj(k_prime))*hbnz2/((zj(k)-zj(k_prime))*(zj(k)-zj(k_prime)));
          end
        
      else
            T3=0; 
      end      
      

        %H = T + V   
        H1 = T1 + T2 + T3;
         if H1 ~= 0
            count=count+1;
            aa(count)=k;
            bb(count)=k_prime;
            cc(count)=H1;
         end
      
    end

     
 end    
  
  otherwise
    disp('Enter either 1 for ellipsoid, 2 for rectangular prism, or 3 for rod.\n')
   end


 H=sparse(aa(aa~=0),bb(bb~=0),cc(cc~=0),len3,len3);
b=size(H); 
count = nnz(H);
bytes = nzmax(H)/10^6;
sprse=count/(b(1)*b(1))*100;
[EVec,E3]=eigs(H,1,'sm');


E2 = diag(E3);
E2 = sort(E2);
Eactual=(h_bar^2/(2*mass*(lz/2)^2))*(2.81718+6.47845*((lz/2)/(lx/2)+0.04372)^1.96084);

perror=((Eactual-E2)/Eactual)*100;

c1=toc;
%--------------------------------------------------------------------------
tic

coulomb=0;

for k=1:len3
    for k_prime=1:len3
        %parametrical equations for xi_prime and yj_prime


       if (x(xi(k))-x(xi(k_prime)))^2+(y(yj(k))-y(yj(k_prime)))^2+(z(zj(k))-z(zj(k_prime)))^2~=0;
          coulomb=coulomb+EVec(k)^2*EVec(k_prime)^2/(sqrt((x(xi(k))-x(xi(k_prime)))^2+(y(yj(k))-y(yj(k_prime)))^2+(z(zj(k))-z(zj(k_prime)))^2));
       end
     end
end  

e = 1.602e-19;
ecdse = 10.6;
enot = 8.854e-12;
Ebulk = 2.805e-19;

Ecoulomb=(coulomb)*e^2/(4*pi*ecdse*enot);
EcA=e^2/(4*pi*ecdse*enot*(lz/2))*(1.77404+1.11755*abs(((lz/2)/(lx/2))-1.00962).^0.85486) ;
perrorEc=((EcA-Ecoulomb)/EcA)*100;

Etot = Ebulk+E2-Ecoulomb;
EtotA= Ebulk+Eactual-EcA;
perrorTot=((EtotA-Etot)/EtotA)*100;

Etot2 = Ebulk+E2; 

lam = (2*pi*h_bar*2.99792458e8/Etot)*1e9;
lam2 = (2*pi*h_bar*2.99792458e8/Etot2)*1e9;

beep
c2=toc;
disp( ['=================================================================='])
disp( ['Comments: sparse with no endpoints. hal server.'])
disp( ['Comments: errors are meaningless for rectangle and rod.'])
disp( [' nx              ny              nz '])
fprintf('%e \t %e \t %e \n\n',nx(a),ny(a),nz(a))
disp( [' lx              ly              lz '])
fprintf('%e \t %e \t %e \n\n',lx,ly,lz)
disp( [' size            sparse %        count '])
fprintf( '%f \t %f \t %f \n\n', b(1), sprse, count)
disp( [' E2              Ecoulomb        Etotal'])
fprintf( '%e \t %e \t %e\n\n', E2, Ecoulomb, Etot)
disp( [' Eactual         EcA             EtotA'])
fprintf( '%e \t %e \t %e \n\n', Eactual, EcA, EtotA)
disp( [' Perror          PerrorC         PerrorTot'])
fprintf( '%f \t %f \t %f \n\n', perror, perrorEc, perrorTot)
disp( [' lam             lam2 '])
fprintf( '%e \t %e \n\n', lam, lam2)
disp( [' time1           time2 '])
fprintf( '%f \t %f \n\n', c1, c2)
diary( 'QDDataSphere101.txt')
clearvars  -except a 

end

