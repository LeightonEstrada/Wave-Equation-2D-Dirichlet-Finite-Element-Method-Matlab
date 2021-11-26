%ecuacion calor 2d elementos finitos 2d
%condiciones newmann
%leighton estrada rayme 04140403 unmsm fcm cc 14.5 
clear all;clc;
%-------------------------------------------------------------------------%
%AJUSTES INICIALES
%-------------------------------------------------------------------------%
tmax=5;%tiempo maximo
m=500;%numero particiones parte temporal
tt=tmax/m;%incremento tiempo
alfa=1;%dato inicial
%alfa=1;
%-------------------------------------------------------------------------%
%MALLADO TRIANGULAR
%-------------------------------------------------------------------------%
n=10;%numero particiones espaciales
h=1/n;ef=4*n^2;pu=(n+1)^2+n^2;%ef puntos %pu=nodos mallas
for j=1:n
    for i=(j-1)*n+1:j*n
        t(4*i-3,1)=i+j-1;t(4*i-2,1)=i+j-1;t(4*i-1,1)=i+j;t(4*i,1)=(n+1)^2+i;
        t(4*i-3,2)=(n+1)^2+i;t(4*i-2,2)=i+j;t(4*i-1,2)=n+1+i+j;t(4*i,2)=n+1+i+j;
        t(4*i-3,3)=n+i+j;t(4*i-2,3)=(n+1)^2+i;t(4*i-1,3)=(n+1)^2+i;t(4*i,3)=n+i+j;
    end
end
for j=1:n+1
    for i=(j-1)*(n+1)+1:j*(n+1)
        p(i,1)=(i-1-(j-1)*(n+1))*h;p(i,2)=(j-1)*h;
    end
end
for j=1:n
    for i=(j-1)*n+1:j*n
        pp(i,1)=(i-1-(j-1)*n)*h+(h/2);pp(i,2)=(j-1)*h+(h/2);
    end
end
p((n+1)^2+1:(n+1)^2+n^2,:)=pp;
%-------------------------------------------------------------------------%
%ELEMENTOS FINITOS 2D
%-------------------------------------------------------------------------%
ind=((p(:,1)>0 & p(:,1)<1)&(p(:,2)>0 & p(:,2)<1));%obtiene puntos interiores caso general dom cuadrado
Np=size(p,1);N=sum(ind);% Np=# de nodos; N=# de nodos interiores
in=zeros(Np,1);in(ind)=(1:N)';
K=zeros(N,N);%matriz de rigidez
M=zeros(N,N);%matriz de masa
for i=1:size(t,1)
j=t(i,1);k=t(i,2);l=t(i,3);
vj=in(j);vk=in(k);vl=in(l);
J=[p(k,1)-p(j,1), p(l,1)-p(j,1); p(k,2)-p(j,2), p(l,2)-p(j,2)];
ar=abs(det(J))/2; ar1=abs(det(J))/24;CC=ar/12; Q=inv(J'*J);
if vj>0 
    K(vj,vj)=K(vj,vj)+ar*sum(sum(Q));
    M(vj,vj)=M(vj,vj)+ar1*2;
end
if vk>0
    K(vk,vk)=K(vk,vk)+ar*Q(1,1);
    M(vk,vk)=M(vk,vk)+ar1*2;
end
if vl>0
    K(vl,vl)=K(vl,vl)+ar*Q(2,2);
    M(vl,vl)=M(vl,vl)+ar1*2;
end
if vj*vk>0
    K(vj,vk)=K(vj,vk)-ar*sum(Q(:,1)); K(vk,vj)=K(vj,vk);
    M(vj,vk)=M(vj,vk)+ar1; M(vk,vj)=M(vj,vk);
end
if vj*vl>0
    K(vj,vl)=K(vj,vl)-ar*sum(Q(:,2)); K(vl,vj)=K(vj,vl);
    M(vj,vl)=M(vj,vl)+ar1; M(vl,vj)=M(vj,vl);
end
if vk*vl>0
    K(vk,vl)=K(vk,vl)+ar*Q(1,2); K(vl,vk)=K(vk,vl);
    M(vk,vl)=M(vk,vl)+ar1; M(vl,vk)=M(vk,vl);
end
end
%-------------------------------------------------------------------------%
%CONDICIONES INICIALES EXPERIMENTAL
%-------------------------------------------------------------------------%
U=zeros(pu,m+1);
V=zeros(pu,m+1);
ff=inline('x.*y.*(x-1).*(y-1)');%membrana vibrante
%ff=inline('5*(x-0.5).*exp(-(25*(x-0.5).^2+25*(y-0.5).^2))');
%ff=inline('5*(x-0.5).*exp(-(25*(x-0.5).^2))+0.*y');
%ff=inline('exp(-(25*(x-0.5).^2+25*(y-0.5).^2))');
U(:,1)=ff(p(:,1),p(:,2));%condicion inicial
%-------------------------------------------------------------------------%
UU=U;VV=V;Ne=Np-N;
indx=~((p(:,1)>0 & p(:,1)<1)&(p(:,2)>0 & p(:,2)<1));%obtiene puntos frontera caso general dom cuadrado
ix=(1:Np)';ii=ix(indx);
for i=1:Ne
    iix(i)=ii(Ne-i+1);
end
for i=iix
    UU(i,:)=[];%elimina filas y columnas de los indices de los puntos frontera
    VV(i,:)=[];
end
%-------------------------------------------------------------------------%
%CRANK NICHOLSON
%-------------------------------------------------------------------------%
gamma=(-alfa/2)*tt;
theta=tt/2;
MK1=(eye(N)-gamma*theta*M^-1*K);
MK2=(eye(N)+gamma*theta*M^-1*K);
MK3=2*gamma*M^-1*K;
for i=2:m+1
    VV(:,i)=MK1^-1*(MK2*VV(:,i-1)+MK3*UU(:,i-1));
    UU(:,i)=theta*(VV(:,i)+VV(:,i-1))+UU(:,i-1);       
end
%adicion a la matriz U los datos dirichlet
U(ind,:)=UU;%esto es importante
V(ind,:)=VV;
%-------------------------------------------------------------------------%
%ANIMACION 3D
%-------------------------------------------------------------------------%
ttt=0:tt:tmax;
for j=1:length(ttt)
    h=sprintf('%5.2f',ttt(j));
    set(gcf,'renderer','zbuffer');
    set(gca,'nextplot','replacechildren');
    caxis manual;caxis([min(min(U)) max(max(U))]);trisurf(t,p(:,1),p(:,2),U(:,j));title({['Membrana Vibrante (Leighton Estrada R. CC FCM UNMSM)'];[num2str(ef),' elementos triangulares y ',num2str(pu),' nodos de malla'];['Tiempo t = ',num2str(h),' s']});xlabel('x');ylabel('y');zlabel('U(x,y;t)');axis([0 1 0 1 min(min(U)) max(max(U))]);
    %view(0,90);shading interp;
    XYZ(j)=getframe;
end
%fin