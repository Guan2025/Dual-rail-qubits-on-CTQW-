clear all
close all
T=2*pi;
J=1;
J2=-sqrt(2);

x=-80:1:80;

for m1=5:30
    tmp1=sqrt(m1^2-16)+m1;
    
%     for m2=6:30
%         tmp2=sqrt(m2^2-16)+m2;
    for k=1:161       
        tmp3(k)=abs(tmp1/2+x(k)/2);
    end
    [mintmp,index]=min(tmp3);
    minvalue(m1-4)=tmp1/2+x(index)/2;
%     end
end

for k1=5:30
      V1=minvalue(k1-4);    
      m1=k1;
      U=V1-sqrt(m1^2-16);
   V2_const=sqrt(4^2-4);
    for i=1:201
        tmp(i)=V2_const+(i-1)/100*0.07-0.07;
        V2=tmp(i);
        
     H101=[0,-1;-1,V2];
     H101_2=[0,-1;-1,-V2*1];
     psi101=[1 0]';
     Uevo1=expm(-1i*H101*T);
     Uevo2=expm(-1i*H101_2*T);
     psi101=Uevo2*Uevo1*psi101;
%    psi101=Uevo1*psi101;
     P(i)=psi101(1)'*psi101(1);
     theta(i)=angle(psi101(1))/pi;
     
     psi011=[0 1]';
     psi011=Uevo2*Uevo1*psi011;
     
        
     Htest1=[V1+V2,J2,J2;J2,U,0; J2,0,U+2*V2];
    U1=expm(-1i*T*Htest1);
    V2=-V2;
    Htest2=[V1+V2,J2,J2;J2,U,0; J2,0,U+2*V2];
     U2=expm(-1i*T*Htest2);
     psi=[1 0 0]';
     psi=U2*U1*psi;
     P111(i)=psi(1)*psi(1)';     
     phi111(k1,i)=angle(psi(1))/pi;
     
     M=eye(8);
    M(6,6)=psi101(1);
    M(7,7)=psi011(2);
%    M(8,8)=P111(i);
    M(8,8)=norm(psi(1));
    F(i)=1/72*(trace(M*M')+abs(trace(M))^2);
    end
    [maxF(k1),index]=max(F');
%    [P111_out(k),index(k)]=max(P111);
    phi111_final(k1)=phi111(k1,index);
    
end


[maxmaxF]=max(maxF);
%phi111(:,index(:))
% The fidelity F>=0.99 when 14<=m<=30
% The generated phases on the state 111 are as follows
phi111_final(14:30)