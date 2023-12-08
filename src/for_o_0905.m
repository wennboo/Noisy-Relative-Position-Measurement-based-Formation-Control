%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%�������ĳ�������KKO�����Ʒ�λ�ǣ�
clc;
clear;
close all;
flag=1;
if(flag)
    format long;
else
    format short;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%ͼ�μ������ǶȵĹ�ϵ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A0=[0 0 0 0 0;
%    1 0 0 0 0;
%    0 1 0 0 0;
%    0 0 1 0 0;
%    0 0 0 1 0];
A0=[0 0 0 0 0 0 0 0 1;
    1 0 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 1;
    0 0 0 1 0 0 0 0 0;
    0 0 0 0 1 0 0 0 1;
    0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 0 1 0 0;
    0 1 0 0 0 0 0 1 0];
N0=length(A0);
A=A0+A0';
A=A>0;%%%%%%%%%%%%%%%%%��ɶԳ���
%A0=A
L=(diag(sum(A0,2))-A0);
M=12000;%��������

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%����ϵ��
delta=0.7;
c=1/max(sum(A0,2))-0.1;
c=0.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sigma0=pi/6;%��λ�ǲ�����������
%sigma0=0.017;
sigma0=pi/6;
alpha2=1;%%%%�޲���
%alpha2=sigma0/sin(sigma0);
sigma1=0.01;%������������

%sigma1=0
sigma2=0.1;%ͨ����������
%sigma2=0
%%%%%%%%%��������ϵ��ȫ������ϵ��λ��
rel_theta0=[pi/2,pi/3,pi/4,pi/5,pi/6,pi/7,pi/8,pi/9,pi/10]';
%%%%%%%%%%%%%%��ʼλ��
im=sqrt(-1);
%p0=[0,1,2,3,4,im,2*im,3*im,4*im];
% p0=[0,0.1,0.2,0.3,0.4,0.1*im,0.2*im,0.3*im,0.4*im];
% p0=p0.';
%%%������ʼֵ
%%%��ʼֵ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Ŀ������q0
q0=3*[-1+im,im,1+im,-1,0,1,-1-im,-im,1-im];
q0=q0.';
p0=0.1*q0;
load ini_p0
p0
hat_theta0=zeros(N0,1);
%%%%%%%%%�ھӽڵ�以�����
Ang0=zeros(N0);
P_r0=zeros(N0);
Dis0=zeros(N0);
Or_0=zeros(N0);
Q_0=zeros(N0);
Err=zeros(N0);
for i0=1:N0
    for j0=1:N0
        P_r0(i0,j0)=p0(j0)-p0(i0);
        Ang0(i0,j0)=phase(P_r0(i0,j0));
        Or_0(i0,j0)=rel_theta0(i0)-rel_theta0(j0);
        Dis0(i0,j0)=norm(P_r0(i0,j0));
        Q_0(i0,j0)=q0(j0)-q0(i0);
        Err(i0,j0)=abs(P_r0(i0,j0))-abs(Q_0(i0,j0));
    end
end
error0=sum(sum(abs(Err)));
err=error0;
ro0=exp(-im*rel_theta0);%%%%%%
Ro0=-rel_theta0*ones(1,N0);
L_ang0=Ro0+Ang0;%%%%%%%���������µĽǶ�
v10=unifrnd (-sigma0,sigma0,N0,N0);%%%%%%%%ƽ���ֲ��Ĳ�������
v20=unifrnd (-sigma1,sigma1,N0,N0);
tem0=exp(im*(L_ang0+v10));
B_b0=(Dis0+v20).*tem0;
%Noi0=zeros(N0);
Ang0=A.*(L_ang0+v10);
%%%%%%%%%%mu_ji�����
B_mu0=Ang0'-Ang0;%%%%%%%%%%%��һ�β������
T_mu0=B_mu0;%%%����mu
H_mu0=B_mu0;%%ñmu
H_v10=zeros(N0,1);
H_mu00=B_mu0;
for i0=1:N0
    for j0=1:N0
        tem0=mod(H_mu00(i0,j0)+2*pi,2*pi)-pi;
        tem1=A0(i0,j0)*(Or_0(i0,j0)-tem0);
        H_v10(i0)=H_v10(i0)+tem1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%�������
t_theta0=hat_theta0-rel_theta0
B_b00=alpha2*B_b0;
H_v100=H_v10;
%B_mu00=B_mu0;
E=eye(N0);
D=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%�ڽӾ���Ŀ�Խ���
for i=1:N0;
D=blkdiag(D,A0(i,:));
end
%%%%%%

%%%%%%%%%%%%%%%%%%%%
%%%%%��ʼֵ
t_theta00=t_theta0;
t_theta=t_theta0;
p00=p0;
%p00(1)=5+im*5;
p=p00;
h_v20=sqrt(sigma2)*randn(N0^2,1);
h_v200=h_v20;
err_o=max(t_theta00)-min(t_theta00);
%%%%%%%%%%%%
flag=1;
if(flag)
for t=1:M
    %%%%%%%%%%%%%%%%%�Ƕ�%%%%%%%%%%%%
    at=c;
    t_theta1=(E-at*L)*t_theta00+at*D*h_v200+at*H_v100;
    t_theta=[t_theta t_theta1];%%%%%%�洢����ֵ
    p1=zeros(N0,1);
    tem_err_o=max(t_theta1)-min(t_theta1);
    err_o=[err_o tem_err_o];
    for m1=1:1:N0
        tem1=0;
        for m2=1:1:N0
            tem0=exp(-im*(t_theta00(m1)+rel_theta0(m1)))*Q_0(m1,m2);
            tem0=A0(m1,m2)*(B_b00(m1,m2)-tem0);
            tem1=tem0+tem1;
        end
        p1(m1)=p00(m1)+at*exp(im*rel_theta0(m1))*tem1;
    end
    p=[p p1];%%%%%%%%%%%�洢��ʵλ��
    P_r1=zeros(N0);
    Dis1=zeros(N0);
    Ang1=zeros(N0);
    Err=zeros(N0);
    for i0=1:N0
        for j0=1:N0
            P_r1(i0,j0)=p1(j0)-p1(i0);
            Ang1(i0,j0)=phase(P_r1(i0,j0));
            Dis1(i0,j0)=norm(P_r1(i0,j0));
            Err(i0,j0)=abs(P_r1(i0,j0))-abs(Q_0(i0,j0));
        end
    end
    error0=sum(sum(abs(Err)));
    err=[err error0];
    L_ang1=Ro0+Ang1;%%%%%%%���������µĽǶ�
    v10=unifrnd (-sigma0,sigma0,N0,N0);%%%%%%%%ƽ���ֲ��Ĳ�������
    %v20=sqrt(sigma1)*randn(N0,N0);
    v20=unifrnd (-sigma1,sigma1,N0,N0);
    tem0=exp(im*(L_ang1+v10));
    B_b00=alpha2*(Dis1+v20).*tem0;%%%%%%%%%%��ʵ���λ�ò���
    Ang1=A.*(L_ang1+v10);
    T_mu1=zeros(N0);
%%%%%%%%%%tilde mu_ji�����
    B_mu1=Ang1'-Ang1;%%%%%%%%%%%��λ��bar_mu
%     for i0=1:N0
%         for j0=1:N0
%             tem0=B_mu1(i0,j0)-B_mu0(i0,j0);
%             if abs(tem0)<pi
%                 T_mu1(i0,j0)=B_mu1(i0,j0);
%             elseif tem0>pi
%                 T_mu1(i0,j0)=B_mu1(i0,j0)-2*pi;
%             else
%                 T_mu1(i0,j0)=B_mu1(i0,j0)+2*pi;
%             end
%                 
%         end
%     end
%     %%%%%%%%%hat muji�����
    H_mu00=B_mu1;
    %%%%%%hat v1�����
    H_v100=zeros(N0,1);
    for i0=1:N0
        for j0=1:N0
            tem0=mod(H_mu00(i0,j0)+2*pi,2*pi)-pi;
            tem1=A0(i0,j0)*(Or_0(i0,j0)-tem0);
            H_v100(i0)=H_v100(i0)+tem1;
        end
    end    
    h_v200=sqrt(sigma2)*randn(N0^2,1);%%%%%��λ�ǹ���ͨ������
    t_theta00=t_theta1;
    p00=p1; 
end
end


p_end=p(:,end)
t_theta_end=t_theta(:,end)

Ang1=zeros(N0);
Dis1=zeros(N0);
Dis2=zeros(N0);
P_end=zeros(N0);
for i0=1:N0
    for j0=1:N0
        P_end(i0,j0)=p_end(j0)-p_end(i0);
        Ang1(i0,j0)=phase(P_end(i0,j0));
        Dis1(i0,j0)=norm(P_end(i0,j0));
        Dis2(i0,j0)=norm(q0(j0)-q0(i0));
    end
end
Dis1
Dis2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%����%%%%%%%%%%%%%%%%%%%%%%
p1=real(p);
p2=imag(p);
z0=real(p_end)
z1=imag(p_end)
%%%%%%%%%������λ��
z2=real(p0);
z3=imag(p0);
scatter(z0,z1,'o','r', 'filled' );
hold on;

scatter(z2,z3,'o', 'filled');
legend('The final position','The initial position');
axis([-1.5,2,-1.5,2]);
axis equal;
for i=1:N0
    plot(p1(i,:),p2(i,:),'k');
    hold on;
end
% hold on;
% z4=real(q0);
% z5=imag(q0);
% scatter(z4,z5);
%axis([3,7,3,6]);
% 
figure;
%%%%%%%%%%%%%%%%%%%%%%
%%%�����
plot(err);
figure;
for i=1:N0
    plot(t_theta(i,:),'k');
    hold on;
end
f_c_err=err;
save c_data f_c_err;
save o_data err_o;




