%%
%~~~~~~~~~~~~~~~~~~~~~~~~~Part1: �������~~~~~~~~~~~~~~~~~~~~~~��~~~~~

clear all;
randn('state',0);
tic
S0 = 36;      %��ʼ��Ʊ�۸�
k  = 40;      %ִ�м۸�
rf  = 0.06;  %�޷�������
q=0;
T  = 1;       %����ʱ��1��
sigma=0.2;     %����������
K  =160;       %ʱ��ڵ���
N  =10^5;    %·������
dt = T/K;      %ʱ�䲽��
%%



%%
%~~~~~~~~~~~~~~Part2:�����Ʊ�۸�St�����ֽ�����Cash flows������~~~~~~~~~~~~~~~~
W=normrnd(0,1,N,K); %�������
S_tn=zeros(N,K);
S_t0=zeros(N,1);
S_t0(:,1)=S0;
St=[S_t0,S_tn]; 
clear('S_t0');
for c = 1 : 1 : K  %����forѭ�������Ʊ�۸����St
    St(:,c+1)=St(:,c).*(exp((rf-q-(sigma^2)*0.5).*dt+sigma*sqrt(dt).*W(:,c)));
end
clear('W');
Cashflows=zeros(N,K);
Cashflows(:,end)=max(k-St(:,end),0); %����Cash flows ����
%%



%%
%~~~~~~~~Part3:����LSM�㷨���㲢ͳ��Stopping time���õ����յ�Cash flows����~~~~~~~~~
for b =  (K-1) : -1 : 1 %��Ӧʱ��b
    Y1=zeros(N,1); %�������
    X1=zeros(N,1);
    basis_functions=zeros(4,N);
    ContinValue=zeros(N,1);
    In_the_money_path = find(St(:,b+1)<k);%�ҵ�b+1ʱ��in the money��·��

    %~~~~~~~~~~~Part3.1:����K-1ʱ��ڵ㣬ֻ�轫Kʱ��Cash flows��������~~~~~~~~~~~~~~
    if b == (K-1)
        Y1=Cashflows(In_the_money_path,(b+1)).*exp(-rf*dt);  %��b+1ʱ�̵�cash flows������Ϊ�ع�������
        X1=St(In_the_money_path,b+1)/S0;  %��ȡb+1ʱ�̵Ĺ�Ʊ�۸����ݣ���Ϊ�ع���Ա���
        A1=polyfit(X1,Y1,6);
        ContinValue=polyval(A1,X1);
        clear('A1');
      
         %������ѧ��ʽ�Ƶ��õ�ϵ�����ʽ�����ɻع鷽�̼���continue_value
        
        exercised_path = find((k-St(In_the_money_path,b+1) ) > ContinValue); %�ҵ�b+1ʱ�̹ɼ۾�����������Ȩ���ŵ�·�� 
        clear('ContinValue');
        Cashflows(In_the_money_path(exercised_path),b) = k-St(In_the_money_path(exercised_path),b+1);  % ����ʱ��cash flows��ӦΪ������Ȩ��ֵ
        Cashflows(In_the_money_path(exercised_path),(b+1):K)=0;  %����ʱ��֮���cash flows��ֵΪ0  
    %~~~~~~~~~Part3.2:����֮ǰ�Ľڵ㣬��Ҫͳ��֮ǰ���һ�γ���cash flows��λ��~~~~~~~~~~~
    else 
        Y1=zeros(size(In_the_money_path));
        for d = (b+1) : 1 : K
            last_cash_flows=find(Cashflows(In_the_money_path,d)~=0);   %��bʱ��֮��cash flows��Ϊ���·��ɸѡ����
            Y1(last_cash_flows,1)=Cashflows(In_the_money_path(last_cash_flows),d).*exp(-rf*(d-b)*dt);  %����Щʱ�̵�cash flows������Ϊ�ع�������
        end
        X1=St(In_the_money_path,b+1)/S0;  %��ȡb+1ʱ�̵Ĺ�Ʊ�۸����ݣ���Ϊ�ع���Ա���
         A1=polyfit(X1,Y1,5);
        ContinValue=polyval(A1,X1);
        clear('A1');
        exercised_path = find((k-St(In_the_money_path,b+1) ) > ContinValue); %�ҵ�b+1ʱ�̹ɼ۾�����������Ȩ���ŵ�·�� 
        clear('ContinValue');
        Cashflows(In_the_money_path(exercised_path),b) = k-St(In_the_money_path(exercised_path),b+1);  % ����ʱ��cash flows��ӦΪ������Ȩ��ֵ
        Cashflows(In_the_money_path(exercised_path),(b+1):K)=0;  %����ʱ��֮���cash flows��ֵΪ0  
    end
end
%%



%%
%~~~~~~~~~~~~~~~~~~~~~~Part4:�������յ���Ȩ��ֵ~~~~~~~~~~~~~~~~~~~~~~~~~  
path_option_value1=zeros(N,1);
for e = 1 : 1 : K  %����ʱ��ڵ�������Ѱ��stopping time
    Stopping_time = find( Cashflows(:,e) ~= 0 );  %�ҵ�֮ǰ�����stopping time
    path_option_value1(Stopping_time,1)=Cashflows(Stopping_time,e).*exp(-rf*e*dt); %����ÿ��cash flows��������
end
path_option_value=sum(path_option_value1);%����õ�N��·������Ȩ�۸�
Final_Option_value=path_option_value/N  %ȡƽ���ĵõ���Ȩ�۸�
error=(Final_Option_value-4.478)/4.478
standard_deviation=std(path_option_value1)/sqrt(N)
toc
%%

