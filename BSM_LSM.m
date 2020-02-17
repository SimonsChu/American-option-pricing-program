%%
%~~~~~~~~~~~~~~~~~~~~~~~~~Part1: 定义变量~~~~~~~~~~~~~~~~~~~~~~·~~~~~

clear all;
randn('state',0);
tic
S0 = 36;      %初始股票价格
k  = 40;      %执行价格
rf  = 0.06;  %无风险利率
q=0;
T  = 1;       %到期时间1年
sigma=0.2;     %隐含波动率
K  =160;       %时间节点数
N  =10^5;    %路径数量
dt = T/K;      %时间步长
%%



%%
%~~~~~~~~~~~~~~Part2:构造股票价格（St）和现金流（Cash flows）矩阵~~~~~~~~~~~~~~~~
W=normrnd(0,1,N,K); %定义变量
S_tn=zeros(N,K);
S_t0=zeros(N,1);
S_t0(:,1)=S0;
St=[S_t0,S_tn]; 
clear('S_t0');
for c = 1 : 1 : K  %利用for循环构造股票价格矩阵St
    St(:,c+1)=St(:,c).*(exp((rf-q-(sigma^2)*0.5).*dt+sigma*sqrt(dt).*W(:,c)));
end
clear('W');
Cashflows=zeros(N,K);
Cashflows(:,end)=max(k-St(:,end),0); %构造Cash flows 矩阵
%%



%%
%~~~~~~~~Part3:利用LSM算法计算并统计Stopping time，得到最终的Cash flows矩阵~~~~~~~~~
for b =  (K-1) : -1 : 1 %对应时刻b
    Y1=zeros(N,1); %定义变量
    X1=zeros(N,1);
    basis_functions=zeros(4,N);
    ContinValue=zeros(N,1);
    In_the_money_path = find(St(:,b+1)<k);%找到b+1时刻in the money的路径

    %~~~~~~~~~~~Part3.1:对于K-1时间节点，只需将K时刻Cash flows进行折现~~~~~~~~~~~~~~
    if b == (K-1)
        Y1=Cashflows(In_the_money_path,(b+1)).*exp(-rf*dt);  %将b+1时刻的cash flows折现作为回归的因变量
        X1=St(In_the_money_path,b+1)/S0;  %获取b+1时刻的股票价格数据，作为回归的自变量
        A1=polyfit(X1,Y1,6);
        ContinValue=polyval(A1,X1);
        clear('A1');
      
         %基于数学公式推导得到系数表达式，并由回归方程计算continue_value
        
        exercised_path = find((k-St(In_the_money_path,b+1) ) > ContinValue); %找到b+1时刻股价矩阵中立刻行权较优的路径 
        clear('ContinValue');
        Cashflows(In_the_money_path(exercised_path),b) = k-St(In_the_money_path(exercised_path),b+1);  % 将此时的cash flows对应为立刻行权价值
        Cashflows(In_the_money_path(exercised_path),(b+1):K)=0;  %将此时刻之后的cash flows赋值为0  
    %~~~~~~~~~Part3.2:对于之前的节点，需要统计之前最近一次出现cash flows的位置~~~~~~~~~~~
    else 
        Y1=zeros(size(In_the_money_path));
        for d = (b+1) : 1 : K
            last_cash_flows=find(Cashflows(In_the_money_path,d)~=0);   %将b时刻之后cash flows不为零的路径筛选出来
            Y1(last_cash_flows,1)=Cashflows(In_the_money_path(last_cash_flows),d).*exp(-rf*(d-b)*dt);  %将这些时刻的cash flows折现作为回归的因变量
        end
        X1=St(In_the_money_path,b+1)/S0;  %获取b+1时刻的股票价格数据，作为回归的自变量
         A1=polyfit(X1,Y1,5);
        ContinValue=polyval(A1,X1);
        clear('A1');
        exercised_path = find((k-St(In_the_money_path,b+1) ) > ContinValue); %找到b+1时刻股价矩阵中立刻行权较优的路径 
        clear('ContinValue');
        Cashflows(In_the_money_path(exercised_path),b) = k-St(In_the_money_path(exercised_path),b+1);  % 将此时的cash flows对应为立刻行权价值
        Cashflows(In_the_money_path(exercised_path),(b+1):K)=0;  %将此时刻之后的cash flows赋值为0  
    end
end
%%



%%
%~~~~~~~~~~~~~~~~~~~~~~Part4:计算最终的期权价值~~~~~~~~~~~~~~~~~~~~~~~~~  
path_option_value1=zeros(N,1);
for e = 1 : 1 : K  %按照时间节点来依次寻找stopping time
    Stopping_time = find( Cashflows(:,e) ~= 0 );  %找到之前计算的stopping time
    path_option_value1(Stopping_time,1)=Cashflows(Stopping_time,e).*exp(-rf*e*dt); %对于每个cash flows进行折现
end
path_option_value=sum(path_option_value1);%计算得到N条路径的期权价格
Final_Option_value=path_option_value/N  %取平均的得到期权价格
error=(Final_Option_value-4.478)/4.478
standard_deviation=std(path_option_value1)/sqrt(N)
toc
%%

