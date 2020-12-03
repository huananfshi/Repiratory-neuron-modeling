%parameters
C = 21;
g_Na_max = 28;
g_K_max = 11.2;
g_NaP_max = 2.4;
gL_K = 2.4;
gL_Na = 0.4;
E_Na = 50;
E_K = -75;
theta_m = -34;
theta_n = -29;
theta_h = -48;
theta_p = -40;
s_m = -5;
s_n = -4;
s_h = 6;
s_p = -6;
tau_n_max = 10;
tau_h_max = 10^4;
%single neuron
T = 4*10^4;
dt = 0.01;
N = T/dt;
n = zeros(1,N);
h = n;
V_m = n;
V_m(1) = -56;
n_max = 1/(1+exp((V_m(1)-theta_n)/s_n));
tau_n = tau_n_max/(cosh((V_m(1)-theta_n)/2/s_n));
h_max = 1/(1+exp((V_m(1)-theta_h)/s_h));
tau_h = tau_h_max/(cosh((V_m(1)-theta_h)/2/s_h));

for j = 2:N
    dn = dt*(n_max-n(j-1))/tau_n;
    dh = dt*(h_max-h(j-1))/tau_h;
    n(j) = dn+n(j-1);
    h(j) = dh+h(j-1);
    m_max = 1/(1+exp((V_m(j-1)-theta_m)/s_m));
    I_Na = g_Na_max*m_max^3*(1-n(j))*(V_m(j-1)-E_Na);
    I_K = g_K_max*n(j)^4*(V_m(j-1)-E_K);
    p_max = 1/(1+exp((V_m(j-1)-theta_p)/s_p));
    I_NaP = g_NaP_max*p_max*h(j)*(V_m(j-1)-E_Na);
    IL_K = gL_K*(V_m(j-1)-E_K);
    IL_Na = gL_Na*(V_m(j-1)-E_Na);
    dV = -dt*(I_NaP+I_Na+I_K+IL_K+IL_Na)/C;
    V_m(j)=V_m(j-1)+dV;
    n_max = 1/(1+exp((V_m(j)-theta_n)/s_n));
    tau_n = tau_n_max/(cosh((V_m(j)-theta_n)/2/s_n));
    h_max = 1/(1+exp((V_m(j)-theta_h)/s_h));
    tau_h = tau_h_max/(cosh((V_m(j)-theta_h)/2/s_h));
end
t = linspace(0,T/1000,N);
subplot(1,2,1)
plot(t,V_m);
title('E_K = 75 mV')
xticks(0:5:40)
yticks(-60:20:0)
xlabel('t (s)')
ylabel('V_m (mV)')
box off
subplot(1,2,2)
t1 = 10*10^3/dt;
t2 = 15*10^3/dt;
plot(t(t1:1:t2),V_m(t1:1:t2));
title('E_K = 75 mV')
xticks(10:1:15)
xlabel('t (s)')
xlim([11 14])
yticks([])
box off
%Pacemaker network
figure(2)
N_c = 50;
s = zeros(1,N_c);
g_syn = normrnd(0.1,0.025,1,N_c);
g_tonic = 0.035;
theta_s = -10;
tau_s = 5;
kr = 1;
s_s = -5;
E_syn = 0;
n = zeros(1,N);
h = n;
V_m = n;
V_m(1) = -56;
V_m_f =V_m;
V_m_f(1) = -60;
n_max = 1/(1+exp((V_m(1)-theta_n)/s_n));
tau_n = tau_n_max/(cosh((V_m(1)-theta_n)/2/s_n));
h_max = 1/(1+exp((V_m(1)-theta_h)/s_h));
tau_h = tau_h_max/(cosh((V_m(1)-theta_h)/2/s_h));

for j = 2:N
    dn = dt*(n_max-n(j-1))/tau_n;
    dh = dt*(h_max-h(j-1))/tau_h;
    n(j) = dn+n(j-1);
    h(j) = dh+h(j-1);
    m_max = 1/(1+exp((V_m(j-1)-theta_m)/s_m));
    I_Na = g_Na_max*m_max^3*(1-n(j))*(V_m(j-1)-E_Na);
    I_K = g_K_max*n(j)^4*(V_m(j-1)-E_K);
    p_max = 1/(1+exp((V_m(j-1)-theta_p)/s_p));
    I_NaP = g_NaP_max*p_max*h(j)*(V_m(j-1)-E_Na);
    IL_K = gL_K*(V_m(j-1)-E_K);
    IL_Na = gL_Na*(V_m(j-1)-E_Na);
    s_max = 1/(1+exp((V_m(j-1)-theta_s)/s_s));
    for i = 1:N_c
        ds = dt*((1-s(i))*s_max-kr*s(i))/tau_s;
        s(i)=ds+s(i);
    end
    I_syn_f = binornd(1,0.5)*g_f*s(N_c)*(V_m_f(j-1)-E_syn);
    dV_f = -dt*(I_Na+I_K+IL_K+IL_Na+I_syn_f)/C;
    V_m_f(j)=V_m_f(j-1)+dV;
    I_syn = sum(g_syn(1:N_c-1).*s(1:N_c-1))*(V_m(j-1)-E_syn);
    I_tonic = g_tonic*(V_m(j-1)-E_syn);
    dV = -dt*(I_NaP+I_Na+I_K+IL_K+IL_Na+I_syn+I_tonic)/C;
    V_m(j)=V_m(j-1)+dV;
    n_max = 1/(1+exp((V_m(j)-theta_n)/s_n));
    tau_n = tau_n_max/(cosh((V_m(j)-theta_n)/2/s_n));
    h_max = 1/(1+exp((V_m(j)-theta_h)/s_h));
    tau_h = tau_h_max/(cosh((V_m(j)-theta_h)/2/s_h));
end
t = linspace(0,T/1000,N);
subplot(2,1,1)
plot(t,V_m);
title('Pacemaker neuron X')
xticks(0:5:40)
yticks(-60:20:0)
xlabel('t (s)')
ylabel('V_m (mV)')
box off
subplot(2,1,2)
plot(t,V_m_f)
xticks(0:5:40)
yticks(-60:20:0)
xlabel('t (s)')
ylabel('V_m (mV)')
title('Follower non-pacemaker neuron X')
box off
%with inhibitory input
figure(3)
T = 2*10^4;
N = T/dt;
n = zeros(1,N);
h = n;
V_m = n;
V_m(1) = -56;
V_m_l = n;
V_m_l(1) = -65;
V_m_f =n;
V_m_f(1) = -60;
n_max = 1/(1+exp((V_m(1)-theta_n)/s_n));
tau_n = tau_n_max/(cosh((V_m(1)-theta_n)/2/s_n));
h_max = 1/(1+exp((V_m(1)-theta_h)/s_h));
tau_h = tau_h_max/(cosh((V_m(1)-theta_h)/2/s_h));
s_l = 0;
for j = 2:N
    dn = dt*(n_max-n(j-1))/tau_n;
    dh = dt*(h_max-h(j-1))/tau_h;
    n(j) = dn+n(j-1);
    h(j) = dh+h(j-1);
    m_max = 1/(1+exp((V_m(j-1)-theta_m)/s_m));
    I_Na = g_Na_max*m_max^3*(1-n(j))*(V_m(j-1)-E_Na);
    I_K = g_K_max*n(j)^4*(V_m(j-1)-E_K);
    p_max = 1/(1+exp((V_m(j-1)-theta_p)/s_p));
    I_NaP = g_NaP_max*p_max*h(j)*(V_m(j-1)-E_Na);
    IL_K = gL_K*(V_m(j-1)-E_K);
    IL_Na = gL_Na*(V_m(j-1)-E_Na);
    s_max = 1/(1+exp((V_m(j-1)-theta_s)/s_s));
    for i = 1:N_c
        ds = dt*((1-s(i))*s_max-kr*s(i))/tau_s;
        s(i)=ds+s(i);
    end
    g_l=binornd(1,0.5)*normrnd(0.1,0.025);
    I_syn_l = g_l*s(N_c)*(V_m_l(j-1)-E_syn);
    I_tonic_l = g_tonic*(V_m(j-1)-E_syn);
    dV_l = -dt*(I_NaP+I_Na+I_K+IL_K+IL_Na-I_syn_l+I_tonic_l)/C;
    V_m_l(j)=V_m_l(j-1)+dV;
    s_l_max = 1/(1+exp((V_m_l(j-1)-theta_s)/s_s));
    ds_l = dt*((1-s_l)*s_l_max-kr*s_l)/tau_s;
    s_l = ds_l+s_l;
    g_f=normrnd(1,0.25,1,2).*binornd(1,0.5,1,2);
    I_syn_f = sum(g_f.*[s(N_c) s_l])*(V_m_f(j-1)-E_syn);
    dV_f = -dt*(I_Na+I_K+IL_K+IL_Na+I_syn_f)/C;
    V_m_f(j)=V_m_f(j-1)+dV;
    I_syn = (sum(g_syn(1:N_c-1).*s(1:N_c-1))-g_l*s_l)*(V_m(j-1)-E_syn);
    I_tonic = g_tonic*(V_m(j-1)-E_syn);
    dV = -dt*(I_NaP+I_Na+I_K+IL_K+IL_Na+I_syn+I_tonic)/C;
    V_m(j)=V_m(j-1)+dV;
    n_max = 1/(1+exp((V_m(j)-theta_n)/s_n));
    tau_n = tau_n_max/(cosh((V_m(j)-theta_n)/2/s_n));
    h_max = 1/(1+exp((V_m(j)-theta_h)/s_h));
    tau_h = tau_h_max/(cosh((V_m(j)-theta_h)/2/s_h));
end
t = linspace(0,T/1000,N);
subplot(3,1,1)
plot(t,V_m);
title('Pacemaker neuron X')
xticks(0:5:20)
yticks(-60:20:0)
xlabel('t (s)')
ylabel('V_m (mV)')
ylim([-70 10])
box off
subplot(3,1,2)
plot(t,V_m_l);
title('Follower pacemaker neuron X')
xticks(0:5:20)
yticks(-60:20:0)
xlabel('t (s)')
ylabel('V_m (mV)')
ylim([-70 10])
box off
subplot(3,1,3)
plot(t,V_m_f)
xticks(0:5:20)
yticks(-60:20:0)
xlabel('t (s)')
ylabel('V_m (mV)')
title('Follower non-pacemaker neuron X')
ylim([-70 10])
box off
figure(4)
t3 = 12*10^3/dt;
t4 = 12.02*10^3/dt;
plot(t(t3:t4),V_m(t3:t4),'k-',t(t3:t4),V_m_f(t3:t4),'r--',t(t3:t4),V_m_l(t3:t4),'b-')
xticks([])
xlabel('t')
yticks(-50:20:10)
ylabel('V_m (mV)')
legend('preBötC','XII_n','pF_L')
box off
