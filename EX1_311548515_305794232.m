close all; clear all;
%%
%loading the data
load('S1.mat'); 
load('S2.mat');
%%
% Qualitative observation of the given signals 
Si=S2; % for generalization of the code
fs = 10000; %sampling rate
dt = 1/fs; %time step
N = length(Si); %number of components in Si
t = dt*(1:N);%time vector
%plot(t,Si);xlabel('sec'); 

%%
%Finding the spike times 
TH =-30; % defining treshold
saTH =(Si>TH); % binary vector-"above treshold"
saTH_diff = diff(saTH); %marking the state changes of SaTH
L2H=find(saTH_diff>0);%crossing the treshold from below (events are the positions of 1)
H2L=find(saTH_diff<0);%crossing the treshold from above (events are the positions of -1) 
for i=1:length(H2L)
    [LM_max(i),LM(i)]= max(Si((L2H(i):H2L(i))));% creating the max spikes value and indexes vectors between H2L & L2H 
    LM(i)= LM(i)+L2H(i)-1;%updating the spike's time vector
end

%% 
%Finding the spike rate per segment
n_segments=length(t)/10000/0.3; % number of segments 
m_segments=(reshape(t,3000,n_segments));% organizing "t" in a matrix wich each column is a segment of 0.3 s 
spiket_comp=ismember(m_segments,LM*dt);% to find how many max values (acoording to LM) are in each column  
SC=sum(spiket_comp,1); % number of spikes in each segment
R=SC/0.2; %firing rate for each segment (0.2s the duration of voltage) 

%% 
%plot
hold on 
plot(t,Si,'b'); % a plot presenting voltage by time
title('\bf Intra-Cellular Recording')%plot title 
xlabel({'\bf Time','(In Seconds)'})% axsis x label
ylabel({'\bf Voltage','(In mv)'})% axsis y label
set(t,Si,'o','MarkerIndices',L2H,'MarkerEdgeColor','g'); % mark the L2H on the first plot 
set(t,Si,'o','MarkerIndices',H2L,'MarkerEdgeColor','r');% mark the H2L on the first plot
set(t,Si,'o','MarkerIndices',LM,'MarkerEdgeColor','k');% mark the LM on the first plot
lgd=legend('Voltage-Time','L2H','H2L','LM (max spikes)'); % plot legend
lgd.FontSize = 10; % font size of legend 
lgd.FontWeight = 'bold'; % bold font for the legend
legend('Location','northeastoutside')% location of legend
seg_step=m_segments(1,:)-dt; % time vector wich seprated into 0.3s parts 
for q=1:length(seg_step) % ceating headlines for the firing rate
    text(2,18,'The Rate Per Segment:')
    text(seg_step(q),10,num2str(R(q)))
end
hold off;






