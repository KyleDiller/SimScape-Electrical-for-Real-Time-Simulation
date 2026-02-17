%Parameters for IEEE 123-node test feeder

ft2km = 0.0003048;

%% General data
Ssub = 5000e3; %Substation Power (VA)
VsubH = 115*1000; % Substation high voltage 
VsubL = 4.16 * 1000; % Substation low voltage 
VxfmH = 4.16*1000; % Transformer high voltage
VxfmL = 0.48*1000; % Transformer low voltage
f=60; %frequency

RonBrk=1e-4; % breaker ON resistance 1e-5 Ohms causes some conditioning warning.

loadcoeff = 1; % it can be modified to make the load more light or eavy 

% %% Configuration line data (in km)
% 
% %configuration 1 
% R1 = [0.4576 0.1560 0.1535; 0.1560 0.4666 0.1580; 0.1535 0.1580 0.4615]/1.60934; % 1.60934 converts ohms per mile into ohms per kilometers
% L1 = [1.0780 0.5017 0.3849; 0.5017 1.0482 0.4236; 0.3849 0.4236 1.0651]/(2*pi*f)/1.60934; % L = X/(2*pi*f)
% CAg1 = (5.6765-1.8319-0.6982)*1e-6/(2*pi*f)/1.60934; % C = B/(2*pi*f)
% CBg1 = (5.9809-1.8319-1.1645)*1e-6/(2*pi*f)/1.60934;
% CCg1 = (5.3971-0.6982-1.1645)*1e-6/(2*pi*f)/1.60934;
% CAB1 = 1.8319*1e-6/(2*pi*f)/1.60934;
% CBC1 = 1.1645*1e-6/(2*pi*f)/1.60934;
% CCA1 = 0.6982*1e-6/(2*pi*f)/1.60934;
% 
% %configuration 2
% R2 = [0.4666 0.1580 0.1560; 0.1580 0.4615 0.1535; 0.1560 0.1535 0.4576]/1.60934; % 1.60934 converts ohms per mile into ohms per kilometers
% L2 = [1.0482 0.4236 0.5017; 0.4236 1.0651 0.3849; 0.5017 0.3849 1.0780]/(2*pi*f)/1.60934;
% CAg2 = (5.9809-1.1645-1.8319)*1e-6/(2*pi*f)/1.60934;
% CBg2 = (5.3971-1.1645-0.6982)*1e-6/(2*pi*f)/1.60934;
% CCg2 = (5.6765-1.8319-0.6982)*1e-6/(2*pi*f)/1.60934;
% CAB2 = 1.1645*1e-6/(2*pi*f)/1.60934;
% CBC2 = 0.6982*1e-6/(2*pi*f)/1.60934;
% CCA2 = 1.8319*1e-6/(2*pi*f)/1.60934;
% 
% %configuration 3
% R3 = [0.4615 0.1535 0.1580; 0.1535 0.4576 0.1560; 0.1580 0.1560 0.4666]/1.60934; % 1.60934 converts ohms per mile into ohms per kilometers
% L3 = [1.0651 0.3849 0.4236; 0.3849 1.0780 0.5017; 0.4236 0.5017 1.0482]/(2*pi*f)/1.60934;
% CAg3 = (5.3971-0.6982-1.1645)*1e-6/(2*pi*f)/1.60934;
% CBg3 = (5.6765-0.6982-1.8319)*1e-6/(2*pi*f)/1.60934;
% CCg3 = (5.9809-1.1645-1.8319)*1e-6/(2*pi*f)/1.60934;
% CAB3 = 0.6982*1e-6/(2*pi*f)/1.60934;
% CBC3 = 1.8319*1e-6/(2*pi*f)/1.60934;
% CCA3 = 1.1645*1e-6/(2*pi*f)/1.60934;
% 
% %configuration 4
% R4 = [0.4615 0.1580 0.1535; 0.1580 0.4666 0.1560; 0.1535 0.1560 0.4576]/1.60934; % 1.60934 converts ohms per mile into ohms per kilometers
% L4 = [1.0651 0.4236 0.3849; 0.4236 1.0482 0.5017; 0.3849 0.5017 1.780]/(2*pi*f)/1.60934;
% CAg4 = (5.3971-1.1645-0.6982)*1e-6/(2*pi*f)/1.60934;
% CBg4 = (5.9809-1.1645-1.8319)*1e-6/(2*pi*f)/1.60934;
% CCg4 = (5.6765-0.6982-1.8319)*1e-6/(2*pi*f)/1.60934;
% CAB4 = 1.1645*1e-6/(2*pi*f)/1.60934;
% CBC4 = 1.8319*1e-6/(2*pi*f)/1.60934;
% CCA4 = 0.6982*1e-6/(2*pi*f)/1.60934;
% 
% %configuration 5
% R5 = [0.4666 0.1560 0.1580; 0.1560 0.4576 0.1535; 0.1580 0.1535 0.4615]/1.60934; % 1.60934 converts ohms per mile into ohms per kilometers
% L5 = [1.0482 0.5017 0.4236; 0.5017 1.0780 0.3849; 0.4236 0.3849 1.0651]/(2*pi*f)/1.60934;
% CAg5 = (5.9809-1.8319-1.1645)*1e-6/(2*pi*f)/1.60934;
% CBg5 = (5.6765-1.8319-0.6982)*1e-6/(2*pi*f)/1.60934;
% CCg5 = (5.3971-1.1645-0.6982)*1e-6/(2*pi*f)/1.60934;
% CAB5 = 1.8319*1e-6/(2*pi*f)/1.60934;
% CBC5 = 0.6982*1e-6/(2*pi*f)/1.60934;
% CCA5 = 1.1645*1e-6/(2*pi*f)/1.60934;
% 
% %configuration 6
% R6 = [0.4576 0.1535 0.1560; 0.1535 0.4615 0.1580; 0.1560 0.1580 0.4666]/1.60934; % 1.60934 converts ohms per mile into ohms per kilometers
% L6 = [1.0780 0.3849 0.5017; 0.3849 1.0651 0.4236; 0.5017 0.4236 1.0482]/(2*pi*f)/1.60934;
% CAg6 = (5.6765-0.6982-1.8319)*1e-6/(2*pi*f)/1.60934;
% CBg6 = (5.3971-0.6982-1.1645)*1e-6/(2*pi*f)/1.60934;
% CCg6 = (5.9809-1.8319-1.1645)*1e-6/(2*pi*f)/1.60934;
% CAB6 = 0.6982*1e-6/(2*pi*f)/1.60934;
% CBC6 = 1.1645*1e-6/(2*pi*f)/1.60934;
% CCA6 = 1.8319*1e-6/(2*pi*f)/1.60934;
% 
% %configuration 7(A C phase)
% R7 = [0.4576 0.1535; 0.1535 0.4615]/1.60934; % 1.60934 converts ohms per mile into ohms per kilometers
% L7 = [1.0780 0.3849; 0.3849 1.0651]/(2*pi*f)/1.60934;
% CAg7 = (5.1154-1.0549)*1e-6/(2*pi*f)/1.60934;
% CCg7 = (5.1704-1.0549)*1e-6/(2*pi*f)/1.60934;
% CAC7 = 1.0549*1e-6/(2*pi*f)/1.60934;
% 
% %configuration 8 (A B phase)
% R8 = [0.4576 0.1535; 0.1535 0.4615]/1.60934; % 1.60934 converts ohms per mile into ohms per kilometers
% L8 = [1.0780 0.3849; 0.3849 1.0651]/(2*pi*f)/1.60934;
% CAg8 = (5.1154-1.0549)*1e-6/(2*pi*f)/1.60934;
% CBg8 = (5.1704-1.0549)*1e-6/(2*pi*f)/1.60934;
% CAB8 = 1.0549*1e-6/(2*pi*f)/1.60934;
% 
% %configuration 9 (Single phase A)
% R9 = 1.3292/1.60934; % 1.60934 converts ohms per mile into ohms per kilometers
% L9 = 1.3475/(2*pi*f)/1.60934; 
% C9 = 4.5193*1e-6/(2*pi*f)/1.60934;
% 
% %configuration 10 (Single phase B)
% R10 = 1.3292/1.60934; % 1.60934 converts ohms per mile into ohms per kilometers
% L10 = 1.3475/(2*pi*f)/1.60934; 
% C10 = 4.5193*1e-6/(2*pi*f)/1.60934;
% 
% %configuration 11 (Single phase C)
% R11 = 1.3292/1.60934; % 1.60934 converts ohms per mile into ohms per kilometers
% L11 = 1.3475/(2*pi*f)/1.60934; 
% C11 = 4.5193*1e-6/(2*pi*f)/1.60934;
% 
% %configuration 12
% R12 = [1.5209 0.5198 0.4924; 0.5198 1.5329 0.5198; 0.4924 0.5198 1.5209]/1.60934; % 1.60934 converts ohms per mile into ohms per kilometers
% L12 = [0.7521 0.2775 0.2157; 0.2775 0.7162 0.2775; 0.2157 0.2775 0.7521]/(2*pi*f)/1.60934;
% CAg12 = (67.2242)*1e-6/(2*pi*f)/1.60934;
% CBg12 = (67.2242)*1e-6/(2*pi*f)/1.60934;
% CCg12 = (67.2242)*1e-6/(2*pi*f)/1.60934;

%% TAP REGULATORS


%%% NOT USED in the model
TR1_KV=1.04;%coefficient to modify the low Voltage on the #1 TAP (V_low = 4.16e3*TR1_KV)****1.04 is the normal****
TR2_KV=1; % #2 TAP
TR3_KV=1; % #3 TAP
TR4_KV=1.048; % #4 TAP

%% Line data & spot load power data
% Dec 2022: remove dependency on Excel file
data_ieee123={
'Node'	'Load'	'Ph-1'	'Ph-1'	'Ph-2'	'Ph-2'	'Ph-3'	'Ph-4';
'*'	'Model'	'kW'	'kVAr'	'kW'	'kVAr'	'kW'	'kVAr';
1	'Y-PQ'	40	20	0	0	0	0;
2	'Y-PQ'	0	0	20	10	0	0;
4	'Y-PR'	0	0	0	0	40	20;
5	'Y-I'	0	0	0	0	20	10;
6	'Y-Z'	0	0	0	0	40	20;
7	'Y-PQ'	20	10	0	0	0	0;
9	'Y-PQ'	40	20	0	0	0	0;
10	'Y-I'	20	10	0	0	0	0;
11	'Y-Z'	40	20	0	0	0	0;
12	'Y-PQ'	0	0	20	10	0	0;
16	'Y-PQ'	0	0	0	0	40	20;
17	'Y-PQ'	0	0	0	0	20	10;
19	'Y-PQ'	40	20	0	0	0	0;
20	'Y-I'	40	20	0	0	0	0;
22	'Y-Z'	0	0	40	20	0	0;
24	'Y-PQ'	0	0	0	0	40	20;
28	'Y-I'	40	20	0	0	0	0;
29	'Y-Z'	40	20	0	0	0	0;
30	'Y-PQ'	0	0	0	0	40	20;
31	'Y-PQ'	0	0	0	0	20	10;
32	'Y-PQ'	0	0	0	0	20	10;
33	'Y-I'	40	20	0	0	0	0;
34	'Y-Z'	0	0	0	0	40	20;
35	'D-PQ'	40	20	0	0	0	0;
37	'Y-Z'	40	20	0	0	0	0;
38	'Y-I'	0	0	20	10	0	0;
39	'Y-PQ'	0	0	20	10	0	0;
41	'Y-PQ'	0	0	0	0	20	10;
42	'Y-PQ'	20	10	0	0	0	0;
43	'Y-Z'	0	0	40	20	0	0;
45	'Y-I'	20	10	0	0	0	0;
46	'Y-PQ'	20	10	0	0	0	0;
47	'Y-I'	35	25	35	25	35	25;
48	'Y-Z'	70	50	70	50	70	50;
49	'Y-PQ'	35	25	70	50	35	20;
50	'Y-PQ'	0	0	0	0	40	20;
51	'Y-PQ'	20	10	0	0	0	0;
52	'Y-PQ'	40	20	0	0	0	0;
53	'Y-PQ'	40	20	0	0	0	0;
55	'Y-Z'	20	10	0	0	0	0;
56	'Y-PQ'	0	0	20	10	0	0;
58	'Y-I'	0	0	20	10	0	0;
59	'Y-PQ'	0	0	20	10	0	0;
60	'Y-PQ'	20	10	0	0	0	0;
62	'Y-Z'	0	0	0	0	40	20;
63	'Y-PQ'	40	20	0	0	0	0;
64	'Y-I'	0	0	75	35	0	0;
65	'D-Z'	35	25	35	25	70	50;
66	'Y-PQ'	0	0	0	0	75	35;
68	'Y-PQ'	20	10	0	0	0	0;
69	'Y-PQ'	40	20	0	0	0	0;
70	'Y-PQ'	20	10	0	0	0	0;
71	'Y-PQ'	40	20	0	0	0	0;
73	'Y-PQ'	0	0	0	0	40	20;
74	'Y-Z'	0	0	0	0	40	20;
75	'Y-PQ'	0	0	0	0	40	20;
76	'D-I'	105	80	70	50	70	50;
77	'Y-PQ'	0	0	40	20	0	0;
79	'Y-Z'	40	20	0	0	0	0;
80	'Y-PQ'	0	0	40	20	0	0;
82	'Y-PQ'	40	20	0	0	0	0;
83	'Y-PQ'	0	0	0	0	20	10;
84	'Y-PQ'	0	0	0	0	20	10;
85	'Y-PQ'	0	0	0	0	40	20;
86	'Y-PQ'	0	0	20	10	0	0;
87	'Y-PQ'	0	0	40	20	0	0;
88	'Y-PQ'	40	20	0	0	0	0;
90	'Y-I'	0	0	40	20	0	0;
92	'Y-PQ'	0	0	0	0	40	20;
94	'Y-PQ'	40	20	0	0	0	0;
95	'Y-PQ'	0	0	20	10	0	0;
96	'Y-PQ'	0	0	20	10	0	0;
98	'Y-PQ'	40	20	0	0	0	0;
99	'Y-PQ'	0	0	40	20	0	0;
100	'Y-Z'	0	0	0	0	40	20;
102	'Y-PQ'	0	0	0	0	20	10;
103	'Y-PQ'	0	0	0	0	40	20;
104	'Y-PQ'	0	0	0	0	40	20;
106	'Y-PQ'	0	0	40	20	0	0;
107	'Y-PQ'	0	0	40	20	0	0;
109	'Y-PQ'	40	20	0	0	0	0;
111	'Y-PQ'	20	10	0	0	0	0;
112	'Y-I'	20	10	0	0	0	0;
113	'Y-Z'	40	20	0	0	0	0;
114	'Y-PQ'	20	10	0	0	0	0;
'Total'	1420 775 915 515 11 55 630;};

len=size(data_ieee123,1);  % first 2 line and last line are not data
for i=3:len-1
    eval(['Pa' num2str(i-2) '=data_ieee123{' num2str(i) ',3}*1000;']);   
    eval(['Qa' num2str(i-2) '=data_ieee123{' num2str(i) ',4}*1000;']); 
    eval(['Pb' num2str(i-2) '=data_ieee123{' num2str(i) ',5}*1000;']); 
    eval(['Qb' num2str(i-2) '=data_ieee123{' num2str(i) ',6}*1000;']); 
    eval(['Pc' num2str(i-2) '=data_ieee123{' num2str(i) ',7}*1000;']); 
    eval(['Qc' num2str(i-2) '=data_ieee123{' num2str(i) ',8}*1000;']); 
    %note that the parameter 'loadcoeff' is included in the block mask
    %parameter and not included here.
end

% % data_line = xlsread('line_data_ieee123node.xls','C4:C121'); %length of all lines (ft.)
% data_load = xlsread('spot_load_data_ieee123node.xls','C5:H89'); % active and reactive power for each spot load
% 
% 
% % i = 1;
% % while i <= length (data_line') % length of every line
% %     eval(['Len' num2str(i) '=data_line(' num2str(i) ',1)*0.0003048;']); %length (km)
% %     i=i+1;
% % end
% 
% i = 1;
% while i <= length (data_load') % power info for each spot load
%     eval(['Pa' num2str(i) '=data_load(' num2str(i) ',1)*1000*loadcoeff;']);
%     eval(['Qa' num2str(i) '=data_load(' num2str(i) ',2)*1000*loadcoeff;']);
%     eval(['Pb' num2str(i) '=data_load(' num2str(i) ',3)*1000*loadcoeff;']);
%     eval(['Qb' num2str(i) '=data_load(' num2str(i) ',4)*1000*loadcoeff;']);
%     eval(['Pc' num2str(i) '=data_load(' num2str(i) ',5)*1000*loadcoeff;']);
%     eval(['Qc' num2str(i) '=data_load(' num2str(i) ',6)*1000*loadcoeff;']);
%     i=i+1;
% end