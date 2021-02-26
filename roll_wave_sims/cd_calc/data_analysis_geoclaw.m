dx = 0.40/5.0;
U = 1.03774;
H = 0.00798;
W = 0.40;

% 10 gauges
M1 = importdata("gauge00001.txt", ' ',3);
M2 = importdata("gauge00002.txt", ' ',3);
M3 = importdata("gauge00003.txt", ' ',3);
M4 = importdata("gauge00004.txt", ' ',3);
M5 = importdata("gauge00005.txt", ' ',3);
M6 = importdata("gauge00006.txt", ' ',3);
M7 = importdata("gauge00007.txt", ' ',3);
M8 = importdata("gauge00008.txt", ' ',3);
M9 = importdata("gauge00009.txt", ' ',3);
M10 = importdata("gauge00010.txt", ' ',3);

m1 = M1.data;
m2 = M2.data;
m3 = M3.data;
m4 = M4.data;
m5 = M5.data;
m6 = M6.data;
m7 = M7.data;
m8 = M8.data;
m9 = M9.data;
m10 = M10.data;

time = m1(:,2)

depth1 = m1(:,3);
depth2 = m2(:,3);
depth3 = m3(:,3);
depth4 = m4(:,3);
depth5 = m5(:,3);
depth6 = m6(:,3);
depth7 = m7(:,3);
depth8 = m8(:,3);
depth9 = m9(:,3);
depth10 = m10(:,3);

fw_1 = 0.50*9.81*(10^(3))*dx*(depth1.^2+depth2.^2+depth3.^2+depth4.^2+depth5.^2); 
fw_2 = 0.50*9.81*(10^(3))*dx*(depth6.^2+depth7.^2+depth8.^2+depth9.^2+depth10.^2);
delta_fw = fw_1-fw_2;
normal_fw = (0.50*(10^(3))*W*(U^2)*H); % take sum by cols

coeff_cd = delta_fw./normal_fw;
cd_time_series = [time, coeff_cd];

csvwrite ('./cd_time_series.csv', cd_time_series)

