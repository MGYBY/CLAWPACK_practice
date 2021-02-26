dx = 42.0/2100;
U = 1.03774;
H = 0.00798;
W = 0.40;

m1 = importdata('./centerline_front',' ');
m2 = importdata('./centerline_rear',' ');
time = importdata('./centerline_time',' ');

delta_fw = sum((0.50*9.81*(10^(3))*dx*(m1.^2-m2.^2)), 2); % take sum by cols
normal_fw = sum((0.50*(10^(3))*W*(U^2)*H), 2); % take sum by cols

cd = delta_fw./normal_fw;
cd_time_series = [time, cd];

csvwrite ('./cd_time_series.csv', cd_time_series)

