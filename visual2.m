filename='sf.csv';
data=csvread(filename);
shape=size(data);
shape(1);
real(1:70000)=0.04887074683528891054717664;
imag(1:70000)=0.000429346057070435914264781;
figure
plot(data(:,1),real.','b'); hold on;
plot(data(:,1),imag.','r'); hold on;
plot(data(:,1),data(:,2),'b'); hold on;
plot(data(:,1),data(:,3),'r')

