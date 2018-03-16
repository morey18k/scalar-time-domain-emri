filename='bhpert.csv';
data=csvread(filename);
figure
plot(data(:,2),data(:,3),'b'); hold on;
plot(data(:,2),data(:,4),'r')