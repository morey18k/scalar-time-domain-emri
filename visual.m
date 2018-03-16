for i=25
    filename=sprintf('dn%02d,01,01mode.csv', i);
    disp(filename)
    data=csvread(filename);
    figure
    plot(data(:,1),data(:,4),'b'); hold on;
    plot(data(:,1),data(:,5),'r')
end
