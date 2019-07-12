data = (webread('http://192.168.1.3/now'));
A = data{2,1};
B = string(A);
C = strsplit(B,',');
keysstr = "Time,Temp S0,Humidity S0,Pressure S0,X Accel S0,Y Accel S0,Z Accel S0,X Gyro S0,Y Gyro S0,Z Gyro S0,X Mag S0,Y Mag S0,Z Mag S0,Temp S1,Humidity S1,Pressure S1,X Accel S1,Y Accel S1,Z Accel S1,X Gyro S1,Y Gyro S1,Z Gyro S1,X Mag S1,Y Mag S1,Z Mag S1,Temp S2,Humidity S2,Pressure S2,X Accel S2,Y Accel S2,Z Accel S2,X Gyro S2,Y Gyro S2,Z Gyro S2,X Mag S2,Y Mag S2,Z Mag S2,Temp S3,Humidity S3,Pressure S3,X Accel S3,Y Accel S3,Z Accel S3,X Gyro S3,Y Gyro S3,Z Gyro S3,X Mag S3,Y Mag S3,Z Mag S3,Temp S4,Humidity S4,Pressure S4,X Accel S4,Y Accel S4,Z Accel S4,X Gyro S4,Y Gyro S4,Z Gyro S4,X Mag S4,Y Mag S4,Z Mag S4,Temp S5,Humidity S5,Pressure S5,X Accel S5,Y Accel S5,Z Accel S5,X Gyro S5,Y Gyro S5,Z Gyro S5,X Mag S5,Y Mag S5,Z Mag S5,Temp S6,Humidity S6,Pressure S6,X Accel S6,Y Accel S6,Z Accel S6,X Gyro S6,Y Gyro S6,Z Gyro S6,X Mag S6,Y Mag S6,Z Mag S6,Temp S7,Humidity S7,Pressure S7,X Accel S7,Y Accel S7,Z Accel S7,X Gyro S7,Y Gyro S7,Z Gyro S7,X Mag S7,Y Mag S7,Z Mag S7";
keys = strsplit(keysstr,',');
dict = containers.Map;
n = 1;
while n <= length(keys)
    dict(keys{n}) = C{n};
    n = n + 1;
end

JJ = 1;
while JJ <= length(keys)
    string(keys{JJ}) + ' = ' + string(dict(keys{JJ}))
    JJ = JJ + 1;
end


%class(C)
% v = table(data, 'VariableNames', {"Time","Temp S0","Humidity S0","Pressure S0","X Accel S0","Y Accel S0","Z Accel S0","X Gyro S0","Y Gyro S0","Z Gyro S0","X Mag S0","Y Mag S0","Z Mag S0","Temp S1","Humidity S1","Pressure S1","X Accel S1","Y Accel S1","Z Accel S1","X Gyro S1","Y Gyro S1","Z Gyro S1","X Mag S1","Y Mag S1","Z Mag S1","Temp S2","Humidity S2","Pressure S2","X Accel S2","Y Accel S2","Z Accel S2","X Gyro S2","Y Gyro S2","Z Gyro S2","X Mag S2","Y Mag S2","Z Mag S2","Temp S3","Humidity S3","Pressure S3","X Accel S3","Y Accel S3","Z Accel S3","X Gyro S3","Y Gyro S3","Z Gyro S3","X Mag S3","Y Mag S3","Z Mag S3","Temp S4","Humidity S4","Pressure S4","X Accel S4","Y Accel S4","Z Accel S4","X Gyro S4","Y Gyro S4","Z Gyro S4","X Mag S4","Y Mag S4","Z Mag S4","Temp S5","Humidity S5","Pressure S5","X Accel S5","Y Accel S5","Z Accel S5","X Gyro S5","Y Gyro S5","Z Gyro S5","X Mag S5","Y Mag S5","Z Mag S5","Temp S6","Humidity S6","Pressure S6","X Accel S6","Y Accel S6","Z Accel S6","X Gyro S6","Y Gyro S6","Z Gyro S6","X Mag S6","Y Mag S6","Z Mag S6","Temp S7","Humidity S7","Pressure S7","X Accel S7","Y Accel S7","Z Accel S7","X Gyro S7","Y Gyro S7","Z Gyro S7","X Mag S7","Y Mag S7","Z Mag S7"})