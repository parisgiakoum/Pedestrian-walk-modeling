%% Checking Martina's work
%% Time
Time = X(:,1,1);
for i=2:215
    Time = [Time;X(:,1,i)];
end
Time = Time(~isnan(Time));

figure;
histogram (Time, 'BinWidth', 0.05, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );

figure;
normplot(Time)
%% Force
Force = X(:,2,1);
for i=2:215
    Force = [Force;X(:,2,i)];
end
Force = Force(~isnan(Force));

figure;
histogram (Force, 'normalization' , 'pdf' );

figure;
normplot(Force)

%% Length
L = X(:,3,1);
for i=2:215
    L = [L;X(:,3,i)];
end
L = L(~isnan(L));

figure;
histogram (L, 'BinWidth', 0.05, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );

figure;
normplot(L)

%% Angle
A= X(:,4,1);
for i=2:215
    A = [A;X(:,4,i)];
end
A=A(~isnan(A));

figure;
histogram (A, 'normalization' , 'pdf');

figure;
normplot(A)
