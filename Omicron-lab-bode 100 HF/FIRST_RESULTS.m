%Program to analyse data collected from Omicron Bode 100
%Author: Refiloe Leanya

  %*****************************Common Mode Data**************************
  %***********************************************************************
  %L1
Comm = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\healthy ' ...
    'motor impedance\csv\L1toGND.csv']);
Comm(:,:,2) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\Fault ' ...
    'maybe 5%\csv\L1toGND.csv']);
Comm(:,:,3) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\10' ...
    '%\csv\L1toGND.csv']);
% L2

Comm(:,:,4) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\healthy ' ...
    'motor impedance\csv\L2toGND.csv']);
Comm(:,:,5) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\Fault ' ...
    'maybe 5%\csv\L2toGND.csv']);
Comm(:,:,6) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\10' ...
    '%\csv\L2toGND.csv']);

%L3

Comm(:,:,7) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\healthy ' ...
    'motor impedance\csv\L3toGND.csv']);
Comm(:,:,8) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\Fault ' ...
    'maybe 5%\csv\L3toGND.csv']);
Comm(:,:,9) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\10' ...
    '%\csv\L3toGND.csv']);

 %*****************************Differential Mode Data**************************
  %***********************************************************************
  %L1
Diff = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\healthy ' ...
    'motor impedance\csv\L1.csv']);
Diff(:,:,2) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\Fault ' ...
    'maybe 5%\csv\L1.csv']);
Diff(:,:,3) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\10' ...
    '%\csv\L1.csv']);
% L2

Diff(:,:,4) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\healthy ' ...
    'motor impedance\csv\L2.csv']);
Diff(:,:,5) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\Fault ' ...
    'maybe 5%\csv\L2.csv']);
Diff(:,:,6) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\10' ...
    '%\csv\L2.csv']);

%L3

Diff(:,:,7) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\healthy ' ...
    'motor impedance\csv\L3.csv']);
Diff(:,:,8) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\Fault ' ...
    'maybe 5%\csv\L3.csv']);
Diff(:,:,9) = readmatrix(['C:\Users\R.Leanya\Documents\UCT\Fourth year\SECOND ' ...
    'SEMESTER\EEE4022S\OMICRON-LAB-BODE 100\First results\10' ...
    '%\csv\L3.csv']);

%%%%%%%%%%%%%%%%%%Calculate Model Parameters
%%%%%%%%%%COMMON MODE
%%%CALCULATING Healthy
CHF = zeros(9,1);
Cg = zeros(9,1);
Ct = zeros(9,1);
Zmin1_fres1 = zeros(9,1);
n=ones(9,1);
k = zeros(9,1);
l = zeros(9,1);
rate = 1;
rate1 = 0;
s = 2*ones(9,1);
fres1 = zeros(3,1);
for ref = [1,2,3,4,5,6,7,8,9]
    while n(ref)<3
        Healthy = Comm(:,:,ref);
        %fprintf('Cg1: %E\n', length(HealthyComm));
        if (Healthy(s(ref),4)<Healthy(s(ref)-1,4) &&  n(ref)==1)
            %fprintf('f: %E\n', Healthy(s(ref),1))
            if Healthy(s(ref),1)<20000
                Cg(ref) = Cg(ref)+(1/(6*2*pi*Healthy(s(ref),4).*Healthy(s(ref),1)));
                %rate = rate-0.0002;
                
                k(ref) = k(ref)+1;
            end
            s(ref) = s(ref)+1;
        else
            %disp('Hello, World!');
            if n(ref)==1
                Zmin1_fres1(ref) = Healthy(s(ref)-1,4);
                fres1(ref) = Healthy(s(ref)-1,1);
                Cg(ref) =Cg(ref)./k(ref);
            end
            n(ref)=2;
         
            s(ref) = s(ref)+1;
            if Healthy(s(ref),1)>600000 && Healthy(s(ref),1)<650000
                %fprintf('f: %E\n', Ct(2))
                Ct(ref) = Ct(ref)+(Cg(ref)-3*Cg(ref).^2*2*pi.*Healthy(s(ref),1).*Healthy(s(ref),4))...
                    ./((6.*Cg(ref)*2*pi.*Healthy(s(ref),1).*Healthy(s(ref),4)-1));
                %rate1 = rate1+0.0013;
                l(ref) = l(ref)+1;
                s(ref) = s(ref)+1;
                
            elseif s(ref)>350
            n(ref) =3;
            end
        end
    end
end




Zmax = zeros(9,1);
Zmin2 = zeros(9,1);
fpmin = zeros(9,1);
nn=ones(9,1);
kk = zeros(9,1);
ll = zeros(9,1);
ss = 2*ones(9,1);
fp1 = zeros(9,1);
for ref = [1,2,3,4,5,6,7,8,9]
    while nn(ref)<3
        Healthy = Diff(:,:,ref);
       
        if (Healthy(ss(ref),4)>Healthy(ss(ref)-1,4) &&  nn(ref)==1)            
            ss(ref) = ss(ref)+1;
            kk(ref) = kk(ref)+1;
        else
            if nn(ref)==1
                Zmax(ref) = Healthy(ss(ref)-1,4);
                fp1(ref) = Healthy(ss(ref)-1,1);
                ss(ref) = ss(ref)+1;
                nn(ref) = 2;
            else
                ss(ref) = ss(ref)+1;
               
            end
            if (Healthy(ss(ref),4)>Healthy(ss(ref)-1,4) &&  nn(ref)==2 && ss(ref)>200) 
                
                Zmin2(ref) = Healthy(ss(ref)-1,4);
                fpmin(ref) = Healthy(ss(ref)-1,1);
                nn(ref)=3;
                 
            end
  
        end
    end
end








%CHF = CHF./l;

Ct = Ct./l;
%Lcomm = abs(1./(12*pi^2.*Ct.*fres1.^2));
Ld = 1./((2*pi.*fp1).^2.*(Cg+Ct));
Lzu = 1./(4*pi^2.*(Cg/2+Ct).*fpmin.^2*100);
Rg =Zmin2;
Re = 3*Zmax;%abs((2*pi*3.*Zmax.*fp1.*Ld)./(2*pi.*fp1.*Ld-3.*Zmax.*(Ld.*(Ct+Cg/2).*(2*pi.*fp1).^2)-3.*Zmax));
f = Comm(:,1,1);
s = 2*pi*f*1j;
Zwg = zeros(401,9);
Zwn = zeros(401,9);
for i = [1,2,3,4,5,6,7,8,9]
    Zwg(:,i) =((s.^2*Ld(i)*(Ct(i)+Cg(i))+Ld(i)/(Re(i)).*s+1))./...
        ...
    (6.*s*Cg(i).*(s.^2.*Ld(i).*(Ct(i)+Cg(i)/2)+s.*Ld(i)./Re(i)+1));


    Zwn(:,i) =((Ld(i).*s))./...
        ...
    (3.*(s.^2.*Ld(i).*(Ct(i)+Cg(i))+s.*Ld(i)./Re(i)+1));
end
fprintf('Cg: %E\n', Cg);
fprintf('Ct: %E\n', Ct);
fprintf('Lzu: %E\n', Lzu);
fprintf('Ld: %E\n', Ld);

%%%%%%%%%%DIFFERENTIAL MODE


%CALCULATING THE MEAN SQUARE DIVIATION
%This section contains code to calculate the squared deviation of 5%
%and 10% fault

%Calculate room mean square diviation
sum_squared_deviation = zeros(1,1,9);
squared_dev = zeros(401,1,9);

for i = 1:401
    deviationfive = Comm(i,4,1) - Comm(i,4,:);
    %fprintf('Deviation: %.2f\n', deviationfive);
    squared_deviation = deviationfive.^2;
    squared_dev(i,1,:) = squared_deviation(1,1,:);

    
    sum_squared_deviation = sum_squared_deviation+squared_deviation;
end
% Calculate the Root Mean Square Deviation
rms_deviation = sqrt(sum_squared_deviation / 401);

% Display the result
fprintf('Root Mean Square Deviation: %.2f\n', rms_deviation);



% Display the moving average
%disp(moving_average);

loglog((Comm(:,1,1)),(Comm(:,4,1)))
grid on, xlabel('Frequency [Hz]', fontName = 'Arial',fontSize = 11)
ylabel('Magnitude [Ω]', fontName = 'Arial',fontSize = 11)
hold on
loglog((Diff(:,1,1)),abs(Zwg(:,1)))
hold on
loglog(Diff(:,1,1),Diff(:,4,1))
hold on
%loglog(Diff(:,1,2),Diff(:,4,2))
%hold on
%loglog((Comm(:,1,2)),(Comm(:,4,2)))
%hold on
%loglog((Comm(:,1,2)),squared_dev(:,1,2))
%hold on
%loglog((Comm(:,1,2)),squared_dev(:,1,3))
loglog(Diff(:,1,1),abs(Zwn(:,1)))
legend('Common Measured','Common simulated','Diff Measured','Diff simulated')
hold off
