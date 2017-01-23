% Calculate the MDA for calibration using 6 parameters
% gamma n beta cf fw nuh/smag nuhfac


%total number of simulation

totaln=1000000;
%Size of the selected subset 
quantity = 300;


% randomly sample within the range below
nvar=7;
gammarange=[0.4 0.8];
nrange=[6 11];
betarange=[0.09 0.25];
cfrange=[0.01 0.1];
fwrange=[0.1 0.4];
smagrange=[0.0001 0.5];
nuhfacrange=[0.0005 1.0];

Varst=[gammarange(1) nrange(1) betarange(1) cfrange(1) fwrange(1) smagrange(1) nuhfacrange(1)];
Vardev=[gammarange(end)-gammarange(1) nrange(end)-nrange(1) betarange(end)-betarange(1) cfrange(end)-cfrange(1) fwrange(end)-fwrange(1) smagrange(end)-smagrange(1) nuhfacrange(end)-nuhfacrange(1)];

[VarstM dummy]=meshgrid(Varst,1:totaln);
[VardevM dummy]=meshgrid(Vardev,1:totaln);

%We need total x nvar uniform random samples
datos=rand(totaln,nvar).*VardevM+VarstM;
% Save files for the RBF?
save('Datos_Allmod.dat','datos','-ascii');

[N,dim]=size(datos);

%Definition of the scalar and directional parameters
scalar=[1 2 3 4 5 6 7];
directional=[];

%Data Parameter
parameters={'gamma'; 'n'; 'beta'; 'cf'; 'fw'; 'smag'; 'nuhfac'};

%Data normalization
minimos=zeros(length(scalar),1);
maximos=zeros(length(scalar),1);
for i=1:length(scalar)
    minimos(i)=min(datos(:,scalar(i)));
    maximos(i)=max(datos(:,scalar(i)));
end

datos_n=zeros(N,dim);
for i=1:length(scalar)
    datos_n(:,scalar(i))=(datos(:,scalar(i))-minimos(i))./(maximos(i)-minimos(i));
end



seed=find(max(datos(:,1))==datos(:,1));
[subset, ncenters] = algoritmo_MaxDiss_MaxMinSimplificado_SinUmbral(seed(1), quantity, datos_n, scalar, directional);

% bmus calculation of all the data
bmus=zeros(N,1);
m=ones(quantity,1);
for i=1:N
    retro=datos_n(i,:);
    retro2=retro(m,:);
    [qerr,bmu]=dist_normalizada(subset,retro2,scalar, directional);
    bmus(i)=bmu;
end
    
%Desnormalization of the results
final=zeros(quantity,dim);
for i=1:length(scalar)
    final(:,scalar(i))=subset(:,scalar(i))*(maximos(i)-minimos(i))+minimos(i);
end

save('MDA_final.dat','final','-ascii');
save('MDA_bmus.dat','bmus','-ascii'); 


for n=2:length(scalar)
    figure
    plot(datos(:,1),datos(:,n),'.g')
    hold on
    plot(final(:,1),final(:,n),'.r');
end

%Frequency of presentation of the centroids based on all the retroanalisis data
acierto=zeros(quantity,1);
for i=1:quantity
    acierto(i)=sum(bmus==i);
end

acierto_cien(1:quantity,1)=acierto/length(bmus)*100;


%Find the selected subset in all the data to identify the dates and the real cases.
pos=zeros(quantity,1);
date=zeros(quantity,4);
for i=1:length(subset)
    mm=ones(length(datos_n),1);
    temp1=subset(i,:);
    temp2=temp1(mm,:);
    temp3=datos_n-temp2;
    temp4=sum(temp3,2);
    temp5=find(temp4==0);  %Position in the reanalysis data
    pos(i)=temp5(1);
  %  date(i,:)=time(pos(i),1:4); %(time = vector time of the reanalysis serie)
end
% Save file
save('MDA_pos.dat','pos','-ascii');




