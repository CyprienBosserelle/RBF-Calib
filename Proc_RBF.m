

                                  %%%% RBF Reconstruction %%%%
 
 
 %Reconstruction of the complete serie of deep water data in a certain position through RBF functions based on 
 %the value of these parameters (scalar and directional) in the propagated cases.
 
%This script calls the functions

%LecturaresultsPropagaciones_Interpolacion_4parameters ( It is necessary to change the file name from where it read the results)
%InterpolacionRBF_Parametros_old ( It is necessary to change the sigma values to be the same as the ones defined in this script)
  %rbfcreate_modificado_matrizA
  %rbfcreate_modificado
  %rbfinterp_modificado
  
%It is necesary to check the function norma !! to define the position of the scalar and directional parameters.
%Norma take into account the distance in the circle of the directional parameters.

clear all
close all
clc

% Select MDA centers
centers=load('MDA_final.dat');   %They are in the real space without any transformation
%reduce to completed simulation so far
centers=centers(1:300,:);



centers_n=zeros(size(centers));
%Normalization of the selected cases with MaxDiss
for v=1:7
    maxH=max(centers(:,v));  minH=min(centers(:,v));
    centers_n(:,v)=(centers(:,v)-minH)./(maxH-minH);
end


% Define the spacing and scaling for the interoloation
totaln=100000;

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

%Normalization of the reanalysis data
datos_n=zeros(size(datos));
%Normalization
for v=1:7
    maxH=max(centers(:,v));  minH=min(centers(:,v));
    datos_n(:,v)=(datos(:,v)-minH)./(maxH-minH);
end


%Sigma values
sigmin=0.10;
sigmax=0.70;

%Quantity of propagated cases
ncases=300;
Propagations=zeros(ncases,3);  %zserror higerror hserror 

fid=fopen('XBoutput.txt');
C=textscan(fid,'%f %f %f %f %f %f', 'headerlines',1);
fclose(fid);

Propagations(:,1)=C{4};
Propagations(:,2)=C{5};
Propagations(:,3)=C{6};

%Parameters to reconstruct through RBF
parameters=3; %zserror higerror hserror 
ndireccion=0;  %Quantity of directional parameters (Dir)
cdireccion=0; %Column of the result file which is a directional parameter
    
[optimal_sigma, results] = InterpolationRBF_Parameters (parameters, ndireccion, cdireccion, centers_n, datos_n, Propagations);

[a b]=min(sum(abs(results),2));

fop=fopen('RBF_results.txt','w');
for n=1:length(results)
    fprintf(fop,'%f %f %f %f %f %f %f %f %f %f\n',datos(n,1),datos(n,2),datos(n,3),datos(n,4),datos(n,5),datos(n,6),datos(n,7),results(n,1),results(n,2),results(n,3));
end
fclose(fop);
% Find the 20 best results
best= (abs(results(:,1))+ abs(results(:,2))+ abs(results(:,3))/2)/3; % half the weight is put on Hs
[dummy,indx]=sort(best);



%indx=find(abs(results(:,1))<0.03 & abs(results(:,2))<0.03 & abs(results(:,3))<0.05);

fop=fopen('RBF_best.txt','w');
for n=1:20%length(indx)
    fprintf(fop,'%f %f %f %f %f %f %f %f %f %f\n',datos(indx(n),1),datos(indx(n),2),datos(indx(n),3),datos(indx(n),4),datos(indx(n),5),datos(indx(n),6),datos(indx(n),7),results(indx(n),1),results(indx(n),2),results(indx(n),3));
end
fclose(fop);

