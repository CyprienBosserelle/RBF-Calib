function [subset, ncenters] = algoritmo_MaxDiss_MaxMinSimplificado_SinUmbral (seed, num, train_n, scalar, directional)  

%Centroids initilization
subset=[];
subset=[train_n(seed,:)];

%Elimination of the selected centroid in the training data 
elementos=1:length(train_n);
dum=find(elementos~=seed);
train_n2=train_n(dum,:);

%The process is repeted until the desired number of centroids are found. 
ncenters=1;

qerr=0;
while ncenters<num 

    [m1,n1]=size(train_n2);
    m=ones(m1,1);
    [m2,n2]=size(subset);
    if m2==1 
        xx2=subset(m,:);
        dultima=distancia_normalizada(train_n2, xx2, scalar, directional);
    else
        xx=subset(end,:);
        xx2=xx(m,:);
        danterior=distancia_normalizada(train_n2,xx2, scalar, directional);
        [dultima,ii]=min([danterior dultima],[],2); 
    end
        
    [qerr,bmu]=max(dultima);
    
    if isnan(qerr)==0
        subset=[subset; train_n2(bmu,:)];
        elementos=1:length(train_n2);
        dum=find(elementos~=bmu(1));
        train_n2=train_n2(dum,:);
        dultima=dultima(dum);
    end
        
    [ncenters,dim]=size(subset);
  
        
end