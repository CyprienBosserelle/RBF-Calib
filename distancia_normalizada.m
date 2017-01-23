function Dist = distancia_normalizada( M, D, escalar, direccional)  

     [m,n]=size(M);   
     dif=zeros(m,n);
     for i=1:length(escalar)
         dif(:,escalar(i))=D(:,escalar(i))-M(:,escalar(i));
     end
     for i=1:length(direccional)
         dif(:,direccional(i))=min(abs(D(:,direccional(i))-M(:,direccional(i))),2*pi-abs(D(:,direccional(i))-M(:,direccional(i))))/pi;
     end
     Dist=sum(dif.^2,2);
    
        
      