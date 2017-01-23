function [optimal_sigma, results] = InterpolationRBF_Parameters (parameters, ndireccion, cdireccion, subset_n, PCs_n, Propagations)

%Sigma values
sigmin=0.1;
sigmax=0.7;


%Calculation of the optimal sigma for each parameter
optimal_sigma=zeros(parameters+ndireccion,1);
results=zeros(length(PCs_n(:,1)),parameters);
par=1;


for j=1:parameters
    time_0 = tic;
    if ismember(j,cdireccion) % Direction variable
        
        % X and Y direction component
        m=length(Propagations(:,1));
        DirX=zeros(m,1);
        DirY=zeros(m,1);
        temp=pi/2-Propagations(:,j)*pi/180;
        dd=find(temp<-pi);
        temp(dd)=2*pi+temp(dd);
        DirX=cos(temp);  %angular component
        DirY=sin(temp);
        optimal_sigma(par) = evaluarmin(subset_n', DirX', sigmin, sigmax);   %DirX
        optimal_sigma(par+1) = evaluarmin(subset_n', DirY', sigmin, sigmax);  %DirY
        time_sigma=toc(time_0)/60;
        
        % DirX serie reconstruction
        coeff_DirX = rbfcreate_modificado(subset_n',DirX','RBFFunction','gaussian','RBFConstant',optimal_sigma(par));
        DirXp = rbfinterp_modificado(PCs_n',coeff_DirX);
        
        % DirY serie reconstruction
        coeff_DirY = rbfcreate_modificado(subset_n',DirY','RBFFunction','gaussian','RBFConstant',optimal_sigma(par+1));
        DirYp = rbfinterp_modificado(PCs_n',coeff_DirY);
        
        %Dir
        Dirp = atan2(DirYp,DirXp)*180/pi; % angle based on matlab criteria
        Dirp = 90-Dirp;
        qq=find(Dirp<0);  Dirp(qq)=Dirp(qq)+360;
        results(:,j)=Dirp;
        time_result=toc(time_0)/60-time_sigma;
        
        par=par+2;
    else
        optimal_sigma(par) = evaluarmin(subset_n', Propagations(:,j)', sigmin, sigmax);
        time_sigma=toc(time_0)/60;
        coeff = rbfcreate_modificado(subset_n',Propagations(:,j)','RBFFunction','gaussian','RBFConstant',optimal_sigma(par));
        results(:,j) = rbfinterp_modificado(PCs_n',coeff);
        par=par+1;
        time_result=toc(time_0)/60-time_sigma;
    end
    disp(['Parameter' num2str(j) ', ( sigma time:' num2str(time_sigma) ',date time:' num2str(time_result) ')'])
    clear time*
end
