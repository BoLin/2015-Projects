%GENETIC ALGORITHM TO SOLVE RATSTRIGIN FUNCTION
clc; clear all; close all
tic
disp('Program for GA is running....')
hh=waitbar(0,'Please wait the program for GA is running.....');
%% function and side constraints definition

f=@(x) 10*2+(x(1)^2-10*cos(2*pi*x(1)))+(x(2)^2-10*cos(2*pi*x(2)));
%here f(x) and f(xp) gives the same  value
%-5<=x1<=5; -5<=x2<=5    These are the bounds which can be converted to
% 0<=x1'<=10 and  0<=x2'<=10; where x1'=x1+5;   x2'=x2+5;
g1=@(xp) -xp(1); g2=@(xp) xp(1)-10; g3=@(xp) -xp(2); g4=@(xp) xp(2)-10;
r=10^6;
PHI=@(xp) f(xp-5)+r*(max(0,g1(xp))^2+max(0,g2(xp))^2+max(0,g3(xp))^2+max(0,g4(xp))^2);
%x=xp-5
%% START  //CREATE AND INITIAL SWARM. HERE CALCULATION IS DONE USING XP....
%%AND LATER CHANGED TO X ORIGINAL VARIABLE
x1=linspace(-5,5,100); x2=linspace(-5,5,100);
[X1,X2]=meshgrid(x1,x2);
objecfun=10*2+(X1^2-10*cos(2*pi*X1))+(X2^2-10*cos(2*pi*X2));
contour(X1,X2,objecfun,5,'r','linewidth',1)
for globitr=1:3

%% create and initial population of 20 populations or chromosomes
num=20; Nmaxiter=100; precision=.001;
rng('shuffle'); var_p=rand(num,2)*10; %10 is the maximum of DV
var_orig=var_p-5;
syms m
bits_num=round(double(solve(2^m==(10-0)/(precision)+1,m)));  %number of bits required for accuracy of 0.01

%create Pseudo-objective function for exterior penalty approach for side
%constraints and then find the fitness
totitern=0; N1=0; N2=0; N3=0; repeat_convrg=15; epstol=1e-10;
disp(' ')
disp('itr     x1          x2          fval')
while 1
    totitern=totitern+1;
    for itr=1:num;  funval(itr)=PHI(var_p(itr,:)); itr=itr+1; end
    Fbar=max(funval);  fitness=Fbar-funval; Fsum=sum(fitness);
    opt_DV(totitern,:)=var_orig(find(funval==min(funval),1),:);  opt_funval(totitern)=f(opt_DV(totitern,:));
    fprintf('\n%d    %f     %f    %f\n', totitern,opt_DV(totitern,1),opt_DV(totitern,2),opt_funval(totitern));
    po=sort(fitness,'descend') ;  %arrange the chromosomes according to the fitness values
    
    X(:,:,totitern)=var_p;
    %%  printing values and plotting the updated figures
    subplot(1,2,1);  scatter(opt_DV(:,1),opt_DV(:,2),'filled','MarkerEdgeColor','b','MarkerFaceColor','r','LineWidth',1.5);
        strtxt1=['Optimum Values untill k= ' num2str(totitern) ' iteration' ]; diprnt1=text(-2,4.1,strtxt1);  set(diprnt1,'Fontsize',11,'Fontweight','bold','Fontangle','italic')
  title('All optimum values of Rastrigin function untill kth iteration using GA',...
    'FontName','Times','Fontsize',14,'FontWeight','bold');
xlabel('x1','FontName','Times','Fontsize',14,'FontWeight','bold');
ylabel('x2','FontName','Times','Fontsize',14,'FontWeight','bold');
    axis([-5 5 -5 5]); grid on
        ax=gca;
set(ax,'FontName','Times','Fontsize',12,'FontWeight','bold','Fontangle','italic');
set(gcf,'color',[1 1 1])
box on
     subplot(1,2,2);        
     scatter(X(:,1,totitern)-5,X(:,2,totitern)-5,'filled','MarkerEdgeColor','b','MarkerFaceColor','r','LineWidth',1.5);
       
       title('Optimization of Rastrigin function during kth iteration using GA',...
    'FontName','Times','Fontsize',14,'FontWeight','bold');
xlabel('x1','FontName','Times','Fontsize',14,'FontWeight','bold');
ylabel('x2','FontName','Times','Fontsize',14,'FontWeight','bold');
     
     strtxt1=['Optimum Values for']; diprnt1=text(2.2,3.7,strtxt1);  set(diprnt1,'Fontsize',12,'Fontweight','bold','Fontangle','italic')
     strtxt2=['iteration, k= ' num2str(totitern)]; diprnt2=text(2.2,3.7-.4,strtxt2);  set(diprnt2,'Fontsize',11,'Fontweight','bold','Fontangle','italic')
    strtxt3=['x1= ' num2str(opt_DV(totitern,1))]; diprnt3=text(2.2,2.8,strtxt3);  set(diprnt3,'Fontsize',11,'Fontweight','bold','Fontangle','italic')
    strtxt4=['x2= ' num2str(opt_DV(totitern,2))]; diprnt4=text(2.2,2.5,strtxt4);  set(diprnt4,'Fontsize',11,'Fontweight','bold','Fontangle','italic')
    strtxt5=['fmin= ' num2str(opt_funval(totitern))]; diprnt5=text(2.2,2.2,strtxt5);  set(diprnt5,'Fontsize',11,'Fontweight','bold','Fontangle','italic')
    timdnw=toc; strtxt6=['Ellapsed time= ' num2str(timdnw) 'secs']; diprnt6=text(2.2,1.9,strtxt6);  set(diprnt6,'Fontsize',11,'Fontweight','bold','Fontangle','italic')
    
    axis([-5 5 -5 5]); grid on
    
    ax=gca;
set(ax,'FontName','Times','Fontsize',12,'FontWeight','bold','Fontangle','italic');
set(gcf,'color',[1 1 1])
box on
    
    
    
    
    %% CONVERGENCE CRITERIA  using the maximum number of iterations, function convergence
    %convergence for maximum number of iterations
    if totitern>=Nmaxiter; disp('Convergence due to Maximum number of iterations!'); break; end
    
    %% convergence with absolute change in objective function
    if totitern>=2
        
        deltafun=opt_funval(totitern)-opt_funval(totitern-1);
        if abs(deltafun)>epstol; N1=0;
        else
            N1=N1+1; if N1>=repeat_convrg; disp('convergence due to absolute change in objective function!'); break; end
        end
        
        %% convergence with relative change in objective function
        if (abs(deltafun)/max(abs(opt_funval(totitern)),10^-5))>epstol; N2=0;
        else
            N2=N2+1; if N2>=repeat_convrg; disp('convergence due to relative change in objective function!'); break; end
        end
        
        %% convergence with change in design variables
        deltax=opt_DV(totitern,:)-opt_DV(totitern-1,:);
        if norm(deltax)>epstol; N3=0;
        else
            N3=N3+1; if N3>=repeat_convrg; disp('convergence due to change in design variables!'); break; end
        end
        
    end
    
    
    %% Evaluate fitness
    
    for kk=1:num
        loc(kk)=find(fitness==po(kk),1); chrom(kk,:)=var_p(loc(kk),:); kk=kk+1;
    end
    
   
    %% creating a rouletter wheel based on rank
    R0=0;
    for jjkk=1:num
        rr(jjkk)=2*(num-jjkk+1)/((num+1)*num)  ;  %here num=ns   number of samples
        if jjkk==1; R(jjkk)=R0+rr(jjkk); else R(jjkk)=R(jjkk-1)+rr(jjkk); end
        jjkk=jjkk+1;
        
    end
    %% SELECTION
    %first generate random number to satisfy the probability if Ps
    count=0;
    while(1)
        count=count+1;
        for jjmm=1:2
            rng('shuffle'); selec_rand=rand(1);
            for it=1:(num-1)
                if it==1; if R0<selec_rand<R(it+1); df=chrom(it+1,:);  break; end
                else if R(it)<selec_rand<R(it+1); df=chrom(it+1,:);  break; end; end
                it=it+1;
            end
            parent(jjmm,:)=df; chrom_num(jjmm)=it;
            jjmm=jjmm+1;
            
        end
        
        %% Converting the design variables to binary with the precision 0.001 to select one pair of parents each time
        chromosome1=strcat(num2str(dec2bin(parent(1,1)/precision,bits_num)),num2str(dec2bin(parent(1,2)/precision,bits_num)));
        chromosome2=strcat(num2str(dec2bin(parent(2,1)/precision,bits_num)),num2str(dec2bin(parent(2,2)/precision,bits_num)));
        
        %% CROSSOVER   here first use probaility for corssover
        rng('shuffle'); cross_rand=rand(1);
        if cross_rand<.98
            rng('shuffle'); L=2*bits_num; cutoff=randi(L-1);
            child1=strcat(num2str(chromosome1(1:cutoff)), num2str(chromosome2(cutoff+1:L)));
            child2=strcat(num2str(chromosome2(1:cutoff)), num2str(chromosome1(cutoff+1:L)));
        end
        %% choosing the best child from the crossover of a pair
        cc1=[bin2dec(child1(1:bits_num)), bin2dec(child1(bits_num+1:L))]*precision;
        cc2=[bin2dec(child2(1:bits_num)), bin2dec(child2(bits_num+1:L))]*precision;
        if f(cc1-5)<f(cc2-5); children(count,:)=child1; else children(count,:)=child2; end
        
        if count>=num;  break; end
    end
    
    %% MUTATION    here mutation is performed with very small probability about 0.01
    N=L; M=num; mutrate=.1; mutat=0; nmuts=round(mutrate*N*(M)); original=children;
    for ic=1:nmuts
        rng('shuffle'); ix=ceil(M*rand(1)); rng('shuffle'); iy=ceil(N*rand(1));
        val_01=num2str(children(ix,iy)); be4=children(ix,:);
        children(ix,iy)=num2str(1-str2num(children(ix,iy)));
        aftr=children(ix,:); mutat=mutat+1; ic=ic+1;
        
    end
    
    %% CONVERT BACK TO DECIMAL NUMBER THEN ORIGINAL VARIABLES X
    for itr=1:num
        xxpp(:,itr)=[bin2dec(children(itr,1:bits_num)), bin2dec(children(itr,bits_num+1:L))]*precision;
        itr=itr+1;
    end
    xp(:,:,totitern)=xxpp'; x=(xp-5);
    check1=xp(:,:,totitern);
    
    %% USING ELLTISM FOR COMBINATION OF PARENTS AND CHILDREN AND SELECT THE BEST
    if totitern==1; combine=var_p;
    elseif totitern>=2; combine=xp(:,:,totitern-1);
    end
    
    size(combine);
    for dupitr=1:length(combine); ellite_fun(dupitr)=f(combine(dupitr,:)-5); dupitr=dupitr+1; end
    ellite=sort(ellite_fun);
    
    for  dupitr22=1:5;
        xp((num-7)+dupitr22,:,totitern)=combine(find(ellite_fun==ellite(dupitr22),1),:);
        check=combine(find(ellite_fun==ellite(dupitr22),1),:);
        f(xp((num-7)+dupitr22,:,totitern)-5);
        dupitr22=dupitr22+1;
    end
    var_p= xp(:,:,totitern) ;%var_p=xp;  %10 is the maximum of DV
    var_orig=var_p-5;
    timdnw=toc; strtxt6=['Ellapsed time= ' num2str(timdnw) 'secs']; diprnt6=text(3.2,1.9,strtxt6);  set(diprnt6,'Fontsize',11,'Fontweight','bold','Fontangle','italic')
    
        waitbar((totitern)*globitr/(Nmaxiter*3)); perc=((totitern)*globitr/(Nmaxiter*3))*100;
    waitbar(perc/100,hh,sprintf('%d %% along for total iteration of %d...',perc,Nmaxiter*3));
    
    
end
    subplot(1,2,1);  scatter(opt_DV(totitern-1,1),opt_DV(totitern-1,2),'filled','MarkerEdgeColor','b','MarkerFaceColor','r','LineWidth',1.5);
diprnt1=text(opt_DV(totitern-1,1)+.1,opt_DV(totitern-1,2)+.1,'fmin');  set(diprnt1,'Fontsize',11,'Fontweight','bold','Fontangle','italic')


lth=length(opt_DV(:,1));
globopt(globitr,:)=opt_DV(totitern,:);
globmin(globitr)=opt_funval(totitern);
[globopt(globitr,1)   globopt(globitr,2)   globmin(globitr)]
if globitr>=3
    
break
end
globitr=globitr+1;
end


close(hh)

disp('x1   x2     fmin')
      globopt
            globmin'

    optimum_design=globopt(find(globmin==min(globmin),1),:)
optimum_value=min(globmin)
     scatter(optimum_design(1),optimum_design(2),'filled','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.5);
diprnt5=text(optimum_design(1)+.1,optimum_design(1)+.1,'fmin');  set(diprnt5,'Fontsize',11,'Fontweight','bold','Fontangle','italic')



toc

