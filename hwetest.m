function pvalue=hwetest(x,varargin)
%HWETEST tests if a population is in the Hardy Weinberg Proportion (HWP) for
%the observed locus.
%The conditional probability, under the Hardy-Weinberg
%equilibrium, to obtain the sample X is compute as described by Howard
%Levene - "On a matching problem arising in genetics".
%Annals of Mathematical Statistics. 1949; 20:91-94.
%             ____ m             __
%         N!   || i=1 Ai!        \
% Pr(f)=---------------------- 2^/_ X(i,j)
%             ____                 j>i
%       (2N)!  ||j>=i X(i,j)!
%
%assuming Ai = the count of i-th allele
%N=genotypes in the matrix.
%
%If the locus is biallelic, the function performs an Exact test, computing
%the p-value of all possible tables and the summing all p-value<=p(observed
%table) and plots a De Finetti's Diagram.
%If the locus is m-allelic (m>2), the function uses a Monte Carlo conventional
%method to evaluate the p-value.
%
%Syntax: pvalue=hwetest(x,verbose,delta,alpha)
%
%Input: X - Genotype matrix. If the locus is biallelic, X is a vector
%           x=[AA AB BB]; else if the locus is m-allelic X is a lower
%           triangular matrix of size=[m m]. If X is not a lower
%           triangular matrix it will be triangularized.
%       VERBOSE (optional)- a logical variable to display more results and comments:
%              0=does not display (default)
%              1=display
%       DELTA and ALPHA (optional)- If Monte Carlo method is used (if locus is more
%           than bi-allelic), it is necessary to evaluate how many times
%           the process must be reiterated to ensure that p-value is
%           within DELTA units of the true one with (1-ALPHA)*100% confidence.
%           (Default DELTA=ALPHA=0.01).
%Output: the probability that the population is in HWP
%        the De Finetti's Diagram if the locus is biallelic
%        if VERBOSE:
%           Polymorphism Information Content (PIC)
%           Matching probability
%           Power of discrimination
%           Power of exclusion
%           Typical Paternity Index
%
%Example:
%          Run hwedemo
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2007) HWtest: a routine to test if a locus is in Hardy
% Weinberg equilibrium (exact test).
% http://www.mathworks.com/matlabcentral/fileexchange/14425

%Input error Handling
p = inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'2d','real','finite','nonnan','nonempty','integer'}));
addOptional(p,'verbose',0, @(x) isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==0 || x==1));
addOptional(p,'delta',0.01, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p,'alpha',0.01, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
parse(p,x,varargin{:});
verbose=p.Results.verbose; alpha=p.Results.alpha; delta=p.Results.delta;
clear p

%check the matrix
[rs,cs]=size(x);
if rs==1 && cs==3
    x=[x(1) 0; x(2) x(3)];
    cs=2;
else
    assert(isequal(rs,cs),'Warning: Matrix must be a square matrix')
end

tr=repmat('-',1,60); %spacer

%If x is not a lower triangular matrix, then make it triangular
if ~isempty(find(triu(x,1),1))
    x=tril(x)+triu(x,1)';
    if verbose
        disp('The genotype X(I,J) is the same that X(J,I);')
        disp('so the matrix must be triangularized.')
        disp('Matrix Triangularization')
        disp(x)
        disp('Press a key to continue'); pause; clc; home
    end
end

N=sum(x(:)); %Total Genotypes
al=sum(x)+sum(x,2)'; %Alleles array

%If there are alleles with null frequency, then delete them from the matrix
z=find(al==0);
if ~isempty(z)
    x(z,:)=[]; %delete the row(s)
    x(:,z)=[]; %delete the column(s)
    al(z)=[]; %delete the allele(s)
    cs=cs-length(z); %decrease the number of columns...
    if verbose
        disp('There is one or more allele with a null frequency, so it will be deleted.')
        disp(['Alleles with null frequency: ' num2str(z)])
        disp('Deleting...')
        disp(x)
        disp(' ')
        disp('Press a key to continue'); pause; clc; home
    end
    clear z
end

fa=al./(2*N); %Allelic frequencies
e=sqrt(prod(fa)/N); %error of frequencies determination
mfa=(fa'*fa); %matrix of binomial expansion of allelic frequencies
xe=round(mfa.*N); %Expected matrix under HWP
xe=tril(xe)+triu(xe,1)'; %triangularize the matrix

if cs==2 %if this is a biallelic locus...
    if verbose
        disp(array2table(x([1,2,4]),'RowNames',{'Observed_Genotypes'},'VariableNames',{'AA','AB','BB'}))
        disp(array2table([fa e],'RowNames',{'Observed_Frequencies'},'VariableNames',{'A_Allele','B_Allele','Error'}))
        disp(array2table(xe([1,2,4]),'RowNames',{'Expected_Genotypes'},'VariableNames',{'AA','AB','BB'}))
    end
    R=min(al); %Rare Allele
    H=mod(R,2):2:R; %all possible heterozygotes
    RR=(R-H)/2; %all possible Rare allele homozygotes
    CC=(2*N-R-H)/2; %all possible Common allele homozygotes
    %Considering that X!=gamma(X+1) and log(X!)=gammaln(X+1) it
    %is possible vectorizing all the test using the gammaln Matlab function.
    %This is about 400x faster than for...end loop
    %Moreover the properties of logarithms must be taken in account:
    %Log(a*b)=log(a)+log(b) and Log(a/b)=log(a)-log(b)
    %Constant factor: (N!C!R!)/[(2N)!]
    KF=(sum(gammaln([N al]+ones(1,3)))-gammaln(2*N+1));
    %Costant array: 2.^H. Considering the properties of logarithms
    %log(2.^H)=H.*log(2)
    KA=(H.*log(2));
    %Construct the matrix of all possible combination of genotypes
    tbs=[CC; H; RR;];
    %Compute the log(x!) by the gammaln function and sum them
    tbsfac=sum(gammaln(tbs+ones(size(tbs))));
    %Finally compute the probability of each combinations
    pg=exp(KF+KA-tbsfac);
    %Sum all pg less or equal than p of the observed table
    p=sum(pg(pg<=pg(H==x(2))));
    %display results
    disp(array2table([length(H),p],'RowNames',{'Hardy_Weinberg_Proportion'},'VariableNames',{'Tables','p_value'}));
    if verbose
        hFig=figure;
        units=get(hFig,'units');
        set(hFig,'units','normalized','outerposition',[0 0 1 1], 'Color', 'white');
        set(hFig,'units',units); clear units
        axis square
        hold on
        set(gca,'Xlim',[-0.05 1.05],'Ylim',[0 1.1])
        plot([-0.05 0],[0 0],'Color','w')
        plot([1 1.05],[0 0],'Color','w')
        plot([0 0.5],[0 1],'Color','k','Linewidth',1)
        plot([0.5 1],[1 0],'Color','k','Linewidth',1)
        passo=0.1:0.1:0.9;
        xa=0.5.*passo; ya=passo;
        xb=xa+0.5; yb=1-passo;
        for I=1:9
            plot([passo(I) xa(I)],[0 ya(I)],'Color',[169 169 169]./255,'Linewidth',1)
            plot([passo(I) xb(I)],[0 yb(I)],'Color',[169 169 169]./255,'Linewidth',1)
            plot([xa(I) 1-xa(I)],[ya(I) ya(I)],'Color',[169 169 169]./255,'Linewidth',1)
        end
        clear passo xa xb ya yb
        if fa(1)>0.5
            yH=-2*fa(1)+2;
        elseif fa(1)<0.5
            yH=2*fa(1);
        else
            yH=1;
        end
        %           Now plot in red the tables outside the Hardy-Weinberg equilibrium
        %           and in green the tables inside the Hardy-Weinberg equilibrium.
        
        y=cumsum(pg); %the sum of all probability is 1...
        %Find the interval containing the 95% (or, if you want, find the points where the 2.5% tails starts)
        xg=linspace(0,1,500);
        h1=plot(xg,2.*(xg.*(1-xg)),'Color','b','Linewidth',2); %plot the Hardy-Weinberg Parabola (in blue) in a ternary plot
        clear xg
        norma=linspace(0,yH,length(H));
        obs=norma(H==x(2));
        h2=plot([fa(1) fa(1)],[0 obs],'Color','r','Linewidth',2);
        xp=-2/5*(-2-0.5*fa(1)+obs); yp=-2*xp+2;
        h3=plot([fa(1) xp],[obs yp],'Color','c','Linewidth',2);
        xp=2/5*(0.5*fa(1)+obs); yp=2*xp;
        h4=plot([fa(1) xp],[obs yp],'Color','m','Linewidth',2);
        h5=plot([fa(1) fa(1)],[norma(find(y<=0.025,1,'last')) norma(find(y>=0.975,1,'first'))],'Color','g','Linewidth',2);
                %Then put a spot in corrispondence of the observed table.
        h6=plot(fa(1),obs,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',6);
        hold off
        set(findobj(hFig, 'type','axes'),'Ycolor','w') %turn-off y axis
        title('De Finetti''s Diagram')
        legend([h1,h2,h3,h4,h5,h6],'HWP parabola','Lenght=AB frequencies','Lenght=BB frequencies','Lenght=AA frequencies','HW equilibrium (alpha=0.05)','Observed genotypes')
        text(-0.05,0,'AA'); text(1.02,0,'BB'); text(0.485,1.02,'AB')
    end
else %if there are more than 2 allelels...
    %Perform Monte Carlo method
    if verbose
        A=cell(1,cs);
        for I=1:cs
            A{I}=strcat('A',num2str(I));
        end
        disp([blanks(21) 'Observed Genotypes']);
        disp(tr)
        disp(array2table(x,'RowNames',A,'VariableNames',A))
        disp([blanks(15) 'Observed Frequencies and Error']);
        disp(tr)
        B=A; B{cs+1}='Error';
        disp(array2table([fa e],'VariableNames',B))
        disp([blanks(16) 'Expected Genotypes under HWP']);
        disp(tr)
        disp(array2table(xe,'RowNames',A,'VariableNames',A))
        clear A B
    end
    %If observed and expected tables are equal don't perform the Monte
    %Carlo simulation.
    if isequal(x,xe)
        tbs=0;
        p=1;
    else
        %Monte Carlo Naive method
        %as described in:
        %Guo SW and Thompson EA. Performing the exact test of
        %Hardy-Weinberg proportion for multiple alleles.
        %Biometrics. 1992; 48: 361-372.
        s=zeros(1,2*N); %gametes array preallocation
        count=1;
        for C=1:cs
            s(count:count+al(C)-1)=C; %construct gametes array
            count=count+al(C);
        end
        %See above for explanation.
        %costant factor N!/(2N!)
        f1=sum([gammaln(N+1) -gammaln(2*N+1)]);
        %costant array prod(alleles!)
        f2=sum(gammaln(al+ones(size(al))));
        f=f1+f2;
        %costant array 2^sum(Heterozygotes)
        f3=sum(sum(tril(x,-1)))*log(2);
        %Compute the log(x!) by the gammaln function and sum them
        f4=sum(sum(gammaln(x+ones(size(x)))));
        pf=exp(f+f3-f4); %p-value of the observed matrix
        %tbs=simulation size to ensure that p-value is within delta units of the
        %true one with (1-alpha)*100% confidence.
        %Psycometrika 1979; Vol.44:75-83.
        tbs=round(((-realsqrt(2)*erfcinv(2-alpha))/(2*delta))^2);
        K=0; %Monte Carlo counter
        for C=1:tbs
            %shuffle the gametes array using the Fisher-Yates shuffle
            %Sattolo's version. This is faster than Matlab RANDPERM: to be
            %clearer: Fisher-Yates is O(n) while Randperm is O(nlog(n))
            for J=2*N:-1:2
                k=ceil((J-1).*rand);
                tmp=s(k);
                s(k)=s(J);
                s(J)=tmp;
            end
            g=zeros(cs,cs); %Construct a new table
            %This cycle is faster than Matlab ACCUMARRAY and allow the
            %costruction of a correct matrix. In fact it is possible that
            %from the permutation some row or column doesn't exist (i.e.
            %the last homozygote cell is empty and so the matrix isn't a
            %square matrix). To avoid this occurrence due to ACCUMARRAY
            %function, you should introduce an IF or SWITCH check, slowing
            %the execution.
            for J=1:N
                A=s(2*J-1); %take the row index
                B=s(2*J); %take the column index
                g(A,B)=g(A,B)+1; %add one to the cell
            end
            g=tril(g)+triu(g,1)'; %triangularize the matrix
            %costant array 2^sum(Heterozygotes)
            f3=sum(sum(tril(g,-1)))*log(2);
            %Compute the log(x!) by the gammaln function and sum them
            f4=sum(sum(gammaln(g+ones(size(g)))));
            pg=exp(f+f3-f4); %compute the p-value of the new matrix
            if pg<=pf
                K=K+1; %update counter
            end
        end
        p=K/tbs; %Monte Carlo p-value
    end
    %display results
    disp('Conventional Monte Carlo Method')
    disp(tr)
    disp(array2table([tbs,p],'RowNames',{'Hardy_Weinberg_Proportion'},'VariableNames',{'Tables','p_value'}));
    fprintf('p-value is within %0.4f units of the true one with %0.4f%% confidence\n',delta,(1-alpha)*100)
    disp(tr)
end
if verbose
    %Computes the observed heterozygosis (Ho), the expected
    %heterozygosis (He) and the inbreeding coefficient (F)
    HO=1-trace(x)/N; %Observed Heterozygosity
    HE=1-trace(xe)/N; %Expected Heterozygosity
    F=1-HO/HE; %Inbreeding coefficient
    %display results
    if F==0
        txt='Perfect fit';
    else
        switch sign(F)
            case 1
                txt='Fewer heterozygotes than expected (inbreeding)';
            case -1
                txt='More heterozygotes than expected (outbreeding)';
        end
    end
    disp(cell2table({HO, HE, F, txt},'VariableNames',{'Obs_Het','Exp_Het','F_coefficient','Comment'}));
    xpic=diag(mfa);
    mpic=xpic*xpic';
    mpic=tril(mpic)+triu(mpic,1)'; %triangularize the matrix
    PIC=1-sum(xpic)-sum(sum(tril(mpic,-1))); %Polymorphism Information Content
    MP=sum(sum((x./N).^2)); %matching probability
    PD=1-MP; %power of discrimination
    OH=1-HO; %observed homozigosity;
    PE=HO^2*(1-2*HO*OH^2); %power of exclusion
    TPI=1/(2*OH); %typical paternity index
    disp(array2table([PIC, MP, PD, PE, TPI ],'VariableNames',{'PIC','MP','PD','PE','TPI'}));
    disp('PIC=Polymorphism Information Content; MP=Matching probability')
    disp('PD=Power of discrimination; PE=Power of exclusion; TPI=Typical Paternity Index')
end

if nargout
    pvalue=p;
end