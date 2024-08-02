function [aa,aax,bb,o,ytrain1,ytest,ratio,ratiov,ratiov1,I]=DBCD60block2fp27gggLS(aa,aax,bb,x,y,C,B,p,lanf,c,qq,h,d)
n1=length(x(1,1,1,:,1));
n2=length(x(1,1,:,1,1));
n3=length(x(1,:,1,1,1));
n4=length(x(:,1,1,1,1));
c=c(1:n4,:,1:n1);
p=p(1:n4,1:n1);
h=h(1:n4,1:n1);
C=C(1:n4,:);
d=d(1:n1);
ytrain1=zeros(n4,n3,n2,n1);
aa1=zeros(n2,n4,n3,n1);aax1=zeros(n2,n4,n3,n1);bb1=zeros(n2,n4,n3,n1);ytest=ones(n4,n3,1,n1);o=zeros(n4,n3,n1);
ratio=permute(repmat(repmat(d,n4,1).*p,1,1,n3),[1 3 2])-c;
t=1;
B=B';
tic
for i=1:n4
    for j=1:n3
        for k=1:n1
            ratiov(1,t)=ratio(i,j,k);
            ratiov1(t,:)=[i j k];
            t=t+1;
        end
    end
end
toc
[~,I]=sort(ratiov,2,'descend');
for ii=1:(t-1)
    i=ratiov1(I(ii),1);
    j=ratiov1(I(ii),2);
    tt=ratiov1(I(ii),3);
    
            tic
            x0=x(i,j,:,1,:);
            x0=permute(x0,[3 5 1 2 4]);
            y0=y(i,j,:,tt);
            y0=permute(y0,[3 1 2 4]);
            n=length(y0);
            [H,k]=kernelg(x0,y0,qq);
            y1=y(i,j,:,1:tt-1)-ytrain1(i,j,:,1:tt-1);
            s0=(sign(-y1)+1)./2;
            y1=y1.*s0;
            y1=permute(y1,[1 3 4 2]);
            s=sign(tt-1);
            y2=y0'+s*sum(y1,3);
            b1=C(i,1)-sum(sum(sum(ytrain1(i,:,:,:),2),3),4);
            b2=B(j,1)-c(i,j,tt)*sum(sum(sum(ytrain1(1:n4,j,:,:),1),3),4);
            if(b1>sum(y2)&&b2>c(i,j,tt)*sum(y2))
                s1=(sign(y2)+1)/2;
                s2=sign(s1);
                f=[-y2.*s2-lanf*c(i,j,tt)*k,y2.*s2+lanf*c(i,j,tt)*k,-lanf*c(i,j,tt)*k];
                ib=zeros(3*n2,1);
                ub3=inf(n2,1);
                ub1=lanf*d(tt)*p(i,tt)*ones(n,1);
                ub2=lanf*h(i,tt).*ones(n,1);
                ub=[ub1
                    ub2
                    ub3];                  

                aeq=[ones(1,n2) -ones(1,n2) ones(1,n2)];
                beq=n2*lanf*c(i,j,tt);
                [a,fval]=quadprog(H,f,[],[],aeq,beq,ib,ub);   
                if(numel(a)==0)
                    aa1(:,i,j,tt)=zeros(n2,1);
                    aax1(:,i,j,tt)=zeros(n2,1);
                    bb1(:,i,j,tt)=zeros(n2,1);
                    o(i,j,tt)=nan;
                    aa(i,j,:,tt)=permute(aa1(:,i,j,tt),[2 3 1 4]);
                    aax(i,j,:,tt)=permute(aax1(:,i,j,tt),[2 3 1 4]);
                    bb(i,j,:,tt)=permute(bb1(:,i,j,tt),[2 3 1 4]);
                    ytrain1(i,j,:,tt)=zeros(1,1,n2,1);
                    ytest(i,j,1,tt)=0;
                else
                    aa1(:,i,j,tt)=a(1:n2,:);
                    aax1(:,i,j,tt)=a(n2+1:2*n2,:);
                    bb1(:,i,j,tt)=a(2*n2+1:3*n2,:);
                    o(i,j,tt)=fval;               
                    aa(i,j,:,tt)=permute(aa1(:,i,j,tt),[2 3 1 4]);
                    aax(i,j,:,tt)=permute(aax1(:,i,j,tt),[2 3 1 4]);
                    bb(i,j,:,tt)=permute(bb1(:,i,j,tt),[2 3 1 4]);
                    ytrain1(i,j,:,tt)=DBCD32bfpggLS(aa,aax,bb,x,y,lanf,c,qq,tt,j,i,ytrain1);
                    ytrain1(i,j,:,tt)=((sign(ytrain1(i,j,:,tt))+1)/2).*ytrain1(i,j,:,tt);%取正数
                end
            elseif(b1>0&&b2>0)
                y21=b1.*ones(1,n2)/n2;
                y22=(b2/c(i,j,tt)).*ones(1,n2)/n2;
                y2=min(y21,y22);
                f=[-y2-lanf*c(i,j,tt)*k,y2+lanf*c(i,j,tt)*k,-lanf*c(i,j,tt)*k];
                ib=zeros(3*n2,1);
                ub3=inf(n2,1);
                ub1=lanf*d(tt)*p(i,tt)*ones(n,1);
                ub2=lanf*h(i,tt).*ones(n,1);
                ub=[ub1
                    ub2
                    ub3];                  

                aeq=[ones(1,n2) -ones(1,n2) ones(1,n2)];
                beq=n2*lanf*c(i,j,tt);
                [a,fval]=quadprog(H,f,[],[],aeq,beq,ib,ub);   
                if(numel(a)==0)
                    aa1(:,i,j,tt)=zeros(n2,1);
                    aax1(:,i,j,tt)=zeros(n2,1);
                    bb1(:,i,j,tt)=zeros(n2,1);
                    o(i,j,tt)=nan;
                    aa(i,j,:,tt)=permute(aa1(:,i,j,tt),[2 3 1 4]);
                    aax(i,j,:,tt)=permute(aax1(:,i,j,tt),[2 3 1 4]);
                    bb(i,j,:,tt)=permute(bb1(:,i,j,tt),[2 3 1 4]);
                    ytrain1(i,j,:,tt)=zeros(1,1,n2,1);
                    ytest(i,j,1,tt)=0;
                else
                    aa1(:,i,j,tt)=a(1:n2,:);
                    aax1(:,i,j,tt)=a(n2+1:2*n2,:);
                    bb1(:,i,j,tt)=a(2*n2+1:3*n2,:);
                    o(i,j,tt)=fval;               
                    aa(i,j,:,tt)=permute(aa1(:,i,j,tt),[2 3 1 4]);
                    aax(i,j,:,tt)=permute(aax1(:,i,j,tt),[2 3 1 4]);
                    bb(i,j,:,tt)=permute(bb1(:,i,j,tt),[2 3 1 4]);
                    ytrain1(i,j,:,tt)=DBCD32bfpggLS(aa,aax,bb,x,y,lanf,c,qq,tt,j,i,ytrain1);
                    ytrain1(i,j,:,tt)=((sign(ytrain1(i,j,:,tt))+1)/2).*ytrain1(i,j,:,tt);%取正数
                end
            else
                ytrain1(i,j,:,tt)=zeros(1,1,n2,1);
                ytest(i,j,1,tt)=0;
            end
            toc 
end