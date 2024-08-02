function [aa,aax,bb,o,x0,y0]=DBCD60block2bpnew15LS(x,y,r,th,p,lanf,c,qq,h,d)
%内点法，包含偏置项，分块，反向算法，算法修改为反向迭代库存更新，删除if语句，考虑订货量和需求不为负数，多步预测,修改b的计算
n1=length(x(1,1,1,:,1));
n2=length(x(1,1,:,1,1));
n3=length(x(1,:,1,1,1));
n4=length(x(:,1,1,1,1));
ytrain1=zeros(n4,n3,n2,n1);
for i=1:n4
    for j=1:n3
        for tt=1:n1
            t=n1-tt+1;
            tic
            x0=x(i,j,:,1,:);
            x0=permute(x0,[3 5 1 2 4]);
            y0=y(i,j,:,t);
            y0=permute(y0,[3 1 2 4]);
            n=length(y0);
            [H,k]=kernelg(x0,y0,qq);
            y1=y(i,j,:,t+1:n1)-ytrain1(i,j,:,t+1:n1);
            s0=(sign(y1)+1)./2;
            y1=y1.*s0;
            y1=permute(y1,[1 3 4 2]);
            s=sign(n1-t);
            y2=y0'+s*sum(y1,3);
            s1=(sign(y2)+1)/2;
            s2=sign(s1);
            f=[-y2.*s2-(r(i,1)+p(i,t)*th(j,1)+lanf*c(i,j,t))*k,y2.*s2+(r(i,1)+p(i,t)*th(j,1)+lanf*c(i,j,t))*k,-(r(i,1)+p(i,t)*th(j,1)+lanf*c(i,j,t))*k];
            ib=zeros(3*n2,1);
            ub3=inf(n2,1);
            ub1=lanf*d(t)*p(i,t)*ones(n,1);
            ub2=lanf*h(i,t).*ones(n,1);
            ub=[ub1
                ub2
                ub3]; 
            aeq=[ones(1,n2) -ones(1,n2) ones(1,n2)];
            beq=n2*(r(i,1)+p(i,t)*th(j,1)+lanf*c(i,j,t));
            [a,fval]=quadprog(H,f,[],[],aeq,beq,ib,ub);            
            aa1(:,i,j,t)=a(1:n2,:);
            aax1(:,i,j,t)=a(n2+1:2*n2,:);
            bb1(:,i,j,t)=a(2*n2+1:3*n2,:);
            o(i,j,t)=fval;           
            aa(i,j,:,t)=permute(aa1(:,i,j,t),[2 3 1 4]);
            aax(i,j,:,t)=permute(aax1(:,i,j,t),[2 3 1 4]);
            bb(i,j,:,t)=permute(bb1(:,i,j,t),[2 3 1 4]);
            ytrain1(i,j,:,t)=DBCD32bnew11ggLS(aa,aax,bb,r,th,x,y,lanf,c,qq,p,t,j,i,ytrain1);
            toc
        end
    end    
end