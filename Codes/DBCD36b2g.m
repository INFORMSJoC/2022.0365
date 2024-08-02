function [ytrain1,b,ytrain]=DBCD36b2g(aa,aax,bb,x,y,lanf,c,qq)
n1=length(x(1,1,1,:,1));
n2=length(x(1,1,:,1,1));
n3=length(x(1,:,1,1,1));
n4=length(x(:,1,1,1,1));
c=c(1:n4,:,:);
x=x(1:n4,:,:,1,:);
tic
yy2=aa-aax+bb-lanf*(permute(repmat(c,1,1,1,n2),[1 2 4 3]));
for k=1:n2
    K=sum(exp(-(repmat(x(:,:,k,:,:),1,1,n2,n1)-repmat(x,1,1,1,n1)).^2./qq),5);
    yy3(:,:,:,:,k)=yy2.*K;
end 
ytrain=sum(yy3,5);
for i=1:n4
    for j=1:n3
        for t=1:n1
            for k=1:n2
                if(t==1)
                    y1=0;
                else
                    y1=y(i,j,k,1:t-1)-ytrain1(i,j,k,1:t-1);
                end
                b0(k)=y(i,j,k,t)+sum(y1)-ytrain(i,j,k,t);
            end
            b(i,j,t)=sum(b0)/n2;
            ytrain1(i,j,:,t)=ytrain(i,j,:,t)+permute(repmat(b(i,j,t),1,1,1,n2),[1 2 4 3]);
        end
    end
end
toc