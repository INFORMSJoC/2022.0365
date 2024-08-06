function [ytrain1,b]=DBCD32bfpgg(aa,aax,bb,x,y,lanf,c,qq,tt,j,i,ytrain0)
n2=length(x(1,1,:,1,1));
yy2=aa(i,j,:,tt)-aax(i,j,:,tt)+bb(i,j,:,tt)-repmat(lanf*c(i,j,tt),1,1,n2);
for k=1:n2
    K=sum(exp(-(repmat(x(i,j,k,1,:),1,1,n2)-x(i,j,:,1,:)).^2./qq),5);
    yy3(:,:,:,:,k)=yy2.*K;
end
ytrain=sum(yy3,5);
if(tt==1)
    y1=0;            
else
    y1=y(i,j,:,1:tt-1)-ytrain0(i,j,:,1:tt-1);
end
b0=permute(y(i,j,:,tt)+sum(y1,4),[3 1 2 4])-ytrain;
b=sum(b0)/n2;
ytrain1=ytrain+b;