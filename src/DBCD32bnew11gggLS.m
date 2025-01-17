function ytrain1=DBCD32bnew11gggLS(aa,aax,bb,r,th,x,y,lanf,c,qq,p,tt,j,i,ytrain0,ytrain2,ytrain3)
n1=length(aa(1,1,1,:));
n2=length(x(1,1,:,1,1));
yy2=aa(i,j,:,tt)-aax(i,j,:,tt)+bb(i,j,:,tt)-repmat(r(i,1),1,1,n2)-repmat(p(i,tt),1,1,n2).*repmat(th(j,1),1,1,n2)-repmat(lanf*c(i,j,tt),1,1,n2);
for k=1:n2
    K=sum(exp(-(repmat(x(i,j,k,tt,:),1,1,n2)-x(i,j,:,tt,:)).^2./qq),5);
    yy3=yy2.*K;
end
ytrain=sum(yy3,5);
for k=1:n2
    if(tt==n1)
        y1=0;            
    else
        y1=y(i,j,k,tt+1:n1)-ytrain0(i,j,k,tt+1:n1);
        s0=(sign(y1)+1)./2;
        y1=y1.*s0;
    end
    o1=max(ytrain2(1,k),0);
    o2=max(ytrain3(1,k),0);
    o3=min(o1,o2);
    o4=min(o3,y(i,j,k,tt)+sum(y1));
    b0(k)=o4-ytrain(1,1,k,1);
end
b=sum(b0)/n2;
ytrain1=ytrain+repmat(b,1,1,n2,1);

