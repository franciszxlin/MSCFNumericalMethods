#homework 4 problem 2
#Crank-Nicolson scheme finite difference scheme
cn<-function(sol, dim, lam)
{
  # Construct C and D matrices
  C=matrix(0, nrow=dim, ncol=dim)
  D=matrix(0, nrow=dim, ncol=dim)
  for (i in 1:dim)
  {
    C[i, i]=1+lam
    D[i, i]=1-lam
  }
  for (i in 1:(dim-1))
  {
    C[i, i+1]=-lam/2
    D[i, i+1]=lam/2
  }
  for (i in 2:dim)
  {
    C[i, i-1]=-lam/2
    D[i, i-1]=lam/2
  }
  for (i in 1:n1)
  {
    un_mat<-matrix(sol[i,2:n2], nrow=dim, ncol=1)
    mm<-solve(C)%*%D
    un1_mat<-mm%*%un_mat
    sol[i+1,2:n2]=as.numeric(un1_mat)
  }
  return(sol)
}
#Parameters: we will discretize the domain into grid
n1=500
n2=700
T=1
X=1
dt=T/n1
dx=X/n2
lam=dt/dx^2
sol=matrix(0, nrow=n1+1, ncol=n2+1)
#Populate the first row
for (i in 2:n2)
{
  sol[1, i]=exp(i*dx)
}
#Solve the PDE
dim=n2-2+1
sol=cn(sol, dim, lam)
head(sol[ ,1:5])
tail(sol[ ,1:5])
#Part ii Calculate rowsums
rowsum<-rowSums(sol)*dx
plot(0:500, rowsum)

