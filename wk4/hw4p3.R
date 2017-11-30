#Euler scheme
euler<-function(sol, dim, dt, dx)
{
  #Construct the discrete interest rate vector
  R=c()
  for (i in 1:dim)
  {
    R=c(R, 5+i*dx)
  }
  #Construct transition matrix
  A=matrix(0, nrow=dim, ncol=dim)
  #Fill diagonal 
  for (i in 1:dim)
  {
    A[i, i]=1-((2.1^2)*R[i]*dt)/(dx^2)
  }
  #Fill upper diagonal
  for (i in 1:(dim-1))
  {
    A[i, i+1]=(sqrt(R[i])*dt)/(2*dx)+((2.1^2/2)*R[i]*dt)/(dx^2)
  }
  #Fill lower diagonal
  for (i in 2:dim)
  {
    A[i, i-1]=-(sqrt(R[i])*dt)/(2*dx)+((2.1^2/2)*R[i]*dt)/(dx^2)
  }
  for (i in 1:n1)
  {
    un_mat=matrix(sol[i,2:n2], nrow=dim, ncol=1)
    un1_mat=A%*%un_mat
    sol[i+1, 2:n2]=as.numeric(un1_mat)
  }
  return(sol)
}
#Parameters
n1=100000
T=1
dt=T/n1
X=15-5
dx=sqrt(dt/0.01)
n2=floor(X/dx)
sol=matrix(0, nrow=n1+1, ncol=n2+1)
# Populate the first row based on the terminal condition of the option
for (i in 2:n2)
{
  sol[1, i]=((5+(i-1)*dx)>=8&(5+(i-1)*dx)<=12)
}
dim=n2-2+1
#Use finite difference to solve the PDE
sol<-euler(sol, dim, dt, dx)
pv=sol[dim(sol)[1],ceiling(dim(sol)[2]/2)]
pv
#part ii part a
sens<-(0.08516807-0.1209332)/0.2
sens

