# Numerical Methods HW3 Problem 2
# dt = dx
n1=100
T=0.5
X=pi
dt=T/n
dx=dt
lam<-dt/(dx^2)
n2=floor(X/dx)
sol=matrix(0, nrow=n1+1, ncol=n2+1)
# populate first row
for (i in 2:n2)
{
  sol[1, i]=sin(i*dx)
}
dim_A<-n2-2+1

# Use explicit method to solve
explicit<-function(sol, dim_A, lam)
{
  # Construct A matrix
  A=matrix(0, nrow=dim_A, ncol=dim_A)
  for (i in 1:dim_A)
  {
    A[i, i]=1-2*lam
  }
  for (i in 1:(dim_A-1))
  {
    A[i, i+1]=lam
  }
  for (i in 2:dim_A)
  {
    A[i, i-1]=lam
  }
  for (i in 1:n1)
  {
    un_mat<-matrix(sol[i,2:n2], nrow=dim_A, ncol=1)
    un1_mat<-A%*%un_mat
    sol[i+1,2:n2]=as.numeric(un1_mat)
  } 
  return(sol)
}


# Use implicit method to solve
implicit<-function(sol, dim_B, lam)
{
  # Construct B matrix
  B=matrix(0, nrow=dim_B, ncol=dim_B)
  for (i in 1:dim_B)
  {
    B[i, i]=1+2*lam
  }
  for (i in 1:(dim_B-1))
  {
    B[i, i+1]=-lam
  }
  for (i in 2:dim_B)
  {
    B[i, i-1]=-lam
  }
  for (i in 1:n1)
  {
    un_mat<-matrix(sol[i,2:n2], nrow=dim_B, ncol=1)
    mm=solve(B)
    un1_mat<-mm%*%un_mat
    sol[i+1,2:n2]=as.numeric(un1_mat)
  }
  return(sol)
}


# Use Crank-Nicolson method to solve
cn<-function(sol, dim_C, lam)
{
  # Construct C and D matrices
  C=matrix(0, nrow=dim_C, ncol=dim_C)
  D=matrix(0, nrow=dim_C, ncol=dim_C)
  for (i in 1:dim_C)
  {
    C[i, i]=1+lam
    D[i, i]=1-lam
  }
  for (i in 1:(dim_C-1))
  {
    C[i, i+1]=-lam/2
    D[i, i+1]=lam/2
  }
  for (i in 2:dim_C)
  {
    C[i, i-1]=-lam/2
    D[i, i-1]=lam/2
  }
  for (i in 1:n1)
  {
    un_mat<-matrix(sol[i,2:n2], nrow=dim_C, ncol=1)
    mm<-solve(C)%*%D
    un1_mat<-mm%*%un_mat
    sol[i+1,2:n2]=as.numeric(un1_mat)
  }
  return(sol)
}

sol11<-explicit(sol, dim_A, lam)
sol12<-implicit(sol, dim_A, lam)
sol13<-cn(sol, dim_A, lam)

# dt = 0.4*dx^2
dx2<-sqrt(dt/0.4)
lam2<-0.4
n2=floor(X/dx2)
sol1=matrix(0, nrow=n1+1, ncol=n2+1)
# populate first row
for (i in 2:n2)
{
  sol1[1, i]=sin(i*dx2)
}
dim_B<-n2-2+1

sol21<-explicit(sol1, dim_B, lam2)
sol22<-implicit(sol1, dim_B, lam2)
sol23<-cn(sol1, dim_B, lam2)
