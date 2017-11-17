# Numerical Methods HW3 Problem 2
# Methods
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


n<-c(100,110,120,130,140)
dx_vec<-as.numeric(length(n))
e11_vec<-as.numeric(length(n))
e12_vec<-as.numeric(length(n))
e13_vec<-as.numeric(length(n))

dx2_vec<-as.numeric(length(n))
e21_vec<-as.numeric(length(n))
e22_vec<-as.numeric(length(n))
e23_vec<-as.numeric(length(n))


for (j in 1:length(n))
{
  # dt = dx
  n1=n[j]
  T=0.5
  X=pi
  dt=T/n1
  dx=dt
  dx_vec[j]=dx
  lam<-dt/(dx^2)
  n2=floor(X/dx)
  sol=matrix(0, nrow=n1+1, ncol=n2+1)
  # populate first row
  for (i in 2:n2)
  {
    sol[1, i]=sin(i*dx)
  }
  dim_A<-n2-2+1
  
  sol11<-explicit(sol, dim_A, lam)
  sol12<-implicit(sol, dim_A, lam)
  sol13<-cn(sol, dim_A, lam)
  
  real_uT1<-numeric(dim_A)
  for (i in 1:dim_A)
  {
    real_uT1[i]<-sin(i*dx)*exp(-0.5)
  }
  
  # Compute errors
  e11_vec[j]=sqrt(sum(abs(sol11[(n1+1),2:(dim_A+1)]-real_uT1)^2*dx))
  e12_vec[j]=sqrt(sum(abs(sol12[(n1+1),2:(dim_A+1)]-real_uT1)^2*dx))
  e13_vec[j]=sqrt(sum(abs(sol13[(n1+1),2:(dim_A+1)]-real_uT1)^2*dx))
  
  # dt = 0.4*dx^2
  dx2<-sqrt(dt/0.4)
  dx2_vec[j]=dx2
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
  
  real_uT2<-numeric(dim_B)
  for (i in 1:dim_B)
  {
    real_uT2[i]<-sin(i*dx2)*exp(-0.5)
  }
  
  # Compute errors
  e21_vec[j]=sqrt(sum(abs(sol21[(n1+1),2:(dim_B+1)]-real_uT2)^2*dx2))
  e22_vec[j]=sqrt(sum(abs(sol22[(n1+1),2:(dim_B+1)]-real_uT2)^2*dx2))
  e23_vec[j]=sqrt(sum(abs(sol23[(n1+1),2:(dim_B+1)]-real_uT2)^2*dx2))
}

plot(dx_vec, e11_vec)
plot(dx_vec, e12_vec, type='b')
plot(dx_vec, e13_vec, type='b')
plot(dx2_vec, e21_vec, type='b')
plot(dx2_vec, e22_vec, type='b')
plot(dx2_vec, e23_vec, type='b')

# Use Monte Carlo method to solve the solution

# dt = dx
n1=100
T=0.5
X=pi
dt=T/n1
dx=dt
dx_vec[j]=dx
lam<-dt/(dx^2)
n2=floor(X/dx)
init_vec=numeric(n2-2+1)
# populate first row
for (i in 1:length(init_vec))
{
  init_vec[i]<-sin(i*dx)
}
# Monte Carlo Method
final_vec=numeric(n2-2+1)
for (i in 1:(n2-2+1))
{
  final_vec[i]<-mean(sin(init_vec[i]+rnorm(10000, 0, sqrt(2*T))))
}

# dt = 0.4*dx^2
dx2<-sqrt(dt/0.4)
dx2_vec[j]=dx2
lam2<-0.4
n2=floor(X/dx2)
init_vec1=numeric(n2-2+1)
# populate first row
# populate first row
for (i in 1:length(init_vec1))
{
  init_vec1[i]<-sin(i*dx2)
}
# Monte Carlo Method
final_vec1=numeric(n2-2+1)
for (i in 1:(n2-2+1))
{
  final_vec1[i]<-mean(sin(init_vec1[i]+rnorm(10000, 0, sqrt(2*T))))
}