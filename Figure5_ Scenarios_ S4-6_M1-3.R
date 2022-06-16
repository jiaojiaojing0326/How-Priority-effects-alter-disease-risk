# AUTHOR:   JING JIAO & MICHAEL H. CORTEZ
# DATE:     February 2022
# PURPOSE:  Simulate relationships between specialist and generalist infectious propagule
#           densities for pairwise combinations of Scenarios S4-S6 and M1-M3

# INPUTS:   Parameters for each host and pathogen

# METHOD:   Step 1: Define equations for ODE model
#           Step 2: Define global simulation parameters
#           Step 3: Compute equilibria for scenarios S1-S3
#           Step 4: Compute equilibria for scenarios M1-M3
#           Step 5: Compute equilibria for pairwise combinations of S and M scenarios
#           Step 6: Plot P1-P2 relationships

# OUTPUTS:  Equilibrium values at different delta2 values (generalist degradation rate)
#           

rm(list=ls()) 

library(deSolve)


#########################################################################
### STEP 1: DEFINE EQUATIONS FOR ODE MODEL

DiffEqns <- function(t, inits, parameters) {
  with(as.list(c(inits, parameters)), {                                           
    dI_F_1 <- beta_F_1 * (N_F - I_F_1 - I_F_2 - I_F_12 - I_F_21) * P1 - beta_F_12 * I_F_1 * P2 - m_F_1 * I_F_1  #Focal host singly-infected by specialist
    dI_F_2 <- beta_F_2 * (N_F - I_F_1 - I_F_2 - I_F_12 - I_F_21) * P2 - beta_F_21 * I_F_2 * P1 - m_F_2 * I_F_2  #Focal host singly-infected by generalist
    dI_A_2 <- beta_A_2 * (N_A - I_A_2) * P2 - m_A_2 * I_A_2 #Alternative host singly-infected by generalist
    dI_F_12 <- beta_F_12 * I_F_1 * P2 - m_F_12 * I_F_12 #Focal host co-infected by specialist first
    dI_F_21 <- beta_F_21 * I_F_2 * P1 - m_F_21 * I_F_21 #Focal host co-infected by specialist second
    dP1 <- chi_F_1 * I_F_1 + chi_F_12 * I_F_12 + chi_F_21 * I_F_21 - (u_F_1*N_F + u_A_1*N_A) * P1 - delta1 * P1  #Specialist infectious propagule density
    dP2 <- phi_A_2 * I_A_2 + phi_F_2 * I_F_2 + phi_F_12 * I_F_12 + phi_F_21 * I_F_21 - (u_F_2*N_F + u_A_2*N_A) * P2 - delta2 * P2  #Generalist infectious propagule density
    list(c(dI_F_1, dI_F_2, dI_A_2, dI_F_12, dI_F_21, dP1, dP2))
  })
}


#########################################################################
### STEP 2: DEFINE GLOBAL SIMULATION PARAMETERS

# Run each simulation for 5000 time units
t_final <- 5000  
times <- seq(0, t_final, by = 1)

# Initial conditions for simulations
initial_conditions <- c(I_F_1 = 0, I_F_2 = 0, I_A_2 = 0, I_F_12 = 0, I_F_21 = 0, P1 = 10, P2 = 10)

# Values of generalist degradation rate
delta2_vals <- c(800000, 600000, 550000, 500000, 450000, 400000, 300000, 200000, 100000, 50000, 25000, 12500, 6000, 3000, 1500, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100)

###############################################

# Parameters shared across simulations
parameters_shared<-c(
  N_F = 200,
  N_A = 200,
  
  beta_F_2 <- 0.01, 
  beta_F_1 <- 0.01, 
  beta_F_12 <- 0.01,
  beta_F_21 <- 0.01,
  beta_A_2 <- 0.01,
  
  m_F_2 <- 0.5,  
  m_A_2 <- 0.5,
  
  phi_F_2 <- 2000,
  phi_F_12 <- 2000,
  phi_F_21 <- 2000,
  phi_A_2 <- 2000,
  
  u_F_1 <- 0.1, 
  u_A_1 <- 0.1,
  u_F_2<- 0.1,
  u_A_2 <- 0.1,
  delta1 <- 100
)


#########################################################################
### Step 3: COMPUTE EQUILIBRIA FOR EACH SCENARIO


## SCENARIO S4 
# Parameter values

chi_F_1_vals = 1500
chi_F_12_vals = 2000
chi_F_21_vals = 3000 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(chi_F_12_vals))
{
  I1_1.eq <- c()
  I1_2.eq <- c()
  I2_2.eq <- c()
  I1_12.eq <- c()
  I1_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      m_F_1 <- 0.5, 
      m_F_21 <- 0.5, 
      m_F_12 <- 0.5, 
      chi_F_1 <- chi_F_1_vals[h], 
      chi_F_12 <- chi_F_12_vals[h],
      chi_F_21 <- chi_F_21_vals[h],
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I1_1.eq[i] <- out[t_final+1, 2]
    I1_2.eq[i] <- out[t_final+1, 3]
    I2_2.eq[i] <- out[t_final+1, 4]
    I1_12.eq[i] <- out[t_final+1, 5]
    I1_21.eq[i] <- out[t_final+1, 6]
    P1.eq[i] <- out[t_final+1, 7]
    P2.eq[i] <- out[t_final+1, 8]
  }
  
  ## save propagule densities
  dd1 <- P1.eq
  record1 <- rbind(record1,dd1)
  
  dd2 <- P2.eq
  record2 <- rbind(record2,dd2)
}


colnames(record1)<- NULL 
colnames(record2)<- NULL

A = which(is.na(record1), arr.ind = TRUE)

B = which(is.na(record2), arr.ind = TRUE)

index = union(A, B)

if(length(index) == 0)
{

ScenarioS4_P1 <- record1
ScenarioS4_P2 <- record2
}else
{
  ScenarioS4_P1 <- record1[-index]
  ScenarioS4_P2 <- record2[-index]
}

##########################
## SCENARIO S5 
# Parameter values

chi_F_1_vals = 2000
chi_F_12_vals = 1500
chi_F_21_vals = 3000 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(chi_F_12_vals))
{
  I1_1.eq <- c()
  I1_2.eq <- c()
  I2_2.eq <- c()
  I1_12.eq <- c()
  I1_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      m_F_1 <- 0.5, 
      m_F_2 <- 0.5,  
      m_F_21 <- 0.5, 
      m_F_12 <- 0.5, 
      chi_F_1 <- chi_F_1_vals[h], 
      chi_F_12 <- chi_F_12_vals[h],
      chi_F_21 <- chi_F_21_vals[h],
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I1_1.eq[i] <- out[t_final+1, 2]
    I1_2.eq[i] <- out[t_final+1, 3]
    I2_2.eq[i] <- out[t_final+1, 4]
    I1_12.eq[i] <- out[t_final+1, 5]
    I1_21.eq[i] <- out[t_final+1, 6]
    P1.eq[i] <- out[t_final+1, 7]
    P2.eq[i] <- out[t_final+1, 8]
  }
  
  ## save propagule densities
  dd1 <- P1.eq
  record1 <- rbind(record1,dd1)
  
  dd2 <- P2.eq
  record2 <- rbind(record2,dd2)
}


colnames(record1)<- NULL 
colnames(record2)<- NULL


A = which(is.na(record1), arr.ind = TRUE)

B = which(is.na(record2), arr.ind = TRUE)

index = union(A, B)

if(length(index) == 0)
{
  
  ScenarioS5_P1 <- record1
  ScenarioS5_P2 <- record2
}else
{
  ScenarioS5_P1 <- record1[-index]
  ScenarioS5_P2 <- record2[-index]
}

##########################
## SCENARIO S6 
# Parameter values

chi_F_1_vals = 3000
chi_F_12_vals = 1500
chi_F_21_vals = 2000 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(chi_F_12_vals))
{
  I1_1.eq <- c()
  I1_2.eq <- c()
  I2_2.eq <- c()
  I1_12.eq <- c()
  I1_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      m_F_1 <- 0.5, 
      m_F_21 <- 0.5, 
      m_F_12 <- 0.5, 
      chi_F_1 <- chi_F_1_vals[h], 
      chi_F_12 <- chi_F_12_vals[h],
      chi_F_21 <- chi_F_21_vals[h],
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I1_1.eq[i] <- out[t_final+1, 2]
    I1_2.eq[i] <- out[t_final+1, 3]
    I2_2.eq[i] <- out[t_final+1, 4]
    I1_12.eq[i] <- out[t_final+1, 5]
    I1_21.eq[i] <- out[t_final+1, 6]
    P1.eq[i] <- out[t_final+1, 7]
    P2.eq[i] <- out[t_final+1, 8]
  }
  
  ## save propagule densities
  dd1 <- P1.eq
  record1 <- rbind(record1,dd1)
  
  dd2 <- P2.eq
  record2 <- rbind(record2,dd2)
}


colnames(record1)<- NULL 
colnames(record2)<- NULL

A = which(is.na(record1), arr.ind = TRUE)

B = which(is.na(record2), arr.ind = TRUE)

index = union(A, B)

if(length(index) == 0)
{
  
  ScenarioS6_P1 <- record1
  ScenarioS6_P2 <- record2
}else
{
  ScenarioS6_P1 <- record1[-index]
  ScenarioS6_P2 <- record2[-index]
}


#########################################################################
### Step 4: COMPUTE EQUILIBRIA FOR SCENARIOS M1-M3

## SCENARIO M1 
# Parameter values

m_F_1_vals = c(0.5, 0.5, 0.5)
m_F_12_vals = c(15, 10, 0.5)
m_F_21_vals = c(15, 15, 15) 


record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(m_F_12_vals))
{
  I1_1.eq <- c()
  I1_2.eq <- c()
  I2_2.eq <- c()
  I1_12.eq <- c()
  I1_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      m_F_1 <- m_F_1_vals[h], 
      m_F_21 <- m_F_21_vals[h], 
      m_F_12 <- m_F_12_vals[h], 
      chi_F_1 <- 2000, 
      chi_F_12 <- 2000,
      chi_F_21 <- 2000,
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I1_1.eq[i] <- out[t_final+1, 2]
    I1_2.eq[i] <- out[t_final+1, 3]
    I2_2.eq[i] <- out[t_final+1, 4]
    I1_12.eq[i] <- out[t_final+1, 5]
    I1_21.eq[i] <- out[t_final+1, 6]
    P1.eq[i] <- out[t_final+1, 7]
    P2.eq[i] <- out[t_final+1, 8]
  }
  
  ## save propagule densities
  dd1 <- P1.eq
  record1 <- rbind(record1,dd1)
  
  dd2 <- P2.eq
  record2 <- rbind(record2,dd2)
}


colnames(record1)<- NULL 
colnames(record2)<- NULL


A = which(is.na(record1), arr.ind = TRUE)

B = which(is.na(record2), arr.ind = TRUE)

index = union(A, B)

if(length(index) == 0)
{
  
  ScenarioM1_P1 <- record1
  ScenarioM1_P2 <- record2
}else
{
  ScenarioM1_P1 <- record1[-index]
  ScenarioM1_P2 <- record2[-index]
}

##########################
## SCENARIO M2 
# Parameter values

m_F_1_vals = c(10, 10, 10)
m_F_12_vals = c(10, 5, 0.5)  
m_F_21_vals = c(15, 15, 15) 


record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(m_F_12_vals))
{
  I1_1.eq <- c()
  I1_2.eq <- c()
  I2_2.eq <- c()
  I1_12.eq <- c()
  I1_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      m_F_1 <- m_F_1_vals[h], 
      m_F_21 <- m_F_21_vals[h], 
      m_F_12 <- m_F_12_vals[h], 
      chi_F_1 <- 2000, 
      chi_F_12 <- 2000,
      chi_F_21 <- 2000,
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I1_1.eq[i] <- out[t_final+1, 2]
    I1_2.eq[i] <- out[t_final+1, 3]
    I2_2.eq[i] <- out[t_final+1, 4]
    I1_12.eq[i] <- out[t_final+1, 5]
    I1_21.eq[i] <- out[t_final+1, 6]
    P1.eq[i] <- out[t_final+1, 7]
    P2.eq[i] <- out[t_final+1, 8]
  }
  
  ## save propagule densities
  dd1 <- P1.eq
  record1 <- rbind(record1,dd1)
  
  dd2 <- P2.eq
  record2 <- rbind(record2,dd2)
}


colnames(record1)<- NULL 
colnames(record2)<- NULL

A = which(is.na(record1), arr.ind = TRUE)

B = which(is.na(record2), arr.ind = TRUE)

index = union(A, B)

if(length(index) == 0)
{
  
  ScenarioM2_P1 <- record1
  ScenarioM2_P2 <- record2
}else
{
  ScenarioM2_P1 <- record1[,-index]
  ScenarioM2_P2 <- record2[,-index]
}


##########################
## SCENARIO M3 
# Parameter values

m_F_1_vals = c(15, 15, 15)
m_F_12_vals = c(8, 3, 0.5)
m_F_21_vals = c(10, 10, 10) 


record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(m_F_12_vals))
{
  I1_1.eq <- c()
  I1_2.eq <- c()
  I2_2.eq <- c()
  I1_12.eq <- c()
  I1_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      m_F_1 <- m_F_1_vals[h], 
      m_F_21 <- m_F_21_vals[h], 
      m_F_12 <- m_F_12_vals[h], 
      chi_F_1 <- 2000, 
      chi_F_12 <- 2000,
      chi_F_21 <- 2000,
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I1_1.eq[i] <- out[t_final+1, 2]
    I1_2.eq[i] <- out[t_final+1, 3]
    I2_2.eq[i] <- out[t_final+1, 4]
    I1_12.eq[i] <- out[t_final+1, 5]
    I1_21.eq[i] <- out[t_final+1, 6]
    P1.eq[i] <- out[t_final+1, 7]
    P2.eq[i] <- out[t_final+1, 8]
  }
  
  ## save propagule densities
  dd1 <- P1.eq
  record1 <- rbind(record1,dd1)
  
  dd2 <- P2.eq
  record2 <- rbind(record2,dd2)
}


colnames(record1)<- NULL 
colnames(record2)<- NULL

A = which(is.na(record1), arr.ind = TRUE)

B = which(is.na(record2), arr.ind = TRUE)

index = union(A, B)

if(length(index) == 0)
{
  
  ScenarioM3_P1 <- record1
  ScenarioM3_P2 <- record2
}else
{
  ScenarioM3_P1 <- record1[,-index]
  ScenarioM3_P2 <- record2[,-index]
}

#########################################################################
### Step 5: COMPUTE EQUILIBRIA FOR PAIRWISE COMBINATIONS OF S AND M SCENARIOS


##########################
## SCENARIO S4 and M1 
# Parameter values

chi_F_1_vals = 1500
chi_F_12_vals = 2000
chi_F_21_vals = 3000 

m_F_1_vals = c(0.5, 0.5, 0.5)
m_F_12_vals = c(15, 10, 0.5)
m_F_21_vals = c(15, 15, 15) 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(m_F_12_vals))
{
  I1_1.eq <- c()
  I1_2.eq <- c()
  I2_2.eq <- c()
  I1_12.eq <- c()
  I1_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      m_F_1 <- m_F_1_vals[h], 
      m_F_21 <- m_F_21_vals[h],
      m_F_12 <- m_F_12_vals[h],
      chi_F_1 <- chi_F_1_vals, 
      chi_F_12 <- chi_F_12_vals,
      chi_F_21 <- chi_F_21_vals,
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I1_1.eq[i] <- out[t_final+1, 2]
    I1_2.eq[i] <- out[t_final+1, 3]
    I2_2.eq[i] <- out[t_final+1, 4]
    I1_12.eq[i] <- out[t_final+1, 5]
    I1_21.eq[i] <- out[t_final+1, 6]
    P1.eq[i] <- out[t_final+1, 7]
    P2.eq[i] <- out[t_final+1, 8]
  }
  
  ## save propagule densities
  dd1 <- P1.eq
  record1 <- rbind(record1,dd1)
  
  dd2 <- P2.eq
  record2 <- rbind(record2,dd2)
}


colnames(record1)<- NULL 
colnames(record2)<- NULL

A = which(is.na(record1), arr.ind = TRUE)

B = which(is.na(record2), arr.ind = TRUE)

index = union(A, B)

if(length(index) == 0)
{
  
  ScenarioS4M1_P1 <- record1
  ScenarioS4M1_P2 <- record2
}else
{
  ScenarioS4M1_P1 <- record1[,-index]
  ScenarioS4M1_P2 <- record2[,-index]
}


##########################
## SCENARIO S4 and M2 
# Parameter values

chi_F_1_vals = 1500
chi_F_12_vals = 2000
chi_F_21_vals = 3000 

m_F_1_vals= c(10, 10, 10)
m_F_12_vals = c(10, 5, 0.5)  
m_F_21_vals = c(15, 15, 15) 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(m_F_12_vals))
{
  I1_1.eq <- c()
  I1_2.eq <- c()
  I2_2.eq <- c()
  I1_12.eq <- c()
  I1_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      m_F_1 <- m_F_1_vals[h], 
      m_F_21 <- m_F_21_vals[h], 
      m_F_12 <- m_F_12_vals[h], 
      chi_F_1 <- chi_F_1_vals, 
      chi_F_12 <- chi_F_12_vals,
      chi_F_21 <- chi_F_21_vals,
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I1_1.eq[i] <- out[t_final+1, 2]
    I1_2.eq[i] <- out[t_final+1, 3]
    I2_2.eq[i] <- out[t_final+1, 4]
    I1_12.eq[i] <- out[t_final+1, 5]
    I1_21.eq[i] <- out[t_final+1, 6]
    P1.eq[i] <- out[t_final+1, 7]
    P2.eq[i] <- out[t_final+1, 8]
  }
  
  ## save propagule densities
  dd1 <- P1.eq
  record1 <- rbind(record1,dd1)
  
  dd2 <- P2.eq
  record2 <- rbind(record2,dd2)
}


colnames(record1)<- NULL 
colnames(record2)<- NULL

A = which(is.na(record1), arr.ind = TRUE)

B = which(is.na(record2), arr.ind = TRUE)

index = union(A, B)

if(length(index) == 0)
{
  
  ScenarioS4M2_P1 <- record1
  ScenarioS4M2_P2 <- record2
}else
{
  ScenarioS4M2_P1 <- record1[,-index]
  ScenarioS4M2_P2 <- record2[,-index]
}

##########################
## SCENARIO S4 and M3 
# Parameter values

chi_F_1_vals = 1500
chi_F_12_vals = 2000
chi_F_21_vals = 3000 

m_F_1_vals = c(15, 15, 15)
m_F_12_vals = c(8, 3, 0.5)
m_F_21_vals = c(10, 10, 10) 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(m_F_12_vals))
{
  I1_1.eq <- c()
  I1_2.eq <- c()
  I2_2.eq <- c()
  I1_12.eq <- c()
  I1_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      m_F_1 <- m_F_1_vals[h], 
      m_F_21 <- m_F_21_vals[h], #fixed
      m_F_12 <- m_F_12_vals[h], #changing
      chi_F_1 <- chi_F_1_vals, 
      chi_F_12 <- chi_F_12_vals,#fixed
      chi_F_21 <- chi_F_21_vals,#changing
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I1_1.eq[i] <- out[t_final+1, 2]
    I1_2.eq[i] <- out[t_final+1, 3]
    I2_2.eq[i] <- out[t_final+1, 4]
    I1_12.eq[i] <- out[t_final+1, 5]
    I1_21.eq[i] <- out[t_final+1, 6]
    P1.eq[i] <- out[t_final+1, 7]
    P2.eq[i] <- out[t_final+1, 8]
  }
  
  ## save propagule densities
  dd1 <- P1.eq
  record1 <- rbind(record1,dd1)
  
  dd2 <- P2.eq
  record2 <- rbind(record2,dd2)
}


colnames(record1)<- NULL 
colnames(record2)<- NULL

A = which(is.na(record1), arr.ind = TRUE)

B = which(is.na(record2), arr.ind = TRUE)

index = union(A, B)

if(length(index) == 0)
{
  
  ScenarioS4M3_P1 <- record1
  ScenarioS4M3_P2 <- record2
}else
{
  ScenarioS4M3_P1 <- record1[,-index]
  ScenarioS4M3_P2 <- record2[,-index]
}

##########################
## SCENARIO S5 and M1 
# Parameter values

chi_F_1_vals = 2000
chi_F_12_vals = 1500
chi_F_21_vals = 3000 

m_F_1_vals = c(0.5, 0.5, 0.5)
m_F_12_vals = c(15, 10, 0.5)
m_F_21_vals = c(15, 15, 15) 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(m_F_12_vals))
{
  I1_1.eq <- c()
  I1_2.eq <- c()
  I2_2.eq <- c()
  I1_12.eq <- c()
  I1_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      m_F_1 <- m_F_1_vals[h], 
      m_F_21 <- m_F_21_vals[h], 
      m_F_12 <- m_F_12_vals[h], 
      chi_F_1 <- chi_F_1_vals, 
      chi_F_12 <- chi_F_12_vals,
      chi_F_21 <- chi_F_21_vals,
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I1_1.eq[i] <- out[t_final+1, 2]
    I1_2.eq[i] <- out[t_final+1, 3]
    I2_2.eq[i] <- out[t_final+1, 4]
    I1_12.eq[i] <- out[t_final+1, 5]
    I1_21.eq[i] <- out[t_final+1, 6]
    P1.eq[i] <- out[t_final+1, 7]
    P2.eq[i] <- out[t_final+1, 8]
  }
  
  ## save propagule densities
  dd1 <- P1.eq
  record1 <- rbind(record1,dd1)
  
  dd2 <- P2.eq
  record2 <- rbind(record2,dd2)
}


colnames(record1)<- NULL 
colnames(record2)<- NULL

A = which(is.na(record1), arr.ind = TRUE)

B = which(is.na(record2), arr.ind = TRUE)

index = union(A, B)

if(length(index) == 0)
{
  
  ScenarioS5M1_P1 <- record1
  ScenarioS5M1_P2 <- record2
}else
{
  ScenarioS5M1_P1 <- record1[,-index]
  ScenarioS5M1_P2 <- record2[,-index]
}

##########################
## SCENARIO S5 and M2 
# Parameter values

chi_F_1_vals = 2000
chi_F_12_vals = 1500
chi_F_21_vals = 3000 

m_F_1_vals = c(10, 10, 10)
m_F_12_vals = c(10, 5, 0.5)  
m_F_21_vals = c(15, 15, 15) 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(m_F_12_vals))
{
  I1_1.eq <- c()
  I1_2.eq <- c()
  I2_2.eq <- c()
  I1_12.eq <- c()
  I1_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      m_F_1 <- m_F_1_vals[h], 
      m_F_21 <- m_F_21_vals[h],
      m_F_12 <- m_F_12_vals[h],
      chi_F_1 <- chi_F_1_vals, 
      chi_F_12 <- chi_F_12_vals,
      chi_F_21 <- chi_F_21_vals,
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I1_1.eq[i] <- out[t_final+1, 2]
    I1_2.eq[i] <- out[t_final+1, 3]
    I2_2.eq[i] <- out[t_final+1, 4]
    I1_12.eq[i] <- out[t_final+1, 5]
    I1_21.eq[i] <- out[t_final+1, 6]
    P1.eq[i] <- out[t_final+1, 7]
    P2.eq[i] <- out[t_final+1, 8]
  }
  
  ## save propagule densities
  dd1 <- P1.eq
  record1 <- rbind(record1,dd1)
  
  dd2 <- P2.eq
  record2 <- rbind(record2,dd2)
}


colnames(record1)<- NULL 
colnames(record2)<- NULL

A = which(is.na(record1), arr.ind = TRUE)

B = which(is.na(record2), arr.ind = TRUE)

index = union(A, B)

if(length(index) == 0)
{
  
  ScenarioS5M2_P1 <- record1
  ScenarioS5M2_P2 <- record2
}else
{
  ScenarioS5M2_P1 <- record1[,-index]
  ScenarioS5M2_P2 <- record2[,-index]
}


##########################
## SCENARIO S5 and M3 
# Parameter values

chi_F_1_vals = 2000
chi_F_12_vals = 1500
chi_F_21_vals = 3000 

m_F_1_vals = c(15, 15, 15)
m_F_12_vals = c(8, 3, 0.5)
m_F_21_vals = c(10, 10, 10) 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(m_F_12_vals))
{
  I1_1.eq <- c()
  I1_2.eq <- c()
  I2_2.eq <- c()
  I1_12.eq <- c()
  I1_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      m_F_1 <- m_F_1_vals[h], 
      m_F_21 <- m_F_21_vals[h], 
      m_F_12 <- m_F_12_vals[h], 
      chi_F_1 <- chi_F_1_vals, 
      chi_F_12 <- chi_F_12_vals,
      chi_F_21 <- chi_F_21_vals,
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I1_1.eq[i] <- out[t_final+1, 2]
    I1_2.eq[i] <- out[t_final+1, 3]
    I2_2.eq[i] <- out[t_final+1, 4]
    I1_12.eq[i] <- out[t_final+1, 5]
    I1_21.eq[i] <- out[t_final+1, 6]
    P1.eq[i] <- out[t_final+1, 7]
    P2.eq[i] <- out[t_final+1, 8]
  }
  
  ## save propagule densities
  dd1 <- P1.eq
  record1 <- rbind(record1,dd1)
  
  dd2 <- P2.eq
  record2 <- rbind(record2,dd2)
}


colnames(record1)<- NULL 
colnames(record2)<- NULL

A = which(is.na(record1), arr.ind = TRUE)

B = which(is.na(record2), arr.ind = TRUE)

index = union(A, B)

if(length(index) == 0)
{
  
  ScenarioS5M3_P1 <- record1
  ScenarioS5M3_P2 <- record2
}else
{
  ScenarioS5M3_P1 <- record1[,-index]
  ScenarioS5M3_P2 <- record2[,-index]
}


##########################
## SCENARIO S6 and M1 
# Parameter values

chi_F_1_vals = 3000
chi_F_12_vals = 1500
chi_F_21_vals = 2000 

m_F_1_vals = c(0.5, 0.5, 0.5)
m_F_12_vals = c(15, 10, 0.5)
m_F_21_vals = c(15, 15, 15) 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(m_F_12_vals))
{
  I1_1.eq <- c()
  I1_2.eq <- c()
  I2_2.eq <- c()
  I1_12.eq <- c()
  I1_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      m_F_1 <- m_F_1_vals[h], 
      m_F_21 <- m_F_21_vals[h], 
      m_F_12 <- m_F_12_vals[h], 
      chi_F_1 <- chi_F_1_vals, 
      chi_F_12 <- chi_F_12_vals,
      chi_F_21 <- chi_F_21_vals,
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I1_1.eq[i] <- out[t_final+1, 2]
    I1_2.eq[i] <- out[t_final+1, 3]
    I2_2.eq[i] <- out[t_final+1, 4]
    I1_12.eq[i] <- out[t_final+1, 5]
    I1_21.eq[i] <- out[t_final+1, 6]
    P1.eq[i] <- out[t_final+1, 7]
    P2.eq[i] <- out[t_final+1, 8]
  }
  
  ## save propagule densities
  dd1 <- P1.eq
  record1 <- rbind(record1,dd1)
  
  dd2 <- P2.eq
  record2 <- rbind(record2,dd2)
}


colnames(record1)<- NULL 
colnames(record2)<- NULL

A = which(is.na(record1), arr.ind = TRUE)

B = which(is.na(record2), arr.ind = TRUE)

index = union(A, B)

if(length(index) == 0)
{
  
  ScenarioS6M1_P1 <- record1
  ScenarioS6M1_P2 <- record2
}else
{
  ScenarioS6M1_P1 <- record1[,-index]
  ScenarioS6M1_P2 <- record2[,-index]
}


##########################
## SCENARIO S6 and M2 
# Parameter values

chi_F_1_vals = 3000
chi_F_12_vals = 1500
chi_F_21_vals = 2000 

m_F_1_vals = c(10, 10, 10)
m_F_12_vals = c(10, 5, 0.5)  
m_F_21_vals = c(15, 15, 15)

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(m_F_12_vals))
{
  I1_1.eq <- c()
  I1_2.eq <- c()
  I2_2.eq <- c()
  I1_12.eq <- c()
  I1_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      m_F_1 <- m_F_1_vals[h], 
      m_F_21 <- m_F_21_vals[h], 
      m_F_12 <- m_F_12_vals[h], 
      chi_F_1 <- chi_F_1_vals, 
      chi_F_12 <- chi_F_12_vals,
      chi_F_21 <- chi_F_21_vals,
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I1_1.eq[i] <- out[t_final+1, 2]
    I1_2.eq[i] <- out[t_final+1, 3]
    I2_2.eq[i] <- out[t_final+1, 4]
    I1_12.eq[i] <- out[t_final+1, 5]
    I1_21.eq[i] <- out[t_final+1, 6]
    P1.eq[i] <- out[t_final+1, 7]
    P2.eq[i] <- out[t_final+1, 8]
  }
  
  ## save propagule densities
  dd1 <- P1.eq
  record1 <- rbind(record1,dd1)
  
  dd2 <- P2.eq
  record2 <- rbind(record2,dd2)
}


colnames(record1)<- NULL 
colnames(record2)<- NULL

A = which(is.na(record1), arr.ind = TRUE)

B = which(is.na(record2), arr.ind = TRUE)


index = union(A, B)

if(length(index) == 0)
{
  
  ScenarioS6M2_P1 <- record1
  ScenarioS6M2_P2 <- record2
}else
{
  ScenarioS6M2_P1 <- record1[,-index]
  ScenarioS6M2_P2 <- record2[,-index]
}


##########################
## SCENARIO S6 and M3 
# Parameter values

chi_F_1_vals = 3000
chi_F_12_vals = 1500
chi_F_21_vals = 2000 

m_F_1_vals = c(15, 15, 15)
m_F_12_vals = c(8, 3, 0.5)
m_F_21_vals = c(10, 10, 10) 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(m_F_12_vals))
{
  I1_1.eq <- c()
  I1_2.eq <- c()
  I2_2.eq <- c()
  I1_12.eq <- c()
  I1_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      m_F_1 <- m_F_1_vals[h], 
      m_F_21 <- m_F_21_vals[h], 
      m_F_12 <- m_F_12_vals[h], 
      chi_F_1 <- chi_F_1_vals, 
      chi_F_12 <- chi_F_12_vals,
      chi_F_21 <- chi_F_21_vals,
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I1_1.eq[i] <- out[t_final+1, 2]
    I1_2.eq[i] <- out[t_final+1, 3]
    I2_2.eq[i] <- out[t_final+1, 4]
    I1_12.eq[i] <- out[t_final+1, 5]
    I1_21.eq[i] <- out[t_final+1, 6]
    P1.eq[i] <- out[t_final+1, 7]
    P2.eq[i] <- out[t_final+1, 8]
  }
  
  ## save propagule densities
  dd1 <- P1.eq
  record1 <- rbind(record1,dd1)
  
  dd2 <- P2.eq
  record2 <- rbind(record2,dd2)
}


colnames(record1)<- NULL 
colnames(record2)<- NULL

A = which(is.na(record1), arr.ind = TRUE)

B = which(is.na(record2), arr.ind = TRUE)


index = union(A, B)

if(length(index) == 0)
{
  
  ScenarioS6M3_P1 <- record1
  ScenarioS6M3_P2 <- record2
}else
{
  ScenarioS6M3_P1 <- record1[,-index]
  ScenarioS6M3_P2 <- record2[,-index]
}


###########################################################################
# STEP 6: PLOT P1-P2 RELATIONSHIPS

par(mfrow = c(4, 4))

par(mar = c(1.85, 1.85, 1.3, 0.5))

plot.new()

plot(
  as.numeric(ScenarioS4_P2[1, 1:ncol(ScenarioS4_P2)]),
  ScenarioS4_P1[1, 1:ncol(ScenarioS4_P1)],
  lwd = 1,
  col = "red",
  type = "l",
  xlab = "",
  ylab = "",
  ylim = c(min(ScenarioS4_P1[1, 1:ncol(ScenarioS4_P1)]), (0 + max(ScenarioS4_P1[1, 1:ncol(ScenarioS4_P1)])))
)
legend(
  500,
  2800,
  cex = 0.8,
  legend = c("(1.5K, 2K, 3K)"),
  col = c("black"),
  bty = "n"
)

plot(
  as.numeric(ScenarioS5_P2[1, 1:ncol(ScenarioS5_P2)]),
  ScenarioS5_P1[1, 1:ncol(ScenarioS5_P1)],
  lwd = 1,
  col = "red",
  type = "l",
  xlab = "",
  ylab = "",
  ylim = c(min(ScenarioS5_P1[1, 1:ncol(ScenarioS5_P1)]), (0 + max(ScenarioS5_P1[1, 1:ncol(ScenarioS5_P1)])))
)
legend(
  1000,
  2800,
  cex = 0.8,
  legend = c("(2K, 1.5K, 3K)"),
  col = c("black"),
  bty = "n"
)

plot(
  as.numeric(ScenarioS6_P2[1, 1:ncol(ScenarioS6_P2)]),
  ScenarioS6_P1[1, 1:ncol(ScenarioS6_P1)],
  lwd = 1,
  col = "red",
  type = "l",
  xlab = "",
  ylab = "",
  ylim = c(min(ScenarioS6_P1[1, 1:ncol(ScenarioS6_P1)]), (0 + max(ScenarioS6_P1[1, 1:ncol(ScenarioS6_P1)])))
)
legend(
  500,
  3300,
  cex = 0.8,
  legend = c("(3K, 1.5K, 2K)"),
  col = c("black"),
  bty = "n"
)

plot(
  as.numeric(ScenarioM1_P2[1, seq(1, ncol(ScenarioM1_P2), 2)]),
  ScenarioM1_P1[1, seq(1, ncol(ScenarioM1_P1), 2)],
  lwd = 1,
  cex = 1.5,
  type = "o",
  pch = 1,
  col = "skyblue",
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioM1_P1[1, seq(1, ncol(ScenarioM1_P1), 2)], ScenarioM1_P1[2, seq(1, ncol(ScenarioM1_P1), 2)], ScenarioM1_P1[3, seq(1, ncol(ScenarioM1_P1), 2)]),
    (3500 + max(
      ScenarioM1_P1[1, seq(1, ncol(ScenarioM1_P1), 2)], ScenarioM1_P1[2, seq(1, ncol(ScenarioM1_P1), 2)], ScenarioM1_P1[3, seq(1, ncol(ScenarioM1_P1), 2)]
    ))
  )
)
points(
  as.numeric(ScenarioM1_P2[2, seq(1, ncol(ScenarioM1_P2), 2)]),
  ScenarioM1_P1[2, seq(1, ncol(ScenarioM1_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "steelblue",
  type = "o",
  pch = 2
)
points(
  as.numeric(ScenarioM1_P2[3, seq(1, ncol(ScenarioM1_P2), 2)]),
  ScenarioM1_P1[3, seq(1, ncol(ScenarioM1_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue4",
  type = "o",
  pch = 0
)
legend(
  700,
  (3000 + max(
    ScenarioM1_P1[1, seq(1, ncol(ScenarioM1_P1), 2)], ScenarioM1_P1[2, seq(1, ncol(ScenarioM1_P1), 2)], ScenarioM1_P1[3, seq(1, ncol(ScenarioM1_P1), 2)]
  )),
  cex = 0.9,
  legend = c("(0.5, 15, 15)", "(0.5, 10, 15)", "(0.5, 0.5, 15)"),
  col = c("skyblue", "steelblue", "blue4"),
  pch = c(1, 2, 0),
  bty = "n"
)

plot(
  as.numeric(ScenarioS4M1_P2[1, seq(1, ncol(ScenarioS4M1_P2), 2)]),
  ScenarioS4M1_P1[1, seq(1, ncol(ScenarioS4M1_P1), 2)],
  lwd = 1,
  cex = 1.5,
  type = "o",
  pch = 1,
  col = "red",
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioS4M1_P1[1, seq(1, ncol(ScenarioS4M1_P1), 2)], ScenarioS4M1_P1[2, seq(1, ncol(ScenarioS4M1_P1), 2)], ScenarioS4M1_P1[3, seq(1, ncol(ScenarioS4M1_P1), 2)]),
    (
      500 + max(ScenarioS4M1_P1[1, seq(1, ncol(ScenarioS4M1_P1), 2)], ScenarioS4M1_P1[2, seq(1, ncol(ScenarioS4M1_P1), 2)], ScenarioS4M1_P1[3, seq(1, ncol(ScenarioS4M1_P1), 2)])
    )
  )
)
points(
  as.numeric(ScenarioS4M1_P2[2, seq(1, ncol(ScenarioS4M1_P2), 2)]),
  ScenarioS4M1_P1[2, seq(1, ncol(ScenarioS4M1_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "black",
  type = "o",
  pch = 2
)
points(
  as.numeric(ScenarioS4M1_P2[3, seq(1, ncol(ScenarioS4M1_P2), 2)]),
  ScenarioS4M1_P1[3, seq(1, ncol(ScenarioS4M1_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 0
)

plot(
  as.numeric(ScenarioS5M1_P2[1, seq(1, ncol(ScenarioS5M1_P2), 2)]),
  ScenarioS5M1_P1[1, seq(1, ncol(ScenarioS5M1_P1), 2)],
  lwd = 1,
  cex = 1.5,
  type = "o",
  pch = 1,
  col = "red",
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioS5M1_P1[1, seq(1, ncol(ScenarioS5M1_P1), 2)], ScenarioS5M1_P1[2, seq(1, ncol(ScenarioS5M1_P1), 2)], ScenarioS5M1_P1[3, seq(1, ncol(ScenarioS5M1_P1), 2)]),
    (
      500 + max(ScenarioS5M1_P1[1, seq(1, ncol(ScenarioS5M1_P1), 2)], ScenarioS5M1_P1[2, seq(1, ncol(ScenarioS5M1_P1), 2)], ScenarioS5M1_P1[3, seq(1, ncol(ScenarioS5M1_P1), 2)])
    )
  )
)
points(
  as.numeric(ScenarioS5M1_P2[2, seq(1, ncol(ScenarioS5M1_P2), 2)]),
  ScenarioS5M1_P1[2, seq(1, ncol(ScenarioS5M1_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "red",
  type = "o",
  pch = 2
)
points(
  as.numeric(ScenarioS5M1_P2[3, seq(1, ncol(ScenarioS5M1_P2), 2)]),
  ScenarioS5M1_P1[3, seq(1, ncol(ScenarioS5M1_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 0
)

plot(
  as.numeric(ScenarioS6M1_P2[1, seq(1, ncol(ScenarioS6M1_P2), 2)]),
  ScenarioS6M1_P1[1, seq(1, ncol(ScenarioS6M1_P1), 2)],
  lwd = 1,
  cex = 1.5,
  type = "o",
  pch = 1,
  col = "red",
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioS6M1_P1[1, seq(1, ncol(ScenarioS6M1_P1), 2)], ScenarioS6M1_P1[2, seq(1, ncol(ScenarioS6M1_P1), 2)], ScenarioS6M1_P1[3, seq(1, ncol(ScenarioS6M1_P1), 2)]),
    (
      500 + max(ScenarioS6M1_P1[1, seq(1, ncol(ScenarioS6M1_P1), 2)], ScenarioS6M1_P1[2, seq(1, ncol(ScenarioS6M1_P1), 2)], ScenarioS6M1_P1[3, seq(1, ncol(ScenarioS6M1_P1), 2)])
    )
  )
)
points(
  as.numeric(ScenarioS6M1_P2[2, seq(1, ncol(ScenarioS6M1_P2), 2)]),
  ScenarioS6M1_P1[2, seq(1, ncol(ScenarioS6M1_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 2
)
points(
  as.numeric(ScenarioS6M1_P2[3, seq(1, ncol(ScenarioS6M1_P2), 2)]),
  ScenarioS6M1_P1[3, seq(1, ncol(ScenarioS6M1_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 0
)

plot(
  as.numeric(ScenarioM2_P2[1, seq(1, ncol(ScenarioM2_P2), 2)]),
  ScenarioM2_P1[1, seq(1, ncol(ScenarioM2_P1), 2)],
  lwd = 1,
  cex = 1.5,
  type = "o",
  pch = 1,
  col = "skyblue",
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioM2_P1[1, seq(1, ncol(ScenarioM2_P1), 2)], ScenarioM2_P1[2, seq(1, ncol(ScenarioM2_P1), 2)], ScenarioM2_P1[3, seq(1, ncol(ScenarioM2_P1), 2)]),
    (3500 + max(
      ScenarioM2_P1[1, seq(1, ncol(ScenarioM2_P1), 2)], ScenarioM2_P1[2, seq(1, ncol(ScenarioM2_P1), 2)], ScenarioM2_P1[3, seq(1, ncol(ScenarioM2_P1), 2)]
    ))
  )
)
points(
  as.numeric(ScenarioM2_P2[2, seq(1, ncol(ScenarioM2_P2), 2)]),
  ScenarioM2_P1[2, seq(1, ncol(ScenarioM2_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "steelblue",
  type = "o",
  pch = 2
)
points(
  as.numeric(ScenarioM2_P2[3, seq(1, ncol(ScenarioM2_P2), 2)]),
  ScenarioM2_P1[3, seq(1, ncol(ScenarioM2_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue4",
  type = "o",
  pch = 0
)
legend(
  600,
  (3000 + max(
    ScenarioM2_P1[1, seq(1, ncol(ScenarioM2_P1), 2)], ScenarioM2_P1[2, seq(1, ncol(ScenarioM2_P1), 2)], ScenarioM2_P1[3, seq(1, ncol(ScenarioM2_P1), 2)]
  )),
  cex = 0.9,
  legend = c("(10, 10, 15)", "(10, 5, 15)", "(10, 0.5, 15)"),
  col = c("skyblue", "steelblue", "blue4"),
  pch = c(1, 2, 0),
  bty = "n"
)

plot(
  as.numeric(ScenarioS4M2_P2[1, seq(1, ncol(ScenarioS4M2_P2), 2)]),
  ScenarioS4M2_P1[1, seq(1, ncol(ScenarioS4M2_P1), 2)],
  lwd = 1,
  cex = 1.5,
  type = "o",
  pch = 1,
  col = "red",
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioS4M2_P1[1, seq(1, ncol(ScenarioS4M2_P1), 2)], ScenarioS4M2_P1[2, seq(1, ncol(ScenarioS4M2_P1), 2)], ScenarioS4M2_P1[3, seq(1, ncol(ScenarioS4M2_P1), 2)]),
    (
      500 + max(ScenarioS4M2_P1[1, seq(1, ncol(ScenarioS4M2_P1), 2)], ScenarioS4M2_P1[2, seq(1, ncol(ScenarioS4M2_P1), 2)], ScenarioS4M2_P1[3, seq(1, ncol(ScenarioS4M2_P1), 2)])
    )
  )
)
points(
  as.numeric(ScenarioS4M2_P2[2, seq(1, ncol(ScenarioS4M2_P2), 2)]),
  ScenarioS4M2_P1[2, seq(1, ncol(ScenarioS4M2_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 2
)
points(
  as.numeric(ScenarioS4M2_P2[3, seq(1, ncol(ScenarioS4M2_P2), 2)]),
  ScenarioS4M2_P1[3, seq(1, ncol(ScenarioS4M2_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 0
)

plot(
  as.numeric(ScenarioS5M2_P2[1, seq(1, ncol(ScenarioS5M2_P2), 2)]),
  ScenarioS5M2_P1[1, seq(1, ncol(ScenarioS5M2_P1), 2)],
  lwd = 1,
  cex = 1.5,
  type = "o",
  pch = 1,
  col = "red",
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioS5M2_P1[1, seq(1, ncol(ScenarioS5M2_P1), 2)], ScenarioS5M2_P1[2, seq(1, ncol(ScenarioS5M2_P1), 2)], ScenarioS5M2_P1[3, seq(1, ncol(ScenarioS5M2_P1), 2)]),
    (
      500 + max(ScenarioS5M2_P1[1, seq(1, ncol(ScenarioS5M2_P1), 2)], ScenarioS5M2_P1[2, seq(1, ncol(ScenarioS5M2_P1), 2)], ScenarioS5M2_P1[3, seq(1, ncol(ScenarioS5M2_P1), 2)])
    )
  )
)
points(
  as.numeric(ScenarioS5M2_P2[2, seq(1, ncol(ScenarioS5M2_P2), 2)]),
  ScenarioS5M2_P1[2, seq(1, ncol(ScenarioS5M2_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 2
)
points(
  as.numeric(ScenarioS5M2_P2[3, seq(1, ncol(ScenarioS5M2_P2), 2)]),
  ScenarioS5M2_P1[3, seq(1, ncol(ScenarioS5M2_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 0
)

plot(
  as.numeric(ScenarioS6M2_P2[1, seq(1, ncol(ScenarioS6M2_P2), 2)]),
  ScenarioS6M2_P1[1, seq(1, ncol(ScenarioS6M2_P1), 2)],
  lwd = 1,
  cex = 1.5,
  type = "o",
  pch = 1,
  col = "blue",
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioS6M2_P1[1, seq(1, ncol(ScenarioS6M2_P1), 2)], ScenarioS6M2_P1[2, seq(1, ncol(ScenarioS6M2_P1), 2)], ScenarioS6M2_P1[3, seq(1, ncol(ScenarioS6M2_P1), 2)]),
    (
      500 + max(ScenarioS6M2_P1[1, seq(1, ncol(ScenarioS6M2_P1), 2)], ScenarioS6M2_P1[2, seq(1, ncol(ScenarioS6M2_P1), 2)], ScenarioS6M2_P1[3, seq(1, ncol(ScenarioS6M2_P1), 2)])
    )
  )
)
points(
  as.numeric(ScenarioS6M2_P2[2, seq(1, ncol(ScenarioS6M2_P2), 2)]),
  ScenarioS6M2_P1[2, seq(1, ncol(ScenarioS6M2_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 2
)
points(
  as.numeric(ScenarioS6M2_P2[3, seq(1, ncol(ScenarioS6M2_P2), 2)]),
  ScenarioS6M2_P1[3, seq(1, ncol(ScenarioS6M2_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 0
)

plot(
  as.numeric(ScenarioM3_P2[1, seq(1, ncol(ScenarioM3_P2), 2)]),
  ScenarioM3_P1[1, seq(1, ncol(ScenarioM3_P1), 2)],
  lwd = 1,
  cex = 1.5,
  type = "o",
  pch = 1,
  col = "skyblue",
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioM3_P1[1, seq(1, ncol(ScenarioM3_P1), 2)], ScenarioM3_P1[2, seq(1, ncol(ScenarioM3_P1), 2)], ScenarioM3_P1[3, seq(1, ncol(ScenarioM3_P1), 2)]),
    (3500 + max(
      ScenarioM3_P1[1, seq(1, ncol(ScenarioM3_P1), 2)], ScenarioM3_P1[2, seq(1, ncol(ScenarioM3_P1), 2)], ScenarioM3_P1[3, seq(1, ncol(ScenarioM3_P1), 2)]
    ))
  )
)
points(
  as.numeric(ScenarioM3_P2[2, seq(1, ncol(ScenarioM3_P2), 2)]),
  ScenarioM3_P1[2, seq(1, ncol(ScenarioM3_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "steelblue",
  type = "o",
  pch = 2
)
points(
  as.numeric(ScenarioM3_P2[3, seq(1, ncol(ScenarioM3_P2), 2)]),
  ScenarioM3_P1[3, seq(1, ncol(ScenarioM3_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue4",
  type = "o",
  pch = 0
)
legend(
  700,
  (3000 + max(
    ScenarioM3_P1[1, seq(1, ncol(ScenarioM3_P1), 2)], ScenarioM3_P1[2, seq(1, ncol(ScenarioM3_P1), 2)], ScenarioM3_P1[3, seq(1, ncol(ScenarioM3_P1), 2)]
  )),
  cex = 0.9,
  legend = c("(15, 8, 10)", "(15, 3, 10)", "(15, 0.5, 10)"),
  col = c("skyblue", "steelblue", "blue4"),
  pch = c(1, 2, 0),
  bty = "n"
)

plot(
  as.numeric(ScenarioS4M3_P2[1, seq(1, ncol(ScenarioS4M3_P2), 2)]),
  ScenarioS4M3_P1[1, seq(1, ncol(ScenarioS4M3_P1), 2)],
  lwd = 1,
  cex = 1.5,
  type = "o",
  pch = 1,
  col = "red",
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioS4M3_P1[1, 1:ncol(ScenarioS4M3_P1)], ScenarioS4M3_P1[2, 1:ncol(ScenarioS4M3_P1)], ScenarioS4M3_P1[3, 1:ncol(ScenarioS4M3_P1)]),
    (
      500 + max(ScenarioS4M3_P1[1, 1:ncol(ScenarioS4M3_P1)], ScenarioS4M3_P1[2, 1:ncol(ScenarioS4M3_P1)], ScenarioS4M3_P1[3, 1:ncol(ScenarioS4M3_P1)])
    )
  )
)
points(
  as.numeric(ScenarioS4M3_P2[2, seq(1, ncol(ScenarioS4M3_P2), 2)]),
  ScenarioS4M3_P1[2, seq(1, ncol(ScenarioS4M3_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 2
)
points(
  as.numeric(ScenarioS4M3_P2[3, seq(1, ncol(ScenarioS4M3_P2), 2)]),
  ScenarioS4M3_P1[3, seq(1, ncol(ScenarioS4M3_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 0
)

plot(
  as.numeric(ScenarioS5M3_P2[1, seq(1, ncol(ScenarioS5M3_P2), 2)]),
  ScenarioS5M3_P1[1, seq(1, ncol(ScenarioS5M3_P1), 2)],
  lwd = 1,
  cex = 1.5,
  type = "o",
  pch = 1,
  col = "red",
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioS5M3_P1[1, seq(1, ncol(ScenarioS5M3_P1), 2)], ScenarioS5M3_P1[2, seq(1, ncol(ScenarioS5M3_P1), 2)], ScenarioS5M3_P1[3, seq(1, ncol(ScenarioS5M3_P1), 2)]),
    (
      500 + max(ScenarioS5M3_P1[1, seq(1, ncol(ScenarioS5M3_P1), 2)], ScenarioS5M3_P1[2, seq(1, ncol(ScenarioS5M3_P1), 2)], ScenarioS5M3_P1[3, seq(1, ncol(ScenarioS5M3_P1), 2)])
    )
  )
)
points(
  as.numeric(ScenarioS5M3_P2[2, seq(1, ncol(ScenarioS5M3_P2), 2)]),
  ScenarioS5M3_P1[2, seq(1, ncol(ScenarioS5M3_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 2
)
points(
  as.numeric(ScenarioS5M3_P2[3, seq(1, ncol(ScenarioS5M3_P2), 2)]),
  ScenarioS5M3_P1[3, seq(1, ncol(ScenarioS5M3_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 0
)

plot(
  as.numeric(ScenarioS6M3_P2[1, seq(1, ncol(ScenarioS6M3_P2), 2)]),
  ScenarioS6M3_P1[1, seq(1, ncol(ScenarioS6M3_P1), 2)],
  lwd = 1,
  cex = 1.5,
  type = "o",
  pch = 1,
  col = "red",
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioS6M3_P1[1, seq(1, ncol(ScenarioS6M3_P1), 2)], ScenarioS6M3_P1[2, seq(1, ncol(ScenarioS6M3_P1), 2)], ScenarioS6M3_P1[3, seq(1, ncol(ScenarioS6M3_P1), 2)]),
    (
      500 + max(ScenarioS6M3_P1[1, seq(1, ncol(ScenarioS6M3_P1), 2)], ScenarioS6M3_P1[2, seq(1, ncol(ScenarioS6M3_P1), 2)], ScenarioS6M3_P1[3, seq(1, ncol(ScenarioS6M3_P1), 2)])
    )
  )
)
points(
  as.numeric(ScenarioS6M3_P2[2, seq(1, ncol(ScenarioS6M3_P2), 2)]),
  ScenarioS6M3_P1[2, seq(1, ncol(ScenarioS6M3_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 2
)
points(
  as.numeric(ScenarioS6M3_P2[3, seq(1, ncol(ScenarioS6M3_P2), 2)]),
  ScenarioS6M3_P1[3, seq(1, ncol(ScenarioS6M3_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue",
  type = "o",
  pch = 0
)
