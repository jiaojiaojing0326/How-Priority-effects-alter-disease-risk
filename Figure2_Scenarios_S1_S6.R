# AUTHOR:   JING JIAO & MICHAEL H. CORTEZ
# DATE:     February 2022
# PURPOSE:  Simulate relationships between specialist and generalist infectious propagule
#           densities for scenarios S1-S6

# INPUTS:   Parameters for each host and pathogen

# METHOD:   Step 1: Define equations for ODE model
#           Step 2: Define global simulation parameters
#           Step 3: Compute equilibria for each scenario
#           Step 4: Plot P1-P2 relationships

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

# Parameters shared across simulations
parameters_shared<-c(
      N_F = 200,
      N_A = 200,

      beta_F_2 <- 0.01, 
      beta_F_1 <- 0.01, 
      beta_F_12 <- 0.01,
      beta_F_21 <- 0.01,
      beta_A_2 <- 0.01,

      m_F_1 <- 0.5, 
      m_F_2 <- 0.5,  
      m_F_12 <- 0.5,
      m_F_21 <- 0.5, 
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


##########################
## SCENARIO S1 
# Parameter values
chi_F_1_vals = c(3000, 3000, 3000)
chi_F_12_vals = c(2000, 2000, 2000)
chi_F_21_vals = c(2000, 1000, 500) 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(chi_F_12_vals))
{
  I_F_1.eq <- c()
  I_F_2.eq <- c()
  I_A_2.eq <- c()
  I_F_12.eq <- c()
  I_F_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      chi_F_1 <- chi_F_1_vals[h], 
      chi_F_12 <- chi_F_12_vals[h],
      chi_F_21 <- chi_F_21_vals[h],
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I_F_1.eq[i] <- out[t_final+1, 2]
    I_F_2.eq[i] <- out[t_final+1, 3]
    I_A_2.eq[i] <- out[t_final+1, 4]
    I_F_12.eq[i] <- out[t_final+1, 5]
    I_F_21.eq[i] <- out[t_final+1, 6]
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

ScenarioS1_P1 <- record1
ScenarioS1_P2<-record2


##########################
## SCENARIO S2 
# Parameter values
chi_F_1_vals = c(2000, 2000, 2000)
chi_F_12_vals = c(2000, 3000, 3500)
chi_F_21_vals = c(1500, 1500, 1500) 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(chi_F_12_vals))
{
  I_F_1.eq <- c()
  I_F_2.eq <- c()
  I_A_2.eq <- c()
  I_F_12.eq <- c()
  I_F_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      chi_F_1 <- chi_F_1_vals[h], 
      chi_F_12 <- chi_F_12_vals[h],
      chi_F_21 <- chi_F_21_vals[h],
      delta2 <- delta2_vals[i]
    )    
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I_F_1.eq[i] <- out[t_final+1, 2]
    I_F_2.eq[i] <- out[t_final+1, 3]
    I_A_2.eq[i] <- out[t_final+1, 4]
    I_F_12.eq[i] <- out[t_final+1, 5]
    I_F_21.eq[i] <- out[t_final+1, 6]
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

ScenarioS2_P1 <- record1
ScenarioS2_P2<-record2

##########################
## SCENARIO S3 
# Parameter values
chi_F_1_vals = c(1500, 1500, 1500)
chi_F_12_vals = c(2500, 3000, 3500)
chi_F_21_vals = c(2000, 2000, 2000) 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(chi_F_12_vals))
{
  I_F_1.eq <- c()
  I_F_2.eq <- c()
  I_A_2.eq <- c()
  I_F_12.eq <- c()
  I_F_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      chi_F_1 <- chi_F_1_vals[h], 
      chi_F_12 <- chi_F_12_vals[h],
      chi_F_21 <- chi_F_21_vals[h],
      delta2 <- delta2_vals[i]
    )    
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I_F_1.eq[i] <- out[t_final+1, 2]
    I_F_2.eq[i] <- out[t_final+1, 3]
    I_A_2.eq[i] <- out[t_final+1, 4]
    I_F_12.eq[i] <- out[t_final+1, 5]
    I_F_21.eq[i] <- out[t_final+1, 6]
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

ScenarioS3_P1 <- record1
ScenarioS3_P2<-record2

##########################
## SCENARIO S4 
# Parameter values
chi_F_1_vals = c(1500, 1500, 1500)
chi_F_12_vals = c(3000, 2500, 2000)
chi_F_21_vals = c(3000, 3000, 3000) 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(chi_F_12_vals))
{
  I_F_1.eq <- c()
  I_F_2.eq <- c()
  I_A_2.eq <- c()
  I_F_12.eq <- c()
  I_F_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      chi_F_1 <- chi_F_1_vals[h], 
      chi_F_12 <- chi_F_12_vals[h],
      chi_F_21 <- chi_F_21_vals[h],
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I_F_1.eq[i] <- out[t_final+1, 2]
    I_F_2.eq[i] <- out[t_final+1, 3]
    I_A_2.eq[i] <- out[t_final+1, 4]
    I_F_12.eq[i] <- out[t_final+1, 5]
    I_F_21.eq[i] <- out[t_final+1, 6]
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

ScenarioS4_P1 <- record1
ScenarioS4_P2<-record2

##########################
## SCENARIO S5 
# Parameter values
chi_F_1_vals = c(2000, 2000, 2000)
chi_F_12_vals = c(1500, 1500, 1500)
chi_F_21_vals = c(2000, 2500, 3000) 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(chi_F_12_vals))
{
  I_F_1.eq <- c()
  I_F_2.eq <- c()
  I_A_2.eq <- c()
  I_F_12.eq <- c()
  I_F_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      chi_F_1 <- chi_F_1_vals[h], 
      chi_F_12 <- chi_F_12_vals[h],
      chi_F_21 <- chi_F_21_vals[h],
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I_F_1.eq[i] <- out[t_final+1, 2]
    I_F_2.eq[i] <- out[t_final+1, 3]
    I_A_2.eq[i] <- out[t_final+1, 4]
    I_F_12.eq[i] <- out[t_final+1, 5]
    I_F_21.eq[i] <- out[t_final+1, 6]
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

ScenarioS5_P1 <- record1
ScenarioS5_P2<-record2

##########################
## SCENARIO S6 
# Parameter values
chi_F_1_vals = c(3000, 3000, 3000)
chi_F_12_vals = c(1500, 1500, 1500)
chi_F_21_vals = c(2000, 2500, 3000) 

record1 <- data.frame(NULL)
record2 <- data.frame(NULL)

for(h in 1:length(chi_F_12_vals))
{
  I_F_1.eq <- c()
  I_F_2.eq <- c()
  I_A_2.eq <- c()
  I_F_12.eq <- c()
  I_F_21.eq <- c()
  P1.eq <- c()
  P2.eq <- c()
  
  for(i in 1:length(delta2_vals))
  {
    parameters <- c(
      parameters_shared,
      chi_F_1 <- chi_F_1_vals[h], 
      chi_F_12 <- chi_F_12_vals[h],
      chi_F_21 <- chi_F_21_vals[h],
      delta2 <- delta2_vals[i]
    )   
    out <- ode(y = initial_conditions, times = times, func = DiffEqns, parms = parameters)
    
    I_F_1.eq[i] <- out[t_final+1, 2]
    I_F_2.eq[i] <- out[t_final+1, 3]
    I_A_2.eq[i] <- out[t_final+1, 4]
    I_F_12.eq[i] <- out[t_final+1, 5]
    I_F_21.eq[i] <- out[t_final+1, 6]
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

ScenarioS6_P1 <- record1
ScenarioS6_P2<-record2

###########################################################################
# STEP 4: PLOT P1-P2 RELATIONSHIPS


par(mfrow = c(3, 2))

par(mar = c(2, 2, 2, 2))

plot(
  as.numeric(ScenarioS1_P2[1, seq(1, ncol(ScenarioS1_P2), 2)]),
  ScenarioS1_P1[1, seq(1, ncol(ScenarioS1_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "skyblue",
  type = "o",
  pch = 1,
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioS1_P1[1, seq(1, ncol(ScenarioS1_P1), 2)], ScenarioS1_P1[2, seq(1, ncol(ScenarioS1_P1), 2)], ScenarioS1_P1[3, seq(1, ncol(ScenarioS1_P1), 2)]),
    (1200 + max(
      ScenarioS1_P1[1, seq(1, ncol(ScenarioS1_P1), 2)], ScenarioS1_P1[2, seq(1, ncol(ScenarioS1_P1), 2)], ScenarioS1_P1[3, seq(1, ncol(ScenarioS1_P1), 2)]
    ))
  )
)
points(
  as.numeric(ScenarioS1_P2[2, seq(1, ncol(ScenarioS1_P2), 2)]),
  ScenarioS1_P1[2, seq(1, ncol(ScenarioS1_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "steelblue",
  type = "o",
  pch = 2
)
points(
  as.numeric(ScenarioS1_P2[3, seq(1, ncol(ScenarioS1_P2), 2)]),
  ScenarioS1_P1[3, seq(1, ncol(ScenarioS1_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue4",
  type = "o",
  pch = 0
)
legend(
  1500,
  (950 + max(
    ScenarioS1_P1[1, seq(1, ncol(ScenarioS1_P1), 2)], ScenarioS1_P1[2, seq(1, ncol(ScenarioS1_P1), 2)], ScenarioS1_P1[3, seq(1, ncol(ScenarioS1_P1), 2)]
  )),
  legend = c("(3K, 2K, 2K)", "(3K, 2K, 1K)", "(3K, 2K, 0.5K)"),
  col = c("skyblue", "steelblue", "blue4"),
  pch = c(1, 2, 0),
  cex = c(1.2, 1.2, 1.2),
  bty = "n"
)


plot(
  as.numeric(ScenarioS4_P2[1, seq(1, ncol(ScenarioS4_P2), 2)]),
  ScenarioS4_P1[1, seq(1, ncol(ScenarioS4_P1), 2)],
  lwd = 1,
  cex = 1.5,
  type = "o",
  pch = 1,
  col = "tan2",
  xlab = "",
  ylab = "",
  ylim = c((
    -500 + min(ScenarioS4_P1[1, seq(1, ncol(ScenarioS4_P1), 2)], ScenarioS4_P1[2, seq(1, ncol(ScenarioS4_P1), 2)], ScenarioS4_P1[3, seq(1, ncol(ScenarioS4_P1), 2)])
  ), (
    100 + max(ScenarioS4_P1[1, seq(1, ncol(ScenarioS4_P1), 2)], ScenarioS4_P1[2, seq(1, ncol(ScenarioS4_P1), 2)], ScenarioS4_P1[3, seq(1, ncol(ScenarioS4_P1), 2)])
  ))
)
points(
  as.numeric(ScenarioS4_P2[2, seq(1, ncol(ScenarioS4_P2), 2)]),
  ScenarioS4_P1[2, seq(1, ncol(ScenarioS4_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "coral2",
  pch = 2,
  type = "o"
)
points(
  as.numeric(ScenarioS4_P2[3, seq(1, ncol(ScenarioS4_P2), 2)]),
  ScenarioS4_P1[3, seq(1, ncol(ScenarioS4_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "red3",
  pch = 0,
  type = "o"
)
legend(
  1000,
  (-1300 + max(
    ScenarioS4_P1[1, seq(1, ncol(ScenarioS4_P1), 2)], ScenarioS4_P1[2, seq(1, ncol(ScenarioS4_P1), 2)], ScenarioS4_P1[3, seq(1, ncol(ScenarioS4_P1), 2)]
  )),
  legend = c("(1.5K, 3K, 3K)", "(1.5K, 2.5K, 3K)", "(1.5K, 2K, 3K)"),
  col = c("tan2", "coral2", "red3"),
  pch = c(1, 2, 0),
  cex = c(1.2, 1.2, 1.2),
  bty = "n"
)


plot(
  as.numeric(ScenarioS2_P2[1, seq(1, ncol(ScenarioS2_P2), 2)]),
  ScenarioS2_P1[1, seq(1, ncol(ScenarioS2_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "skyblue",
  type = "o",
  pch = 1,
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioS2_P1[1, seq(1, ncol(ScenarioS2_P1), 2)], ScenarioS2_P1[2, seq(1, ncol(ScenarioS2_P1), 2)], ScenarioS2_P1[3, seq(1, ncol(ScenarioS2_P1), 2)]),
    (2000 + max(
      ScenarioS2_P1[1, seq(1, ncol(ScenarioS2_P1), 2)], ScenarioS2_P1[2, seq(1, ncol(ScenarioS2_P1), 2)], ScenarioS2_P1[3, seq(1, ncol(ScenarioS2_P1), 2)]
    ))
  )
)
points(
  as.numeric(ScenarioS2_P2[2, seq(1, ncol(ScenarioS2_P2), 2)]),
  ScenarioS2_P1[2, seq(1, ncol(ScenarioS2_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "steelblue",
  pch = 2,
  type = "o"
)
points(
  as.numeric(ScenarioS2_P2[3, seq(1, ncol(ScenarioS2_P2), 2)]),
  ScenarioS2_P1[3, seq(1, ncol(ScenarioS2_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue4",
  pch = 0,
  type = "o"
)
legend(
  1800,
  (1900 + max(
    ScenarioS2_P1[1, seq(1, ncol(ScenarioS2_P1), 2)], ScenarioS2_P1[2, seq(1, ncol(ScenarioS2_P1), 2)], ScenarioS2_P1[3, seq(1, ncol(ScenarioS2_P1), 2)]
  )),
  legend = c("(2K, 2K, 1.5K)", "(2K, 3K, 1.5K)", "(2K, 3.5K, 1.5K)"),
  col = c("skyblue", "steelblue", "blue4"),
  pch = c(1, 2, 0),
  cex = c(1.2, 1.2, 1.2),
  bty = "n"
)

plot(
  as.numeric(ScenarioS5_P2[1, seq(1, ncol(ScenarioS5_P2), 2)]),
  ScenarioS5_P1[1, seq(1, ncol(ScenarioS5_P1), 2)],
  lwd = 1,
  cex = 1.5,
  type = "o",
  pch = 1,
  col = "tan2",
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioS5_P1[1, seq(1, ncol(ScenarioS5_P1), 2)], ScenarioS5_P1[2, seq(1, ncol(ScenarioS5_P1), 2)], ScenarioS5_P1[3, seq(1, ncol(ScenarioS5_P1), 2)]),
    (1600 + max(
      ScenarioS5_P1[1, seq(1, ncol(ScenarioS5_P1), 2)], ScenarioS5_P1[2, seq(1, ncol(ScenarioS5_P1), 2)], ScenarioS5_P1[3, seq(1, ncol(ScenarioS5_P1), 2)]
    ))
  )
)
points(
  as.numeric(ScenarioS5_P2[2, seq(1, ncol(ScenarioS5_P2), 2)]),
  ScenarioS5_P1[2, seq(1, ncol(ScenarioS5_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "coral2",
  pch = 2,
  type = "o"
)
points(
  as.numeric(ScenarioS5_P2[3, seq(1, ncol(ScenarioS5_P2), 2)]),
  ScenarioS5_P1[3, seq(1, ncol(ScenarioS5_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "red3",
  pch = 0,
  type = "o"
)
legend(
  1000,
  (1500 + max(
    ScenarioS5_P1[1, seq(1, ncol(ScenarioS5_P1), 2)], ScenarioS5_P1[2, seq(1, ncol(ScenarioS5_P1), 2)], ScenarioS5_P1[3, seq(1, ncol(ScenarioS5_P1), 2)]
  )),
  legend = c("(2K, 1.5K, 2K)", "(2K, 1.5K, 2.5K)", "(2K, 1.5K, 3K)"),
  col = c("tan2", "coral2", "red3"),
  pch = c(1, 2, 0),
  cex = c(1.2, 1.2, 1.2),
  bty = "n"
)


plot(
  as.numeric(ScenarioS3_P2[1, seq(1, ncol(ScenarioS3_P2), 2)]),
  ScenarioS3_P1[1, seq(1, ncol(ScenarioS3_P1), 2)],
  lwd = 1,
  cex = 1.5,
  type = "o",
  col = "skyblue",
  pch = 1,
  xlab = "",
  ylab = "",
  ylim = c((
    -1500 + min(ScenarioS3_P1[1, seq(1, ncol(ScenarioS3_P1), 2)], ScenarioS3_P1[2, seq(1, ncol(ScenarioS3_P1), 2)], ScenarioS3_P1[3, seq(1, ncol(ScenarioS3_P1), 2)])
  ), (
    350 + max(ScenarioS3_P1[1, seq(1, ncol(ScenarioS3_P1), 2)], ScenarioS3_P1[2, seq(1, ncol(ScenarioS3_P1), 2)], ScenarioS3_P1[3, seq(1, ncol(ScenarioS3_P1), 2)])
  ))
)
points(
  as.numeric(ScenarioS3_P2[2, seq(1, ncol(ScenarioS3_P2), 2)]),
  ScenarioS3_P1[2, seq(1, ncol(ScenarioS3_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "steelblue",
  pch = 2,
  type = "o"
)
points(
  as.numeric(ScenarioS3_P2[3, seq(1, ncol(ScenarioS3_P2), 2)]),
  ScenarioS3_P1[3, seq(1, ncol(ScenarioS3_P1), 2)],
  lwd = 1,
  cex = 1.5,
  col = "blue4",
  type = "o",
  pch = 0
)
legend(
  1000,
  (-1800 + max(
    ScenarioS3_P1[1, seq(1, ncol(ScenarioS3_P1), 2)], ScenarioS3_P1[2, seq(1, ncol(ScenarioS3_P1), 2)], ScenarioS3_P1[3, seq(1, ncol(ScenarioS3_P1), 2)]
  )),
  legend = c("(1.5K, 2.5K, 2K)", "(1.5K, 3K, 2K)", "(1.5K, 3.5K, 2K)"),
  col = c("skyblue", "steelblue", "blue4"),
  pch = c(1, 2, 0),
  cex = c(1.2, 1.2, 1.2),
  bty = "n"
)


plot(
  as.numeric(ScenarioS6_P2[1, seq(1, ncol(ScenarioS6_P2), 2)]),
  ScenarioS6_P1[1, seq(1, ncol(ScenarioS6_P1), 2)],
  lwd = 1,
  type = "o",
  pch = 1,
  cex = 1.5,
  col = "tan2",
  xlab = "",
  ylab = "",
  ylim = c(
    min(ScenarioS6_P1[1, seq(1, ncol(ScenarioS6_P1), 2)], ScenarioS6_P1[2, seq(1, ncol(ScenarioS6_P1), 2)], ScenarioS6_P1[3, seq(1, ncol(ScenarioS6_P1), 2)]),
    (200 + max(
      ScenarioS6_P1[1, seq(1, ncol(ScenarioS6_P1), 2)], ScenarioS6_P1[2, seq(1, ncol(ScenarioS6_P1), 2)], ScenarioS6_P1[3, seq(1, ncol(ScenarioS6_P1), 2)]
    ))
  )
)
points(
  as.numeric(ScenarioS6_P2[2, seq(1, ncol(ScenarioS6_P2), 2)]),
  ScenarioS6_P1[2, seq(1, ncol(ScenarioS6_P1), 2)],
  lwd = 1,
  col = "coral2",
  pch = 2,
  type = "o",
  cex = 1.5
)
points(
  as.numeric(ScenarioS6_P2[3, seq(1, ncol(ScenarioS6_P2), 2)]),
  ScenarioS6_P1[3, seq(1, ncol(ScenarioS6_P1), 2)],
  lwd = 1,
  col = "red3",
  pch = 0,
  type = "o",
  cex = 1.5
)
legend(
  800,
  (100 + max(
    ScenarioS6_P1[1, seq(1, ncol(ScenarioS6_P1), 2)], ScenarioS6_P1[2, seq(1, ncol(ScenarioS6_P1), 2)], ScenarioS6_P1[3, seq(1, ncol(ScenarioS6_P1), 2)]
  )),
  legend = c("(3K, 1.5K, 2K)", "(3K, 1.5K, 2.5K)", "(3K, 1.5K, 3K)"),
  col = c("tan2", "coral2", "red3"),
  pch = c(1, 2, 0),
  cex = c(1.2, 1.2, 1.2),
  bty = "n"
)
