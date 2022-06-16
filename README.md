# How-Priority-effects-alter-disease-risk
R code for the publication "Exploring how a generalist pathogen and within-host priority effects alter the risk of being infected by a specialist pathogen"
R code for all the figures in "Exploring how a generalist pathogen and within-host priority effects alter the risk of being infected by a specialist pathogen"
by Jiao & Cortez 2022

The names of all files are arranged by:

FigureX: indicates the figure label in main text or Appendix (e.g., Figure2 corresponds to the Figure 2 in the main text, FigureS1 is the Figure S1 in Appendix)

Scenarios_SX_SY: indicates the simulated relationship between generalist and specialist infectious propagule densities at equilibrium (P1-P2 relationship) when within-host priority effects exhibit scenarios X to Y in Shedding according to Table 2
                 e.g., Scenario_S1_S6 represents simulated P1-P2 relationship at equlibrium at scenario 1-6 in Shedding according to Table 2

Scenarios_MX_MY: indicates the simulated relationship between generalist and specialist infectious propagule densities at equilibrium (P1-P2 relationship) when within-host priority effects exhibit scenarios X to Y in Mortality according to Table 2
                  e.g., Scenario_M1_M6 represents simulated P1-P2 relationship at equlibrium at scenario 1-6 in Mortality according to Table 2

Scenarios_SX_MX: indicates the simulated relationship between generalist and specialist infectious propagule densities at equilibrium (P1-P2 relationship) when within-host priority effects exhibit scenarios X in Shedding and Y in Mortality according to Table 2
			e.g., Scenario_S4-6_M1-3 represents simulated P1-P2 relationship at equlibrium at combined scenarios with scenario 4-6 in Shedding and 1-3 in Mortality according to Table 2


Note: The scenario 1-3 in Shedding (i.e.,S1-3) or in Mortality (i.e., M1-3) are cases with first arriver advantage for within-host priority effects;
      The scenario 4-6 in Shedding (i.e.,S4-6) or in Mortality (i.e., M4-6) are cases with second arriver advantage for within-host priority effects.


Sym: indicate that both pathogens exhibit same type of within-host priority effects in Shedding and Mortality at focal host species (or both pathogens show symmetric priority effects in Shedding and Mortality at focal host)
     e.g., Sym after "S4-6_M1-3" means that both specialist and generalist pathogen exhibit second arriver advantage (scenario 4-6) in Shedding and first arriver advantage (scenario 1-3) in Mortality at focal host "F"

Asym: indicate that two pathogens exhibit different types of within-host priority effects in Shedding and Mortality at focal host species (or, two pathogens show asymmetric priority effects in Shedding and Mortality at focal host)
     e.g., Asym after "S4-6_M1-3" means that specialist pathogen exhibit second arriver advantage (scenario 4-6) in Shedding but first arriver advantage (scenario 1-3) in Mortality at focal host "F", while generalist pathogen
     exhibits first arriver advantage (scenario 1-3) in Shedding but second arriver advantage (scenario 4-6) in Mortality at focal host "F".

Note: in all the Sym or Asym cases, generalist pathogen is assumed not to exhibit any within-host priority effects in alternative host "A", i.e., all shedding or mortality at all infected classes are equal for generalist pathogen at alternative host species "A". 



 

