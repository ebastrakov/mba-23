# A script to build a pH-Ctot diagram for uranyl complexes paper
# by Migdisov et al. (2023)
# As of 16 Sept 2023 uses a development version of CHNOSZ from GitHub

library(CHNOSZ)
OBIGT(no.organics = TRUE)
iNew <- add.OBIGT("data/mhkf-fitted-results-UO3-new.csv")

exp <- read.csv("data/kalintsev-2021.csv")
exp <- exp[exp$T == 200,]
# Note median(exp[exp$NaCl == 1,]$log_mU) is -5.16
# Note median(exp[exp$NaCl == 1,]$log_mU) is -5.44

chem_subsystem <- c("H","O","Cl","S","C","Na")
icr <- retrieve("U", chem_subsystem, state = "cr")
iaq <- retrieve("U", chem_subsystem, state = "aq")

basis(c("UO2+2", "H2O", "H+", "O2", "Cl-", "SO4-2", "CO2", "Na+"))

m_tot <- 1
u_tot <- -4.5

T <- 200
P <- "Psat"
NaCl <- NaCl(m_tot = m_tot, T = T, P = P) # , pH = pH)
IS = NaCl$IS

basis("Na+", log10(NaCl$m_Na))
basis("Cl-", log10(NaCl$m_Cl))
basis("SO4-2", -12)

species(icr)
species(iaq, u_tot, add = TRUE)

# Define basis species for mosaic diagram
# The first one in each group must be one of the species defined by basis()
bases <- c("CO2", "HCO3-", "CO3-2")
#bases <- c("SO4-2", "HSO4-", "HS-", "H2S")
## Use this to change both C and S species
#bases <- list(c("CO2", "HCO3-"), c("SO4-2", "HSO4-", "HS-", "H2S"))

# pH-Carbonate diagram for CO2 (aq) and HCO3- dominance fields
args <- list(bases, c(c(2,10), 600), c(c(-4,0.5), 600), T, P, IS)
names(args) <- c("bases", "pH", "CO2", "T", "P", "IS")
# mosaic() can now change basis species that are axis variables on a diagram 20230809
a <- do.call(mosaic, args)
diagram(a$A.species, cex = par("cex"), cex.names = 1, cex.axis = 1.5) #, main = "Mosaic diagram for CO2 (aq) and HCO3-")
diagram(a$A.bases, add=TRUE, fill = NA, lty=2, lwd=1, col="red", col.names="red", cex.names=1)
# points(exp$pH, log10(exp$mCO3), cex = 1.5, pch = 19, col = "#3C8DBC")
# points(c(4.77), log10(c(0.3)), cex = 1.5, pch = 19, col = "red")
# abline(v=-subcrt(c("CO2","H2O","HCO3-","H+"),c(-1,-1,1,1),T=T)$out$logK,lty=2, col = "blue")
# abline(h=log10(0.3), col = "red")
dTP <- describe.property(c("T", "P"), c(T, P))
dIS <- describe.property("IS", IS, digits = 1)
dNaCl<- substitute(y == x~mol~kg^-1, list(y="NaCl", x=m_tot))
dM <- substitute(y == x~mol~kg^-1, list(y="U", x=formatC(10^u_tot, format = "e", digits = 2)))
dL <- c(dTP, dIS, dNaCl, dM)
legend("bottomright", as.expression(dL), bty="y")


