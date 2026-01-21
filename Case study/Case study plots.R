# Figure 19----

# Figure 19 displays the density function and the histogram of the response variable.

# 1. Load necessary libraries

suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(tidyverse)))

# 2. Import database

Data <- read_xlsx(path = "~/Trabajo de grado/Database - Case study 1.xlsx")

# 2.1 Select the response variable

Data <- Data %>%
  select(CODE, # Country code
         GR6096, # Average growth rate of Gross Domestic Product (GDP) per capita between 1960 and 1996
         GDPCH60L, # GDP per capita in 1960 (logaritmic scale)
         LIFE060, # Life expectancy in 1960
         P60, # Primary school enrollment rate in 1960
         SAFRICA, # Dummy variable for Sub-Sahara Africa
         LAAM, # Dummy variable for Latin America
         EUROPE, # Dummy variable for Europe
         EAST, # Dummy variable for Asia
         SCOUT, # Dummy variable for Outward orientation
         DPOP6090, # Growth rate of population between 1960 and 1990
         H60, # Higher education enrollment rate in 1960
         YRSOPEN, # Number of years an economy has been open between 1950 and 1994
         REVCOUP, # Number of military coups and revolutions
         WARTORN, # Dummy variable for countries that have been involved in war any time between 1960 and 1990
         PRIGHTS, # Index of political rights
         CIV72, # Index of civil liberties
         ABSLATIT, # Absolute latitude
         AVELF, # Index of ethnolinguistic fractionalization
         PRIEXP70, # Fraction of primary exports in total exports in 1970
         RERD, # Real exchange rate distortions
         BRIT, # Dummy variable for former British colonies
         SPAIN, # Dummy variable for former Spanish colonies
         BUDDHA, # Percentage of the population that follows the Buddhist religion
         CATH00, # Percentage of the population that follows the Catholic religion
         CONFUC, # Percentage of the population that follows Confucian religion
         HINDU00, # Percentage of the population that follows the Hindu religion
         MUSLIM00, # Percentage of the population that follows the Muslim religion
         PROT00, # Percentage of the population that follows the Protestant religion
         MINING, # Percentage of the GDP in mining
         ECORG, # Index of degree in which economies favor capitalist forms of production
         OTHFRAC, # Percentage of the population speaking foreign language
         ENGFRAC # Percentage of the population speaking english
  )

# 2.2 Remove observations with some missing values

Data <- Data[complete.cases(Data), ]

y <- Data$GR6096 # Set the response variable

# 3. Plot histogram

hist(x = y, freq = FALSE, xlim = c(-0.06, 0.08), 
     ylab = "Densidad", main = "",
     col = alpha("grey", 0.3), cex.label = 1.5, cex.axis  = 1.5)
# Density plot using a smooth kernel density plot
lines(density(x = y), lwd = 2, col = "darkorange1")

# Figure 20----

# Figure 20 displays the density function and the histogram of the response variable.

# 1. Load necessary libraries

suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(tidyverse)))

# 2. Import database

Data <- read_csv("~/Trabajo de grado/Database - Case study 2.txt")

# 2.1 Select the response variable

y <- Data$Wage # Set the response variable

# 3. Plot histogram

hist(x = y, freq = FALSE,
     ylab = "Densidad", main = "", xlim = c(10, 18), ylim = c(0, 1),
     col = alpha("grey", 0.3), cex.label = 1.5, cex.axis  = 1.5)
# Density plot using a smooth kernel density plot
lines(density(x = y), lwd = 2, col = "indianred1")
