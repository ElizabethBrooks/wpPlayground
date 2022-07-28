2+2
200 * 5 / 2

southbendpop <- 102037
mishawakapop <- 49245

#How many people live in South Bend and Mishawaka?
southbendpop + mishawakapop

# Activity 1: Create a new object called "elkhartpop" for Elkhart's population total of 52257. Next, calculate how many people live in South Bend, Mishawaka, and Elkhart.

#Next, we'll move on to making an object called vectors. A vector strings together objects that are the same type (for example, numbers)

socialmediafollowers <- c(550, 700, 1000)
socialmediafollowers
mean(socialmediafollowers)
sum(socialmediafollowers)

names(socialmediafollowers) <- c("TikTok", "Instagram", "Twitter")
socialmediafollowers

socialmediafollowers[2]

### end of part 1.

#setwd("/Users/bamflappy/Repos/rPlayground/data/")

ships <- read.csv("PiratesShip.csv")

head(ships)

View(ships)

shipCapacity <- ships$CrewCapacity
shipCapacity
mean(shipCapacity)

ships$CrewCapacity
mean(ships$CrewCapacity)

summary(ships)

#install.packages("ggplot2")
library(ggplot2)

ggplot(data = ships, aes(x = CrewCapacity)) +
  geom_bar()

ggplot(data = ships, aes(x = Sails, y = CrewCapacity)) +
  geom_point()

## Activity 2:  Try to change variables around and produce at least 1 new “geom_point” visualization. Do you notice any trends? Or are there no definite trends?
#ex:
ggplot(data = ships, aes(x = MaidenYear, y = CrewCapacity)) +
  geom_point()

ggplot(data = ships, aes(x = Sails, y = CrewCapacity, color = MaidenYear)) +
  geom_point()

## Wes Anderson Section

#install.packages("wesanderson")
library(wesanderson)

names(wes_palettes)

zis_colors <- wes_palette("Zissou1", type = "discrete")
zis_colors

ggplot(data = ships, aes(x = Sails, y = CrewCapacity, color = MaidenYear)) +
  geom_point() +
  scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous"))
