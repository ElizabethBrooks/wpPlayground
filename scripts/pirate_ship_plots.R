#### Colorful R Plots with Wes Anderson Palettes & ggplot2 - Pirate Ship Features

# import the pirate ship data from the PiratesShip.csv file
ships <- read.csv("PiratesShip.csv")

# look at the column names of the ship data
colnames(ships)

# look at the first few lines of ship data
head(ships)

# view the ship data set in a new tab
View(ships)

### One Dimension of Data with Vectors - CrewCapacity or Sails
# create a plot for each dimension of CrewCapacity or Sails

# load the ggplot2 library
library(ggplot2)

# look at the CrewCapacity vector (dimension) of the ships data
ships$CrewCapacity

# plot only the CrewCapacity dimension of the ship data
ggplot(data = ships, aes(x = CrewCapacity)) +
  geom_bar()

# look at the Sails vector (dimension) of the ships data
ships$Sails

# plot only the Sails dimension of the ship data
ggplot(data = ships, aes(x = Sails)) +
  geom_bar()

# look at the MaidenYear vector (dimension) of the ships data
ships$MaidenYear

# plot only the MaidenYear dimension of the ship data
ggplot(data = ships, aes(x = MaidenYear)) +
  geom_bar()

### Two Dimensions of Data with Dataframes - CrewCapacity & Sails
# combine the two dimensions into one plot with geom_point

# first, create a scatter plot
ggplot(data = ships, aes(x = Sails, y = CrewCapacity)) +
  geom_point()

### Three Dimensions of Data with Colors - CrewCapacity & Sails & MaidenYear
# add a third dimension to your plot with color

# color the scatter plot by MaidenYear
ggplot(data = ships, aes(x = Sails, y = CrewCapacity, color = MaidenYear)) +
  geom_point()

### Fun Colors - Wes Anderson Palette
# create a fun colorful plot using the wesanderson color palette

# load the wesanderson package
library(wesanderson)

# checkout the color palette options
names(wes_palettes)

# look at the info for the wes_palette function
?wes_palette

# look at the info for the scale_color_gradientn function
?scale_color_gradientn

# color the scatter plot using wes_palettes and scale_color_gradientn
ggplot(data = ships, aes(x = Sails, y = CrewCapacity, color = MaidenYear)) +
  geom_point() +
  scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous"))

### Bonus Exercises - Adjusting Plot Appearance

# add axis and plot titles
ggplot(data = ships, aes(x = Sails, y = CrewCapacity, color = MaidenYear)) +
  geom_point() +
  scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous")) +
  labs(title = "Pirate Ship Crew Capacity by Sails and Maiden Year", 
       x ="Number of Sails", 
       y = "Crew Capacity")

# retrieve the vector of colors associated with Zissou1
(zis_colors <- wes_palette("Zissou1", type = "discrete"))

# adjust the axis and plot title colors
ggplot(data = ships, aes(x = Sails, y = CrewCapacity, color = MaidenYear)) +
  geom_point() +
  scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous")) +
  labs(title = "Pirate Ship Crew Capacity by Sails and Maiden Year", 
       x ="Number of Sails", 
       y = "Crew Capacity") +
  theme(
    plot.title = element_text(color = zis_colors[1], size = 14, face = "bold.italic"),
    axis.title.x = element_text(color = zis_colors[4], size = 14, face = "bold"),
    axis.title.y = element_text(color = zis_colors[5], size = 14, face = "bold")
  )

# center the plot title
ggplot(data = ships, aes(x = Sails, y = CrewCapacity, color = MaidenYear)) +
  geom_point() +
  scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous")) +
  labs(title = "Pirate Ship Crew Capacity by Sails and Maiden Year", 
       x ="Number of Sails", 
       y = "Crew Capacity") +
  theme(
    plot.title = element_text(color = zis_colors[1], size = 14, face = "bold.italic", hjust = 0.5),
    axis.title.x = element_text(color = zis_colors[4], size = 14, face = "bold"),
    axis.title.y = element_text(color = zis_colors[5], size = 14, face = "bold")
  )

### Saving Plots - ggsave

# check out the info for the ggsave function
?ggsave

# check out the info for the last_plot function
?last_plot

# save the last plot using the ggsave function
ggsave("plots/bonus/ship_plot_crew_sails_year_centerTitle.png", plot = last_plot())
