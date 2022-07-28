#### Fun R Plots with ggplot2 - Pirate Features

# import the pirate data from the PiratesPirate.csv file
pirates <- read.csv("PiratesPirate.csv")

# look at the column names of the pirates data
colnames(pirates)

# look at the first few lines of pirates data
head(pirates)

# view the pirates data set in a new tab
View(pirates)

### One Dimension of Data with Vectors - Teeth or Limbs
# create a plot for each dimension of Teeth or Limbs

# load the ggplot2 library
library(ggplot2)

# look at the Teeth vector (dimension) of pirate data
pirates$Teeth

# plot only the Teeth dimension of the pirates data
ggplot(data = pirates, aes(x = Teeth)) +
  geom_bar()

# look at the Limbs vector (dimension) of pirate data
pirates$Limbs

# plot only the Limbs dimension of the pirates data
ggplot(data = pirates, aes(x = Limbs)) +
  geom_bar()

### Two Dimensions of Data with Dataframes - Teeth & Limbs
# explore the relationship between the numbers of Teeth and Limbs
# combine the two dimensions into one plot with geom_point

# create a scatter plot comparing the two dimensions of Limbs vs Teeth
ggplot(data = pirates, aes(x = Limbs, y = Teeth)) +
  geom_point()

# modify the scatter plot to include the count of pirates at each point
ggplot(data = pirates, aes(x = Limbs, y = Teeth)) +
  geom_count()

# check the data for the outlier pirate with 2.5 Limbs
pirates[pirates$Limbs == "2.5", ]

# change to a boxplot to better describe the data at each point
ggplot(data = pirates, aes(x = Limbs, y = Teeth)) +
  geom_boxplot()

# check the info for the factor function
?factor

# use the factor function to make the different numbers of Limbs into categories
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth)) +
  geom_boxplot()

### Three Dimensions of Data with Colors - Teeth & Limbs & Origin
# explore the relationship between the numbers of Teeth and Limbs by Origin
# add another dimension to your plots with color

# first, change the color of the boxes
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth)) +
  geom_boxplot(color="red", fill="orange", alpha=0.2)

# next, color the boxes by Origin 
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth, fill = Origin)) +
  geom_boxplot(alpha=0.3) +
  theme(legend.position="none")

# create a set of boxplots with one for each Origin of pirates
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth)) +
  geom_boxplot(color="darkgreen", fill="orange", alpha=0.2) +
  facet_wrap(~ Origin)

### Bonus Exercises - Adjusting Plot Appearance

# change the appearance of the individual plot titles
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth)) +
  geom_boxplot(color = "darkgreen", fill = "orange", alpha = 0.2) +
  facet_wrap(~ Origin) +
  theme(strip.background = element_rect(colour = "black", fill = "white", 
                                        size = 1.5, linetype = "solid"))

# change the axis titles and add a plot title
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth)) +
  geom_boxplot(color = "darkgreen", fill = "orange", alpha = 0.2) +
  facet_wrap(~ Origin) +
  theme(strip.background = element_rect(colour = "black", fill = "white", 
                                        size = 1.5, linetype = "solid")) +
  labs(title = "Comparison of Pirate Teeth and Limb Numbers by Origin", 
       x ="Number of Limbs", 
       y = "Number of Teeth")

# adjust the colors of the axis and plot titles
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth)) +
  geom_boxplot(color = "darkgreen", fill = "orange", alpha = 0.2) +
  facet_wrap(~ Origin) +
  theme(strip.background = element_rect(colour = "black", fill = "white", 
                                        size = 1.5, linetype = "solid")) +
  labs(title = "Comparison of Pirate Teeth and Limb Numbers by Origin", 
       x ="Number of Limbs", 
       y = "Number of Teeth") +
  theme(
    plot.title = element_text(color = "red", size = 14, face = "bold.italic"),
    axis.title.x = element_text(color = "blue", size = 14, face = "bold"),
    axis.title.y = element_text(color = "purple", size = 14, face = "bold")
  )

# center the plot title
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth)) +
  geom_boxplot(color = "darkgreen", fill = "orange", alpha = 0.2) +
  facet_wrap(~ Origin) +
  theme(strip.background = element_rect(colour = "black", fill = "white", 
                                        size = 1.5, linetype = "solid")) +
  labs(title = "Comparison of Pirate Teeth and Limb Numbers by Origin", 
       x ="Number of Limbs", 
       y = "Number of Teeth") +
  theme(
    plot.title = element_text(color = "red", size = 14, face = "bold.italic", hjust = 0.5),
    axis.title.x = element_text(colo  = "blue", size = 14, face = "bold"),
    axis.title.y = element_text(color = "purple", size = 14, face = "bold")
  )

### Saving Plots - ggsave

# check out the info for the ggsave function
?ggsave

# check out the info for the last_plot function
?last_plot

# save the last plot using the ggsave function
ggsave("pirates_plot_teeth_limbs_origin.png", plot = last_plot())
