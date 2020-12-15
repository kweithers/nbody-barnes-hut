library(tidyverse)
library(gganimate)
setwd("~/Documents/advancedProgramming/homework-3-kweithers")

df <- read_delim('particles.dat',delim = ' ',skip=1,col_names = c('t','x','y'))
df['Particle'] = as.factor(rep(1:16,1000))
p<-ggplot(df,aes(x=x,y=y,color=Particle)) +
  geom_point() +
  xlim(-.2,1.2) +
  ylim(-.2,1.2) +
  geom_vline(xintercept = c(0,1)) +
  geom_hline(yintercept = c(0,1))
plot(p)

anim <- p + 
  transition_states(t,transition_length = 3, state_length = 3) + 
  ggtitle('Barnes-Hut Approximation of n-Body Problem', subtitle = 'Timestep {frame} of {nframes}')
animate(anim,nframes=200)
# anim_save('16particles.gif')