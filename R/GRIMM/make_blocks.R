source('./R/GRIMM/functions/make_grimm_command.R')


# Make anchors
system(make_grimm_command(), wait = T)

# Make blocks
system(make_grimm_command(A=F), wait=T)

rm(list=ls())
