source('./R/GRIMM/functions/make_grimm_command.R')


# Make anchors
# system(make_grimm_command(), wait = T)

# Make blocks
# system(make_grimm_command(A=F, input='/R/GRIMM/output_data/anchors/unique_coords.txt'), wait=T)
system(make_grimm_command(A=F, input='/R/Clean/output_data/GRIMM.txt'), wait=T)

# rm(list=ls())
