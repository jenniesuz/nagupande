#*******************host kill and tsetse count data****************************
dat <- read.csv("nagupande_data_R.csv") # read in dataset
names(dat)
# columns
# nagupande: mean monthly G. morsitans counts at Nagupande
# lusulu: mean monthly G. morsitans counts at Lusulu control site
# animal columns: numbers killed each month

dat$time <- seq(15,69*30,30) # create column of estimated times (days) taking the middle of each month
