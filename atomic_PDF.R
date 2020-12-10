
library(tidyverse)

# sample gold NP coordinates are given in meters, convert to Angstrom when reading the file

read_delim("C:/Users/new/Downloads/3589atoms_4.08A.txt", " ", col_names = F) %>% rename(x=1,y=2,z=3) %>% mutate_all(as.double) %>% mutate_all(~.*1e10) %>% as.matrix()-> atomic_arrangement
N <- nrow(atomic_arrangement)


# calculate interatomic separations rij

x <- atomic_arrangement[,1]
y <- atomic_arrangement[,2]
z <- atomic_arrangement[,3]

comb_id <- combn(N,2)

as.vector(comb_id)[seq(1,length(comb_id),2)]->iii
as.vector(comb_id)[seq(2,length(comb_id),2)]->jjj

xxi <- x[iii]
xxj <- x[jjj]
yyi <- y[iii]
yyj <- y[jjj]
zzi <- z[iii]
zzj <- z[jjj]
rij_all <- sqrt((xxi-xxj)^2+(yyi-yyj)^2+(zzi-zzj)^2)

# calculate F(Q)

Qmin=1 
Qmax=4*pi
sigma=0.0 # correlated broadening factor for the atom pair (debyer-waller)

hist(rij_all,breaks = seq(0,50,1e-3),plot = F)->rij_bin
rij_bin$breaks[1:length(rij_bin$breaks)-1]-> rij_bin_b
rij_bin$counts-> rij_bin_c

Q_range <- seq(Qmin,Qmax,1e-2)
F_Q <- 0

for (i in 1:length(Q_range)) {
  q_part <- exp(-1/2*sigma^2*Q_range[i]^2)*rij_bin_c*sin(Q_range[i]*rij_bin_b)/rij_bin_b
  F_Q[i] <- 1/N*(2*sum(q_part,na.rm = TRUE))
}

# fourier transformation of F(Q)

G_r <- 0
tst_r <- seq(0,50,0.0001)
dQ <- 1e-2
for (i in 1:length(tst_r)) {
  G_r[i] <- 2/pi*sum(F_Q*sin(Q_range*tst_r[i])*dQ)
}

# plot Gr
tibble(r=tst_r,G_r=G_r) %>% ggplot(aes(r,G_r))+geom_line()+scale_x_continuous(limits = c(0,50),n.breaks = 10)



































