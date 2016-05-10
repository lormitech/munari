
#
# Copyright (c) 2009-2010, 2014, 2016
# Author: Lorenzo Milazzo
#
# All rights reserved.
#



#
# last modified: 10.05.16
#

#
# vers. 0.4.0 (128)
#



#
#              -- Multifractal System (MFS) Analyser --
#


#--

# Munari root directory
munari_root = "<path>/munari/"

# I/O processing
# reading input file
data = scan(paste(munari_root, "working_dir/singularity_spectrum.dat", sep=""), list(x=0, y=0, z=0, w=0))


#--

# -- box sizes --

#
# box size: l=2, l=4, l=8, l=16, l=32, l=64
#

# no. of box sizes
no_box_sizes = 6

l_values = c(data$y[1:no_box_sizes])
l_values
l_values_log = log(l_values)
l_values_log


# -- order moments --

#
# order moment (q), -5 <= q <= 5
#

# no. of order moments
no_q = 21

q_values_t = c(data$x)
q_values = rep(0, times=no_q)
for (i in 1:no_q) {
  q_values[i] = q_values_t[1 + (i-1)*no_box_sizes];
}
q_values



#--

#
# 'alpha(q)' (singularity strength):
#             alpha(q) = lim, l->0 a_q_numerator/a_q_denominator
# where
#      a_q_numerator   = Sum{mu_i(q,l) * ln(P_i(l))}
#      a_q_denominator = ln(l)
#

a_q_num_t = c(data$z)
a_q_num_t

# setting infinite values to NaN values 
for (i in 1:length(a_q_num_t)) {
  if(is.infinite(a_q_num_t[i])==TRUE){
    a_q_num_t[i] = NaN
  }
}

# slopes of the curves (a_q_numerator vs. a_q_denominator)
# for different values of q
a_q = rep(0, times=no_q)

lbound = 1
ubound = no_box_sizes
for (i in 1:no_q) {
  a_q_num = rep(0, times=no_box_sizes)
  j = 1
  for (w in lbound:ubound) {
    a_q_num[j] = a_q_num_t[w];
    j = j + 1
  }
  print(a_q_num)
  # if there is degradation of the scaling,
  # applying a correction to obtain a straight line
  #a_q_num[1]<- NaN
  a_q_num[2]<- NaN
  print(a_q_num)

  # evaluating the regression coefficients
  lm_data_1 = lm(a_q_num ~ l_values_log)
  a_q[i] = lm_data_1$coefficients[2]
  lbound = lbound + no_box_sizes
  ubound = ubound + no_box_sizes  
}
a_q



#--

#
# 'f(q)' (Hausdorff dimension of normalized measures):
#             f(q) = lim, l->0 f_q_numerator/f_q_denominator
# where
#      f_q_numerator   = Sum{mu_i(q,l) * ln(mu_i(q,l))}
#      f_q_denominator = ln(l)
#

f_q_num_t = c(data$w)
f_q_num_t
# (check the data to see if there is degradation of the scaling,
# expecially for negative values of q; if there is degradation of
# the scaling, the scaling range is smaller)
#

# slopes of the curves (f_q_numerator vs. f_q_denominator)
# for different values of q
f_q = rep(0, times=no_q)

lbound = 1
ubound = no_box_sizes
for (i in 1:no_q) {
  f_q_num = rep(0, times=no_box_sizes)
  j = 1
  for (w in lbound:ubound) {
    f_q_num[j] = f_q_num_t[w];
    j = j + 1
  }
  print(f_q_num)
  # if there is degradation of the scaling,
  # applying a correction to obtain a straight line
  #f_q_num[1]<- NaN
  f_q_num[2]<- NaN
  print(f_q_num)

  # evaluating the regression coefficients
  lm_data_2 = lm(f_q_num ~ l_values_log)
  f_q[i] = lm_data_2$coefficients[2]
  lbound = lbound + no_box_sizes
  ubound = ubound + no_box_sizes  
}
f_q



#--

# I/O processing
# writing output files

#
# plotting (a_q_numerator vs. a_q_denominator) and
# regression line for q=-5.0
#

lbound = 1
ubound = no_box_sizes
a_q_num = rep(0, times=no_box_sizes)
j = 1
for (i in lbound:ubound) {
  a_q_num[j] = a_q_num_t[i];
  j = j + 1
}
a_q_num
lm_data_1 = lm(a_q_num ~ l_values_log)
lm_data_1


xrange = c(0,7)
yrange = c(-25,0)
postscript(paste(munari_root, "working_dir/a_q_q01.ps", sep=""))
plot(l_values_log, a_q_num, type="o", xlim=xrange, ylim=yrange, xlab="log(l)", ylab="a_q_num")
abline(lm_data_1, lty=2)
dev.off()



#
# plotting (f_q_numerator vs. f_q_denominator) and
# regression line for q=-5.0
#

lbound = 1
ubound = no_box_sizes
f_q_num = rep(0, times=no_box_sizes)
j = 1
for (i in lbound:ubound) {
  f_q_num[j] = f_q_num_t[i];
  j = j + 1
}
f_q_num

lm_data_2 = lm(f_q_num ~ l_values_log)
lm_data_2

xrange = c(0,7)
yrange = c(-15,5)
postscript(paste(munari_root, "working_dir/f_q_q01.ps", sep=""))
plot(l_values_log, f_q_num, type="o", xlim=xrange, ylim=yrange, xlab="log(l)", ylab="f_q_num")
abline(lm_data_2, lty=2)
dev.off()



#
# plotting slopes of the curves:
#   a) (a_q_numerator vs. a_q_denominator)
#   b) (f_q_numerator vs. f_q_denominator)
# for different values of q
#

postscript(paste(munari_root, "working_dir/af_q.ps", sep=""))
xrange = c(-6.5,5.5)
yrange = c(-1,8.0)
# line type: 0=blank, 1=solid, 2=dashed, 3=dotted,
#            4=dotdash, 5=longdash, 6=twodash
l_types = c(1, 2)
plot(xrange, yrange, type="n", xlab="q", ylab="alpha, f")
lines(q_values, a_q, lty=l_types[1])
lines(q_values, f_q, lty=l_types[2])
legend(x="topright", y=NULL, c("alpha", "f"), cex=1, lty=l_types)
dev.off()



#
# plotting singularity spectrum 'f(alpha)'
#

xrange = c(2.0,7.0)
yrange = c(-1,4.0)
postscript(paste(munari_root, "working_dir/singularity_spectrum.ps", sep=""))
plot(a_q, f_q, type="o", xlim=xrange, ylim=yrange, xlab="alpha", ylab="f")
title("singularity spectrum - f(alpha)")
dev.off()
