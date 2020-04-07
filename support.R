# Construct fictitious initial (S,E,I,R) = (1, 0, I0, 0) state at t=Start.
# We find the I0 that generates the user-specified host fraction (E+I) at hosts_guess_t.
get_I0 <- function(hosts_guess, hosts_guess_t, start_time, dt, seir_model, seir_model_params){
#print(">>>hosts_guess, hosts_guess_t>>>")
#print(hosts_guess)
#print(hosts_guess_t)
#print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")


  projected_hosts_error <- function(I_start){
    trajectory <- lsoda(c(S=1, E=0, I=I_start, R=0), times = seq(start_time, hosts_guess_t, by=dt), seir_model, seir_model_params, rtol=1e-7, atol=1e-7)
    last_index = dim(trajectory)[1]
    (trajectory[last_index, 3] + trajectory[last_index, 4]) - hosts_guess
  }
  uniroot(projected_hosts_error, c(0,hosts_guess), tol = .Machine$double.eps^0.5)$root
}

# Interface to the built-in filter() function: truncates the weighted
# sum (rather than generating NAs) whenever the initial segment of 
# input_series history is shorter than the vector of filter weights.
zp_filter <- function(input_series, linear_filter){
  zp_input_series <- c(rep(0, length(linear_filter)-1), input_series)
  zp_output_series <- 
    filter(zp_input_series, linear_filter, method='convolution', sides = 1)
  zp_output_series[-(1:(length(linear_filter)-1))]
}
