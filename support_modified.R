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
    stats::filter(zp_input_series, linear_filter, method='convolution', sides = 1)
  zp_output_series[-(1:(length(linear_filter)-1))]
}

create_df <- function(start_t,infectious_period=2.9,latent_period=5.2,R0_pre,R0_post,R0_post2,CV,mu,hosts_count_guess,ICU_rate){
  Start = start_t
  hosts_guess_t = 1     # user will supply guess re number of hosts at this time
  
  intervention_t = 13 # 1st intervention
  intervention_t2 = 27 # 2nd intervention
  End = 60
  dt = 1/20
  today = as.numeric(difftime(Sys.Date(),as.Date("2020-3-1")))
  
  
  
  # Other model parameters (empirical)
  
  # User may also want to experiment with these values, but
  # app should make clear that the default settings are well grounded
  # in the scientific literature
  #infectious_period = 2.9 # from Goh epidemic simulator
  #latent_period = 5.2     # from epidemic simulator (will need better source)
  
  #Seir model parameters based on above assumptions
  alpha = 1 / latent_period
  gamma = 1 / infectious_period
  beta1 = gamma * R0_pre
  beta2 = gamma * R0_post
  beta3 = gamma * R0_post2
  parameter_list = c (alpha, beta1, beta2, gamma, beta3, intervention_t, intervention_t2)
  
  # Compute initial SEIR state
  N = 4921500    # population of Republic (Wikipedia)
  # Number of hosts (exposed + infectious) at hosts_guess_t / Population
  hosts_guess = hosts_count_guess / N
  
  # compute I0 such that SEIR state (S=1, E=0, I=I0, R=0) at t=Start
  # will evolve, under SEIR dynamics, to a state with E+I = hosts_guess
  # at t = hosts_guess_t
  I0 = get_I0(hosts_guess, hosts_guess_t, Start, dt, seir_model, parameter_list)
  initial_values = c(S = 1, E = 0, I = I0, R = 0)
  
  # Solve.
  timepoints = seq (Start, End, by=dt)
  output = lsoda (initial_values, timepoints, seir_model, parameter_list, rtol=1e-7, atol=1e-7)
  
  # Estimate ICU admissions based on user-inputted ICU_rate, k, mu and CV.
  # filter width: number of mus to scan back. User doesn't need to play with this.
  # Compute ICU admissions per timestep, and hence also cumulative admissions.
  # Admissions rate now will be proportional to a weighted average of E values
  # over the past k*mu days.
  k = 4
  filter_weights <- ICU_rate*N*dt*dgamma(seq(0, k*mu, dt), shape=1/CV^2, rate=1/(mu * CV^2))
  admissions <- zp_filter(alpha*dt*output[,'E'], filter_weights)
  cumulative_admissions <- cumsum(admissions)
  
  # HSE data from date of first intervention + 3.
  hse = read.csv('hse_data.csv')
  
  # Number datapoints in the csv file.
  icupoints = nrow(hse)
  
  # Get per day increases on ICU admissions.
  hse = cbind(day=((intervention_t+3):(intervention_t+3+icupoints-1)),hse)
  daily_icu = c()
  daily_icu[2:icupoints] = hse$icu[2:icupoints] - hse$icu[1:(icupoints-1)]
  daily_icu[1] = hse$icu[1]
  hse = cbind(hse,daily_icu)
  
  # Get per day deaths increases.
  daily_deaths = c()
  daily_deaths[2:icupoints] = hse$deaths[2:icupoints] - hse$deaths[1:(icupoints-1)]
  daily_deaths[1] = hse$deaths[1]
  hse = cbind(hse,daily_deaths)
  
  
  # Plots
  disp = c(3,4)
  dashed_size = 0.4
  df = data.frame(time=output[,1],output[,disp]*N, ICU=cumulative_admissions)
  df = melt(df,id.vars = 'time', variable.name = 'series')
  df2 = data.frame(time=output[,1], ICU= cumulative_admissions)
  df2 = melt(df2,id.vars = 'time', variable.name = 'series')
  return(list(df,df2,hse))
}