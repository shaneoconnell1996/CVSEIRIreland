
# Updated to include arguments intervention_t and intervention_t2 to make function
# standalone (& not dependant on any globals etc.)
# 
seir_model = function (t, state_values, parameters)
{
  S = state_values [1]        # susceptibles
  E = state_values [2]        # exposed
  I = state_values [3]        # infectious
  R = state_values [4]        # recovered
  alpha = parameters[1]
  if(t < parameters[6]) {  # intervention_t
    beta = parameters[2]
  } else if (t < parameters[7]) {  # intervention_t2
    beta = parameters[3]
  } else {
    beta = parameters[5]
  }
  gamma = parameters[4]
  dS = (-beta * S * I)
  dE = (beta * S * I) - (alpha * E)
  dI = (alpha * E) - (gamma * I)
  dR = (gamma * I)
  results = c (dS, dE, dI, dR)
  list (results)
}
