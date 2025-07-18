num_analisys = function( draws, burn = 0, lags = 1, names, digits, hdp = FALSE ){
  if( !require(coda) ) install.packages("coda")
  
  ############### Numeric Analysis
  Draws = draws
  N = ncol( Draws ) 
  if( burn != 0 ) Draws = Draws[, -c( 1:burn )]
  jumps = seq(1, N - burn, by = lags)
  Draws = Draws[, jumps ]
  
  chain = coda::as.mcmc( t( Draws ) )
  CD = coda::geweke.diag( chain )
  N_eff = coda::effectiveSize( chain )
  IF = ncol( Draws ) / N_eff
  mc_error = apply( chain, 
                    MARGIN = 2, 
                    FUN = sd) / sqrt( N_eff )
  theta_hat = apply( chain, MARGIN = 2, FUN = mean )
  theta_sd  = apply( chain, MARGIN = 2, FUN = sd )
  theta_median = apply( chain, MARGIN = 2, FUN = median )
  if( hdp ){
    i = coda::HPDinterval( chain )
    theta_min = i[, 1]
    theta_max = i[, 2]
  }else{
    theta_min = apply( chain, MARGIN = 2, FUN = quantile, probs = c(0.025) )
    theta_max = apply( chain, MARGIN = 2, FUN = quantile, probs = c(0.975) )
  }
  
  data = matrix(
    c(theta_hat,
      theta_sd,
      theta_min,
      theta_median,
      theta_max,
      CD$z,
      IF,
      N_eff,
      mc_error), nrow = nrow( Draws ), byrow = FALSE
  )
  row.names( data ) = names
  if( hdp ){
    colnames( data ) = c( 'media', 'sd', 'HPD.min', '50%','HPD.max', 'CD', 'IF', 'n_eff', 'MC erro')
  }else{
    colnames( data ) = c( 'media', 'sd', '2.5%', '50%','97.5%', 'CD', 'IF', 'n_eff', 'MC erro')
  }
  
  return( round( data, digits ) )
}