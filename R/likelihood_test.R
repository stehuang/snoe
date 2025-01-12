
# function: extract response variable from model formula
# input: model formula from GAM
# output: name of response variable
getResponseFromFormula = function(formula) {
  # compatible with GAM, but not lm
  if (attr(terms(as.formula(formula)), which = 'response'))
    return(all.vars(formula)[1])
  else
    return(NULL)
}


# function: check correlation between response vs smooth term; if near 1 for both directions, then relation is linear
# input: two models and their response variables
# output: TRUE/FALSE if components are both linear
check_linear_relation <- function(m1, m2, m1_response, m2_response){
  m1_predictor <- paste('s(', m2_response, ')', sep="")
  m2_predictor <- paste('s(', m1_response, ')', sep="")
  if(class(m1)[1] == 'list'){
    # extract fitted smooth terms
    # m1_smooth_pred_1 <- predict(m1[[1]], type='terms')[, m1_predictor]
    # m1_smooth_pred_2 <- predict(m1[[2]], type='terms')[, m1_predictor]
    # m2_smooth_pred_1 <- predict(m2[[1]], type='terms')[, m2_predictor]
    # m2_smooth_pred_2 <- predict(m2[[2]], type='terms')[, m2_predictor]
    m1_smooth_pred <- c(predict(m1[[1]], type='terms')[, m1_predictor], predict(m1[[2]], type='terms')[, m1_predictor])
    m2_smooth_pred <- c(predict(m2[[1]], type='terms')[, m2_predictor], predict(m2[[2]], type='terms')[, m2_predictor])
    # extract original values to predict
    m1_y <- c(m1[[1]]$y, m1[[2]]$y)
    m2_y <- c(m2[[1]]$y, m2[[2]]$y)
  }else{
    # when we only do one fold training of models
    m1_smooth_pred <- predict(m1, type='terms')
    m2_smooth_pred <- predict(m2, type='terms')
    if(ncol(m1_smooth_pred)>1){
      m1_smooth_pred <- predict(m1, type='terms')[, m1_predictor]
    }
    if(ncol(m2_smooth_pred)>1){
      m2_smooth_pred <- predict(m2, type='terms')[, m2_predictor]
    }

    # m1_smooth_pred <- predict(m1, type='terms')[, m1_predictor]
    # m2_smooth_pred <- predict(m2, type='terms')[, m2_predictor]
    m1_y <- m1$y
    m2_y <- m2$y
  }
  # arbitrary threshold
  if(abs(cor(m1_smooth_pred, m1_y)) > 0.975 & abs(cor(m2_smooth_pred, m2_y)) > 0.975){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

# function: get model df
# input: model of gam/lm or vector of values
# output: degrees of freedom
# https://stats.stackexchange.com/questions/346379/calculating-total-estimated-degrees-of-freedom-for-a-gam
get_model_df <- function(mod){
  if(class(mod)[1]=="gam"){
    return(sum(mod$edf1))
  }else if(class(mod)[1]=="lm"){
    return(length(coef(mod)))
  }else{
    return(1) # node has no parents
  }
}


# function: estimate the standard error of samples from residuals
# input: model object
# output: estimated error
est_normsd <- function(mod){
  est_sd <- sqrt(sum((mod$y-mod$fitted.values)^2)/mod$df.residual)
  return(est_sd)
}


# function: estimate the log-lik from models constructed under one edge direction
# input: model for child node, model for parent node, test data, parent node, child node
# output: list containing (1) individual log-likelihood of p(ch | pa) and (2) p(pa |pa(pa))
gam_ll <- function(m_main, m_prior, test_data, pa_node, ch_node){
  # get predicted values for conditional prob
  pred_val_main <- predict(m_main, test_data)
  if(class(m_main)[1] == "lm"){
    residual_se_main <- summary(m_main)$sigma
    ll_pa <- dnorm(test_data[,ch_node], mean = pred_val_main, sd = residual_se_main, log = TRUE)
  }else{ # this is for GAM
    # residual_se_main <- sqrt(m_main$scale)
    residual_se_main <- est_normsd(m_main)
    ll_pa <- dnorm(test_data[,ch_node], mean = pred_val_main, sd = residual_se_main, log = TRUE)
  }

  # get predicted values for prior prob
  if(class(m_prior)[1]=="numeric"){ # when X has no parents
    ll_prior <- dnorm(test_data[, pa_node], mean = mean(test_data[, pa_node]), sd = sd(test_data[, pa_node]), log = TRUE)
  }else if(class(m_prior)[1] == "lm"){
    pred_val_prior <- predict(m_prior, test_data)
    residual_se_prior <- summary(m_prior)$sigma
    ll_prior <- dnorm(test_data[,pa_node], mean = pred_val_prior, sd = residual_se_prior, log = TRUE)
  }else{
    pred_val_prior <- predict(m_prior, test_data)
    # residual_se_prior <- sqrt(m_prior$scale)
    residual_se_prior <- est_normsd(m_prior)
    ll_prior <- dnorm(test_data[,pa_node], mean = pred_val_prior, sd = residual_se_prior, log = TRUE)
  }
  return(cbind(ll_pa, ll_prior))
}



# adj_ll_test <- function(m_null, m_full){
#   df_diff_x <- get_model_df(m_null) - get_model_df(m_full)
#   test_stat <- as.numeric(-2*(logLik(m_null)-logLik(m_full)))
#   var_pval_x <- pchisq(test_stat, df=max(df_diff_x, 1), lower.tail = FALSE)
#   return(var_pval_x)
# }


#' Likelihood test for edge orientation (5/5 train test split)
#'
#' With estimated regression models, compute the log-lik under each direction
#' Execute likelihood test to determine true causal direction
#'
#' @param m_XtoY model under X->Y; Y as function of PA_Y and X
#' @param m_YtoX model under Y->X; X as function of PA_X and Y
#' @param m_x model under X->Y; X as function of PA_X
#' @param m_y model under Y->X; Y as function of PA_Y
#' @param test_data data used to estimate log-lik
#'
#' @return list containing test results; 'pval_ll_diff' is the p-value, and 'llratio' indicates edge direction
#' @examples
#' get_edge_dir(gam_XtoY[[1]], gam_YtoX[[1]], gam_XtoY[[2]], gam_YtoX[[2]], test_data = test_data)
#' @export
get_edge_dir <- function(m_XtoY, m_YtoX, m_x, m_y, test_data){
  start_time <- Sys.time()
  node_Y <- getResponseFromFormula(m_XtoY$formula)
  node_X <- getResponseFromFormula(m_YtoX$formula)
  # 1. check if nodes are linear in fitted models
  # if so for both, then relation is linear --> can't identify edge direction
  are_linear <- check_linear_relation(m_XtoY, m_YtoX, m1_response=node_Y, m2_response=node_X)
  if(are_linear){
    eval_time <- as.numeric(difftime(time1 = Sys.time(), time2 = start_time, units = "secs"))
    rval <- list(pval_ll_diff = 1, llratio = 0, eval_time = eval_time) # to indicate that edge directions are equivalent
  }
  # get log likelihoods of both directions
  ll_XtoY_ind <- gam_ll(m_XtoY, m_x, test_data, pa_node=node_X, ch_node=node_Y)
  ll_YtoX_ind <- gam_ll(m_YtoX, m_y, test_data, pa_node=node_Y, ch_node=node_X)
  ll_XtoY <- ll_XtoY_ind[,1] + ll_XtoY_ind[,2]
  ll_YtoX <- ll_YtoX_ind[,1] + ll_YtoX_ind[,2]

  # compare individual log-likelihoods
  ind_ll_diff <- ll_XtoY-ll_YtoX

  # calculate variance in log-likelihood ratios
  nmis <- sum(is.na(ll_XtoY))
  n <- NROW(ll_XtoY) - nmis
  omega.hat.2 <- var(ind_ll_diff , na.rm = TRUE)


  # if variance is very small, no need to conduct hypothesis test
  sd_mean_ratio <- abs(sqrt(omega.hat.2)/mean(ind_ll_diff, na.rm = TRUE))
  if(is.nan(sd_mean_ratio) || sd_mean_ratio >= n*1.5){
    eval_time <- as.numeric(difftime(time1 = Sys.time(), time2 = start_time, units = "secs"))
    rval <- list(omega = omega.hat.2, # estimated variance
                 pval_ll_diff = 1, # variance test p-value
                 llratio = sum(ind_ll_diff),
                 eval_time = eval_time
    )
    return(rval)
  }
  # or, if ratio of variances is small
  min_sd <- min(sqrt(var(ll_YtoX, na.rm=TRUE)), sqrt(var(ll_XtoY, na.rm=TRUE)))
  if(min_sd > 0){
    sd_ratio <- sqrt(omega.hat.2)/min_sd
    if(sd_ratio < 0.003){
      eval_time <- as.numeric(difftime(time1 = Sys.time(), time2 = start_time, units = "secs"))
      rval <- list(omega = omega.hat.2, # estimated variance
                   pval_ll_diff = 1, # variance test p-value
                   llratio = sum(ind_ll_diff),
                   eval_time = eval_time
      )
      return(rval)
    }
  }



  # run test
  ll_t_test_res <- t.test(ll_XtoY, ll_YtoX, paired=TRUE, mu=0) # t-test

  ## Calculate likelihood ratio; Eq (6.4)
  lr <- sum(ll_XtoY-ll_YtoX)
  # lr <- mean(ind_ll_diff)


  # STEP 2: testing if models fit equally
  # only valid when f!=g
  teststat <- (1/sqrt(n)) * lr/sqrt(omega.hat.2)

  if(is.nan(teststat)){
    pLRTA <- 0.5
    pLRTB <- 0.5
  }else{
    # pLRTA <- pnorm(teststat, lower.tail=FALSE)
    # pLRTB <- pnorm(teststat)
    pLRTA <- pt(teststat, df=n-1,lower.tail=FALSE)
    pLRTB <- pt(teststat, df=n-1)
  }

  eval_time <- as.numeric(difftime(time1 = Sys.time(), time2 = start_time, units = "secs"))
  rval <- list(
    mean_ll_diff = unname(ll_t_test_res$estimate), # mean of diff in individual log likelihood estimates
    stat_ll_diff = unname(ll_t_test_res$statistic), # test statistic
    sd_ll_diff = unname(ll_t_test_res$stderr), # sample error
    pval_ll_diff = unname(ll_t_test_res$p.value), # pval from t-test
    llratio = lr,  # total difference in log likelihood
    sample_var =  omega.hat.2,
    ll_values = c(sum(ll_XtoY_ind[,1]), sum(ll_XtoY_ind[,2]),
                  sum(ll_YtoX_ind[,1]), sum(ll_YtoX_ind[,2])),
    var_values = c(omega.hat.2,
                   (n-1)/n * var(ll_XtoY_ind[,1] - ll_YtoX_ind[,2], na.rm = TRUE),
                   (n-1)/n * var(ll_XtoY_ind[,2] - ll_YtoX_ind[,1], na.rm = TRUE)),
    eval_time = eval_time
  )
  return(rval)
}


# function: determining if the variance in the diff of log-liks are significant (i.e. too small); if so, can skip edge dir test and declare models as equivalent
# method: compares the ratio of the variances var(joint)/min(var(X), var(Y)) to threshold
# input: variance joint conditional loglik p(X, Y), conditional marginal loglik p(X), conditional marginal loglik p(Y)
# output: true/false is ratio if below threshold
ll_var_test <- function(joint_ll_var, ll_x, ll_y, tol=1e-3){
  min_sd <- min(sqrt(var(ll_x, na.rm=TRUE)), sqrt(var(ll_y, na.rm=TRUE)))
  # print(c(min_sd, joint_ll_var))
  if(min_sd > 0){
    sd_ratio <- sqrt(joint_ll_var)/min_sd
    if(sd_ratio < tol) return(TRUE) # variance is too small
  }
  return(FALSE)
}

#' Likelihood test for edge orientation (2fold approach)
#'
#' With estimated regression models, compute the log-lik under each direction
#' Execute likelihood test to determine true causal direction
#' Performs test twice, with train/test exchanging functions
#'
#' @param set_XtoY_1 model 1 under X->Y; Y as function of PA_Y and X, X as function of PA_X
#' @param set_YtoX_1 model 1 under Y->X; X as function of PA_X and Y, Y as function of PA_Y
#' @param set_XtoY_2 model 2 under X->Y; Y as function of PA_Y and X, X as function of PA_X
#' @param set_YtoX_2 model 2 under Y->X; X as function of PA_X and Y, Y as function of PA_Y
#' @param train_data1 data used to estimate log-lik for set 2
#' @param train_data2 data used to estimate log-lik for set 1
#'
#' @return list containing test results; 'pval_ll_diff' is the p-value, and 'llratio' indicates edge direction
#' @examples
#' get_edge_dir_2fold(gam_XtoY_m1, gam_YtoX_m1, gam_XtoY_m2, gam_YtoX_m2,
#' train_data1=train_data, train_data2=test_data)
#' @export
get_edge_dir_2fold <- function(set_XtoY_1, set_YtoX_1, set_XtoY_2, set_YtoX_2, train_data1, train_data2){
  start_time <- Sys.time()
  # saving all the models individually
  m_XtoY_1 <- set_XtoY_1[[1]]
  m_Xp_1 <- set_XtoY_1[[2]]
  m_YtoX_1 <- set_YtoX_1[[1]]
  m_Yp_1 <- set_YtoX_1[[2]]
  m_XtoY_2 <- set_XtoY_2[[1]]
  m_Xp_2 <- set_XtoY_2[[2]]
  m_YtoX_2 <- set_YtoX_2[[1]]
  m_Yp_2 <- set_YtoX_2[[2]]
  # get name of node from models
  node_Y <- getResponseFromFormula(m_XtoY_1$formula)
  node_X <- getResponseFromFormula(m_YtoX_1$formula)
  # 1. check if nodes are linear in fitted models
  # if so for both, then relation is linear --> can't identify edge direction
  are_linear <- check_linear_relation(list(m_XtoY_1, m_XtoY_2), list(m_YtoX_1, m_YtoX_2),
                                      m1_response=node_Y, m2_response=node_X)
  if(are_linear){
    eval_time <- as.numeric(difftime(time1 = Sys.time(), time2 = start_time, units = "secs"))
    rval <- list(pval_ll_diff = 1, llratio = 0, eval_time = eval_time) # to indicate that edge directions are equivalent
  }
  # get ind log likelihoods of both directions
  ll_XtoY_ind_1 <- gam_ll(m_XtoY_1, m_Xp_1, train_data2, pa_node=node_X, ch_node=node_Y)
  ll_YtoX_ind_1 <- gam_ll(m_YtoX_1, m_Yp_1, train_data2, pa_node=node_Y, ch_node=node_X)
  ll_XtoY_ind_2 <- gam_ll(m_XtoY_2, m_Xp_2, train_data1, pa_node=node_X, ch_node=node_Y)
  ll_YtoX_ind_2 <- gam_ll(m_YtoX_2, m_Yp_2, train_data1, pa_node=node_Y, ch_node=node_X)
  # get log likelihood from each each
  ll_XtoY_1 <- ll_XtoY_ind_1[,1] + ll_XtoY_ind_1[,2]
  ll_YtoX_1 <- ll_YtoX_ind_1[,1] + ll_YtoX_ind_1[,2]
  ll_XtoY_2 <- ll_XtoY_ind_2[,1] + ll_XtoY_ind_2[,2]
  ll_YtoX_2 <- ll_YtoX_ind_2[,1] + ll_YtoX_ind_2[,2]
  # combine to get loglik on full dataset
  ll_XtoY <- c(ll_XtoY_1, ll_XtoY_2)
  ll_YtoX <- c(ll_YtoX_1, ll_YtoX_2)

  # compare individual log-likelihoods
  ind_ll_diff <- ll_XtoY-ll_YtoX

  # calculate variance in log-likelihood ratios
  nmis <- sum(is.na(ll_XtoY))
  n <- NROW(ll_XtoY) - nmis
  omega.hat.2 <- (n-1)/n * var(ind_ll_diff, na.rm = TRUE)
  ll_diff_var_1 <- var(ll_XtoY_1 - ll_YtoX_1, na.rm = TRUE)
  ll_diff_var_2 <- var(ll_XtoY_2 - ll_YtoX_2, na.rm = TRUE)
  # print(c(sum(ll_XtoY_1), sum(ll_YtoX_1), sum(ll_XtoY_2), sum(ll_YtoX_2)))

  # if variance is very small, no need to conduct hypothesis test
  # or, if ratio of variances is small
  var_test_res_1 <- ll_var_test(ll_diff_var_1, ll_XtoY_1, ll_YtoX_1, tol=0.003)
  var_test_res_2 <- ll_var_test(ll_diff_var_2, ll_XtoY_2, ll_YtoX_2, tol=0.003)

  temp1 <- t.test(ll_XtoY_1, ll_YtoX_1, paired=TRUE)
  temp2 <- t.test(ll_XtoY_2, ll_YtoX_2, paired=TRUE)

  # check for 3 different conditions here; first one is for the overall variance
  lr <- sum(ll_XtoY-ll_YtoX)
  sample_mu <- lr/n
  sample_var <- (ll_diff_var_1 + ll_diff_var_2)/2
  # overall_var <- ((sqrt(ll_diff_var_1) < 0.005) + (sqrt(ll_diff_var_2) < 0.005)) > 0
  # overall_var <- abs(sqrt(sample_var)/sample_mu) >= (n/2)*1.5
  overall_var <- (abs(sqrt(ll_diff_var_1)/mean(ll_XtoY_1 - ll_YtoX_1)) >= (n/2)*1.5) + (abs(sqrt(ll_diff_var_2)/mean(ll_XtoY_2 - ll_YtoX_2)) >= (n/2)*1.5) == 2
  # print(c(abs(sqrt(ll_diff_var_1)/mean(ll_XtoY_1 - ll_YtoX_1)), abs(sqrt(ll_diff_var_2)/mean(ll_XtoY_2 - ll_YtoX_2)), n/2, overall_var,var_test_res_1, var_test_res_2))
  # if all are true, then reject
  if(overall_var | (var_test_res_1 & var_test_res_2)){
    eval_time <- as.numeric(difftime(time1 = Sys.time(), time2 = start_time, units = "secs"))
    rval <- list(omega = omega.hat.2, # estimated variance
                 pval_ll_diff = 1, # variance test p-value
                 llratio = sum(ind_ll_diff),
                 eval_time = eval_time
    )
    return(rval)
  }

  # adjusted t-test
  # accounts for merging samples w/ diff distribution parameters
  # lr <- sum(ll_XtoY-ll_YtoX)
  # sample_mu <- lr/n
  # sample_var <- (ll_diff_var_1*(length(ll_XtoY_1)-1) + ll_diff_var_2*(length(ll_XtoY_2)-1))/(length(ll_XtoY_1)+length(ll_XtoY_2)-2)
  t_hat <- (sample_mu)/sqrt(sample_var/n)
  t_test_pval <- 1 - abs(pt(t_hat, df=n-1) - pt(-t_hat, df=n-1))

  eval_time <- as.numeric(difftime(time1 = Sys.time(), time2 = start_time, units = "secs"))
  rval <- list(
    mean_ll_diff = sample_mu, # mean of diff in individual log likelihood estimates
    stat_ll_diff = t_hat, # test statistic
    sd_ll_diff = sqrt(sample_var), # sample error
    pval_ll_diff = t_test_pval, # pval from t-test
    llratio = lr,  # total difference in log likelihood
    sample_var =  omega.hat.2,
    ll_values = c(sum(ll_XtoY_1), sum(ll_YtoX_1), sum(ll_XtoY_2), sum(ll_YtoX_2)),
    var_values = c(sample_var, ll_diff_var_1, ll_diff_var_2, var(ll_XtoY), var(ll_YtoX)),
    eval_time = eval_time,
    test1 = temp1$p.value,
    test2 = temp2$p.value
  )
  return(rval)
}



#### OBSOLETE FUNCTION
# function: determining if the variance in the diff of log-liks are significant (i.e. too small); if so, can skip edge dir test and declare models as equivalent
# method: compares the ratio of the variances var(joint)/min(var(X), var(Y)) to threshold
# input: variance joint conditional loglik p(X, Y), conditional marginal loglik p(X), conditional marginal loglik p(Y)
# output: true/false is ratio if below threshold
ll_var_test <- function(joint_ll_var, ll_x, ll_y, tol=1e-3){
  min_sd <- min(sqrt(var(ll_x, na.rm=TRUE)), sqrt(var(ll_y, na.rm=TRUE)))
  # print(c(min_sd, joint_ll_var))
  if(min_sd > 0){
    sd_ratio <- sqrt(joint_ll_var)/min_sd
    if(sd_ratio < tol) return(TRUE) # variance is too small
  }
  return(FALSE)
}

# function: conducting the edge dir test using 2-folds CV with SEPARATE t-tests
# input: pretrained models and their training datasets; the training datasets serve as the other folds' test dataset
# return: p-value from t-test (and other metrics), run-time
get_edge_dir_2fold_2test <- function(set_XtoY_1, set_YtoX_1, set_XtoY_2, set_YtoX_2, train_data1, train_data2, check_both_folds=TRUE){
  start_time <- Sys.time()
  # saving all the models individually
  m_XtoY_1 <- set_XtoY_1[[1]]
  m_Xp_1 <- set_XtoY_1[[2]]
  m_YtoX_1 <- set_YtoX_1[[1]]
  m_Yp_1 <- set_YtoX_1[[2]]
  m_XtoY_2 <- set_XtoY_2[[1]]
  m_Xp_2 <- set_XtoY_2[[2]]
  m_YtoX_2 <- set_YtoX_2[[1]]
  m_Yp_2 <- set_YtoX_2[[2]]
  # get name of node from models
  node_Y <- getResponseFromFormula(m_XtoY_1$formula)
  node_X <- getResponseFromFormula(m_YtoX_1$formula)
  # 1. check if nodes are linear in fitted models
  # if so for both, then relation is linear --> can't identify edge direction
  are_linear <- check_linear_relation(list(m_XtoY_1, m_XtoY_2), list(m_YtoX_1, m_YtoX_2),
                                      m1_response=node_Y, m2_response=node_X)
  if(are_linear){
    eval_time <- as.numeric(difftime(time1 = Sys.time(), time2 = start_time, units = "secs"))
    rval <- list(pval_ll_diff = 1, llratio = sum(ind_ll_diff), eval_time = eval_time) # to indicate that edge directions are equivalent
  }
  # get ind log likelihoods of both directions
  ll_XtoY_ind_1 <- gam_ll(m_XtoY_1, m_Xp_1, train_data2, pa_node=node_X, ch_node=node_Y)
  ll_YtoX_ind_1 <- gam_ll(m_YtoX_1, m_Yp_1, train_data2, pa_node=node_Y, ch_node=node_X)
  ll_XtoY_ind_2 <- gam_ll(m_XtoY_2, m_Xp_2, train_data1, pa_node=node_X, ch_node=node_Y)
  ll_YtoX_ind_2 <- gam_ll(m_YtoX_2, m_Yp_2, train_data1, pa_node=node_Y, ch_node=node_X)
  # get log likelihood from each each
  ll_XtoY_1 <- ll_XtoY_ind_1[,1] + ll_XtoY_ind_1[,2]
  ll_YtoX_1 <- ll_YtoX_ind_1[,1] + ll_YtoX_ind_1[,2]
  ll_XtoY_2 <- ll_XtoY_ind_2[,1] + ll_XtoY_ind_2[,2]
  ll_YtoX_2 <- ll_YtoX_ind_2[,1] + ll_YtoX_ind_2[,2]
  # combine to get loglik on full dataset
  ll_XtoY <- c(ll_XtoY_1, ll_XtoY_2)
  ll_YtoX <- c(ll_YtoX_1, ll_YtoX_2)

  # compare individual log-likelihoods
  ind_ll_diff <- ll_XtoY-ll_YtoX

  # calculate variance in log-likelihood ratios
  nmis <- sum(is.na(ll_XtoY))
  n <- NROW(ll_XtoY) - nmis
  omega.hat.2 <- (n-1)/n * var(ind_ll_diff, na.rm = TRUE)
  ll_diff_var_1 <- var(ll_XtoY_1 - ll_YtoX_1, na.rm = TRUE)
  ll_diff_var_2 <- var(ll_XtoY_2 - ll_YtoX_2, na.rm = TRUE)
  # print(c(sum(ll_XtoY_1), sum(ll_YtoX_1), sum(ll_XtoY_2), sum(ll_YtoX_2)))

  # if ratio of variances is small , no need to conduct hypothesis test
  var_test_res_1 <- ll_var_test(ll_diff_var_1, ll_XtoY_1, ll_YtoX_1, tol=0.001)
  var_test_res_2 <- ll_var_test(ll_diff_var_2, ll_XtoY_2, ll_YtoX_2, tol=0.001)
  # or, if variance is very small, no need to conduct hypothesis test
  overall_var <- (abs(sqrt(ll_diff_var_1)/mean(ll_XtoY_1 - ll_YtoX_1)) >= (n/2)*1.5) + (abs(sqrt(ll_diff_var_2)/mean(ll_XtoY_2 - ll_YtoX_2)) >= (n/2)*1.5) > 0
  # print(c(abs(sqrt(ll_diff_var_1)/mean(ll_XtoY_1 - ll_YtoX_1)), abs(sqrt(ll_diff_var_2)/mean(ll_XtoY_2 - ll_YtoX_2)), n/2, overall_var, var_test_res_1, var_test_res_2))
  # if all are true, then reject
  if(overall_var | (var_test_res_1 | var_test_res_2)){
    eval_time <- as.numeric(difftime(time1 = Sys.time(), time2 = start_time, units = "secs"))
    rval <- list(omega = omega.hat.2, # estimated variance
                 pval_ll_diff = 1, # variance test p-value
                 llratio = sum(ind_ll_diff),
                 eval_time = eval_time
    )
    return(rval)
  }


  # t-test on fold 1
  ll_res_large_res <- t.test(ll_XtoY_1, ll_YtoX_1, paired=TRUE, mu=0) # t-test
  # t-test on fold 2
  ll_res_small_res <- t.test(ll_XtoY_2, ll_YtoX_2, paired=TRUE, mu=0) # t-test
  # print(ll_res_large_res)
  # print(ll_res_small_res)
  # logic: under conservative approach for both folds to be stat sig, consider the largest p-val; otherwise consider the smaller one
  # identify the larger value here
  if(ll_res_large_res$p.value < ll_res_small_res$p.value){
    ll_res_large_res_temp <- ll_res_large_res
    ll_res_large_res <- ll_res_small_res
    ll_res_small_res <- ll_res_large_res_temp
  }

  # identify which element to return
  if(check_both_folds){ # take the larger value
    ll_t_test_res <- ll_res_large_res
  }else{ # take the smaller value
    ll_t_test_res <- ll_res_small_res
  }

  eval_time <- as.numeric(difftime(time1 = Sys.time(), time2 = start_time, units = "secs"))
  rval <- list(
    mean_ll_diff = unname(ll_t_test_res$estimate), # mean of diff in individual log likelihood estimates
    stat_ll_diff = unname(ll_t_test_res$statistic), # test statistic
    sd_ll_diff = unname(ll_t_test_res$stderr), # sample error
    pval_ll_diff = unname(ll_t_test_res$p.value), # pval from t-test
    llratio = unname(ll_t_test_res$estimate)*length(ll_XtoY_1),  # total difference in log likelihood
    ll_values = c(sum(ll_XtoY_1), sum(ll_YtoX_1), sum(ll_XtoY_2), sum(ll_YtoX_2)),
    var_values = c(omega.hat.2, ll_diff_var_1, ll_diff_var_2, var(ll_XtoY), var(ll_YtoX)),
    eval_time = eval_time
  )
  return(rval)
}
