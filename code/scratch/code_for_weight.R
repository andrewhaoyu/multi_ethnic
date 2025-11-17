#choose the prediction algorithms
SL.libray <- c(
  "SL.glmnet",
  "SL.ridge"
)
sl = SuperLearner(Y = y_tun, X = prs_tun_clean, family = gaussian(),
                  # For a real analysis we would use V = 10.
                  # V = 3,
                  SL.library = SL.libray)
sl_fit = sl
#algorithm weight
alg_weights <- sl_fit$coef
#glmnet
glmnet_obj = sl_fit$fitLibrary$SL.glmnet$object
best_lambda <- glmnet_obj$lambda[which.min(glmnet_obj$cvm)]
glmnet_coefs <- coef(glmnet_obj, s = best_lambda)
#ridge
ridge_coefs = sl_fit$fitLibrary$SL.ridge_All$bestCoef
#final
final_coefs <- alg_weights[1] * glmnet_coefs + alg_weights[2] * ridge_coefs

