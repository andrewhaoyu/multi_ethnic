Meta = function(coef_vec,se_vec){
  var_vec = se_vec^2
  meta_var = (sum(1/var_vec))^-1
  meta_coef = meta_var*sum(coef_vec/var_vec)
  meta_se = sqrt(meta_var)
  return(c(meta_coef,meta_se))
}
#load KG.ID
setwd("/data/zhangh24/multi_ethnic/data/AABC_data")
load("BC_AFR_overall_valid_KGID.rdata")
KG.ID = sum.data.update$KG.ID
#process overall data
load("BC_AFR_overall_train.rdata")
sum.data.train = sum.data
load("BC_AFR_overall_valid.rdata")
sum.data.valid = sum.data
all.equal(sum.data.train$ALT,
          sum.data.valid$ALT)
coef_mat = cbind(sum.data.train$beta.ALT,
                 sum.data.valid$beta.ALT)
se_mat = cbind(sum.data.train$SE,
               sum.data.valid$SE)
N = nrow(sum.data.train)
coef_meta = rep(0,N)
se_meta = rep(0,N)
for(i in 1:N){
  if(i%%1000==0){print(i)}
  coef_vec = coef_mat[i,]
  se_vec = se_mat[i,]
  result.temp = Meta(coef_vec,se_vec)
  coef_meta[i] = result.temp[1]
  se_meta[i] = result.temp[2]
}

sum.data.meta = sum.data.train
sum.data.meta$beta.ALT = coef_meta
sum.data.meta$SE = se_meta
sum.data.meta$KG.ID = KG.ID
#save(sum.data.meta,file = "BC_AFR_overall_meta.rdata")
sum.data.meta$N_Cases = 6533+2702
sum.data.meta$N_Controls = 7065+3119
sum.data.meta = sum.data.meta %>% 
select(colnames(sum.data.meta)[1:13])
save(sum.data.meta,file = "BC_AFR_overall_meta.rdata")



#setwd("/data/zhangh24/multi_ethnic/data/AABC_data")
load("BC_AFR_ERpos_train.rdata")
sum.data.train = sum.data
load("BC_AFR_ERpos_valid.rdata")
sum.data.valid = sum.data
all.equal(sum.data.train$ALT,
          sum.data.valid$ALT)
coef_mat = cbind(sum.data.train$beta.ALT,
                 sum.data.valid$beta.ALT)
se_mat = cbind(sum.data.train$SE,
               sum.data.valid$SE)
N = nrow(sum.data.train)
coef_meta = rep(0,N)
se_meta = rep(0,N)
for(i in 1:N){
  if(i%%1000==0){print(i)}
  coef_vec = coef_mat[i,]
  se_vec = se_mat[i,]
  result.temp = Meta(coef_vec,se_vec)
  coef_meta[i] = result.temp[1]
  se_meta[i] = result.temp[2]
}

sum.data.meta = sum.data.train
sum.data.meta$beta.ALT = coef_meta
sum.data.meta$SE = se_meta
sum.data.meta$KG.ID = KG.ID
save(sum.data.meta,file = "BC_AFR_ERpos_meta.rdata")









load("BC_AFR_ERneg_train.rdata")
sum.data.train = sum.data
load("BC_AFR_ERneg_valid.rdata")
sum.data.valid = sum.data
all.equal(sum.data.train$ALT,
          sum.data.valid$ALT)
coef_mat = cbind(sum.data.train$beta.ALT,
                 sum.data.valid$beta.ALT)
se_mat = cbind(sum.data.train$SE,
               sum.data.valid$SE)
N = nrow(sum.data.train)
coef_meta = rep(0,N)
se_meta = rep(0,N)
for(i in 1:N){
  if(i%%10000==0){print(i)}
  coef_vec = coef_mat[i,]
  se_vec = se_mat[i,]
  result.temp = Meta(coef_vec,se_vec)
  coef_meta[i] = result.temp[1]
  se_meta[i] = result.temp[2]
}

sum.data.meta = sum.data.train
sum.data.meta$beta.ALT = coef_meta
sum.data.meta$SE = se_meta
sum.data.meta$KG.ID = KG.ID
save(sum.data.meta,file = "BC_AFR_ERneg_meta.rdata")

