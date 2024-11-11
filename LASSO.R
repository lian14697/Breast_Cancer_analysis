# -------------- ★★★ Final-LASSO 2023-4-25--------------- #
# dev.off()
cat("\014"); rm(list = ls());  library(glmnet); library(caTools);library(caret); library(ROCR)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); seed1 = 0


# ===========================================================================================
# ------------- 1) Divide the data set randomly -------------------
mydata   <- read.csv("ComBat_ICC_linear_source_All.csv"); ex = 0; names(mydata)[1]<-"p_ID"
mydata <- mydata[c(which(mydata$batch != 3)), ]
# val_method = 'load'; val_list = c("Combat_data2.csv"); ex = 1
val_method = 'Define'; val_list = c('DUKE'); ex=1
seed1 = 725; split = sample.split(mydata$label, SplitRatio = 0.7,set.seed(seed1))
train_data = subset(mydata, split == TRUE); test_data = subset(mydata, split == FALSE)


cc = 17; Uni_p = 0.05; corr=0.9; alpha=1

# ## ------ Training set data upsampled --------------------------
Train_data = train_data
if(sum(Train_data$label) > nrow(Train_data)/2){train_data$label = 1-train_data$label}
# --- Calculated usage ---
before = sum(train_data$label)/length(train_data$label);  uptime = round(sum(train_data$label==0)/sum(train_data$label))-1
# --- Perform upsampling ---
if (uptime>0) {Yang = train_data[train_data$label == 1,]; for (up in 1:uptime) {train_data = rbind(train_data, Yang)}; after = sum(train_data$label)/length(train_data$label); before; after; uptime}
# --- Yin-yang inversion ---
if(sum(Train_data$label) > nrow(Train_data)/2){train_data$label = 1-train_data$label}
# ------------ Batch - single factor rank sum test --------------
Uni_list = data.frame()
Index0 = which(train_data$label==0);  Index1 = which(train_data$label==1)
for (j in ((cc+1):ncol(train_data))){Uni_list[j-cc, 1] = names(train_data)[j]
Uni_list[j-cc, 2] = wilcox.test(train_data[Index0, j], train_data[Index1, j], alternative = "two.sided")$statistic
Uni_list[j-cc, 3] = wilcox.test(train_data[Index0, j], train_data[Index1, j], alternative = "two.sided")[["p.value"]]}
names(Uni_list) = c('Feature', 'Statistic', 'P.value')

Sig.Num = length(which(Uni_list$P.value < Uni_p))
print(paste0('The number of significant features of the single factor rank sum test is:', Sig.Num))
Sig_uni_list = Uni_list[which(Uni_list$P.value < Uni_p), ]
train_data = train_data[, c(1:cc, cc+as.numeric(row.names(Sig_uni_list)))]




# # ----------- pearson Correlation analysis ---------------
# Feature pairs with |r| > 0.6 were selected, and then in each of these pairs, the feature with larger mean absolute correlation was removed.
x = train_data[,(cc+1):ncol(train_data)]
highcorr <- findCorrelation(cor(x),cutoff = corr,names = FALSE)
if (length(highcorr)==0){highcorr=0} 
train_data <- train_data[, -(cc+highcorr)]
Corr.num=ncol(train_data)-cc; print(paste0('Correlation analysis of residual features:', Corr.num))





# ------------- Data standardization (x-mean(x))/std(x) --------------------
y = train_data$label; x = scale(train_data[, -c(1:cc)], center = TRUE, scale = TRUE)


# ------------- LASSO-logistic --------------
fit=glmnet(x,y,family='binomial',type.measure="auc", alpha=alpha, intercept = FALSE)
pdf("min.pdf");    myfit2 <- cv.glmnet(x,y,family='binomial', type.measure="auc", alpha=alpha, intercept = FALSE, nfolds=5, set.seed(100));  plot(myfit2); abline(v=log(c(myfit2$lambda.min, myfit2$lambda.1se)),lty="dashed"); print(myfit2); dev.off()
pdf("lambda.pdf"); plot(fit, xvar = "lambda", label = TRUE); abline(v=log(c(myfit2$lambda.min,myfit2$lambda.1se)),lty="dashed"); dev.off()
# ------------- Deterministic truncation parameter -------------
log_lambda = myfit2; coe <- coef(fit, s = myfit2$lambda.min); act_index <- which(coe != 0); act_coe <- coe[act_index]; num = length(row.names(coe)[act_index]); lasso.method = 'min'
if (num > 20){coe <- coef(fit, s = myfit2$lambda.1se); act_index <- which(coe != 0); act_coe <- coe[act_index]; num = length(row.names(coe)[act_index]); lasso.method = '1se'}
# ------------- Determine standardized parameters -----------
LASSO_data = train_data[row.names(coe)[act_index][1:num]] 
Scale_LASSO = scale(LASSO_data, center=TRUE,scale=TRUE)


# ------------ Finishing characteristic coefficient -------------------
library(data.table)
get_coe <- function(the_fit,the_lamb){Coefficients <- coef(the_fit, s = the_lamb);Active.Index <- which(Coefficients != 0); Active.Coefficients <- Coefficients[Active.Index]; re <- data.frame(rownames(Coefficients)[Active.Index],Active.Coefficients); re <- data.table('Features'=rownames(Coefficients)[Active.Index], 'Coefficient'=Active.Coefficients); re$OR <- exp(re$Coefficient); return(re[order(OR)])}
if(lasso.method == 'min') {log_lambda = log(myfit2$lambda.min); results = get_coe(myfit2,myfit2$lambda.min)}
if(lasso.method == '1se') {log_lambda = log(myfit2$lambda.1se); results = get_coe(myfit2,myfit2$lambda.1se)}

info = paste0("seed1=", seed1, ', Up=', uptime, ', method=',lasso.method, ", log(λ)=",round(log_lambda,3), ', corr=',corr,', α=',alpha, ', num=', num)
write.csv(results, file = paste0(info, ".csv"), row.names =FALSE, quote =TRUE)



# ============ Model application =======================
model_app = function(data, cohort){
  # ------------------ radscore ------------------
  LASSO_data = data[row.names(coe)[act_index][1:num]] 
  scale_LASSO = scale(LASSO_data, center = attr(Scale_LASSO,'scaled:center'), scale = attr(Scale_LASSO,'scaled:scale'))
  radscore = t(act_coe[1:num] %*% t(scale_LASSO)) 

  # ----------------------- AUC ------------------------
  library(ROCR); pred <- prediction(1/(1+exp(-radscore)), data$label);
  auc <- round(unlist(slot(performance(pred,'auc'),'y.values')), 3)
  
  ## ---------------------- Calibration curve -----------------------
  pdf(paste0('AUC = ', auc, ",Calibration curve ", cohort, ".pdf"));  pacman::p_load("rms"); y_pred <- 1/(1+exp(-radscore)); val.prob(y_pred,as.numeric(data$label),m=10,cex=0.5); dev.off()
  ## ---------------------- Fall chart -----------------------
  library(ggplot2); patient = 1:length(radscore); pp <- data.frame(radscore, label=data$label); pp = pp[order(pp$radscore),]; pp <- data.frame(patient, pp)
  pdf(paste0("Waterfall plot ", cohort, ".pdf"), height=3, width=8); Group <- as.factor(pp$label); ggplot(pp, aes(x=patient, y=radscore,fill=Group)) + geom_bar(stat="identity") + labs(x="Patients", y="Radscore")+theme(panel.background = element_rect(fill = "transparent",colour = NA), panel.border = element_rect(linetype = "solid", fill = NA)); dev.off()
  
  # ----------------Result saving--------------------
  df1 <- cbind(data[, 1:cc], radscore);   write.csv(df1, file = paste0('LASSO得分_', cohort, '.csv'), row.names =FALSE)
  df2 <- cbind(data[, 1:cc], LASSO_data); write.csv(df2, file = paste0('LASSO特征_', cohort, '.csv'), row.names =FALSE)
  return(paste0(cohort, ' AUC = ', auc))
}


## ------------------------ External verification -----------------------------
if (ex==1){for (k in c(1:length(val_list))) {
  if (val_method == 'load')   {val_data = read.csv(val_list[k])}
  if (val_method == 'Define') {val_data = mydata[c(which(mydata$group==val_list[k])), ]}
  print(model_app(val_data, val_list[k]))}
}
## ------------------------ Training verification -----------------------------
model_app(train_data, 'train'); model_app(test_data,  'test')


info; print('finish')
