# -------------- ZEN-ComBat BY:Shelby 2023-4-5---------------- #
# dev.off()
cat("\014"); rm(list = ls()); options(warn = -1); options(digits=3); cc0=0
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); library(sva)


# ============================== 1.Read data ====================================================
# ---------------------- 1.Data sorting -----------------------------
filename0 = "ICC_linear_source_All.csv"; cc = 17
Groupname = c('ZJPH', 'GDPH', 'DUKE'); Batch = c(1, 2, 3)
mydata <- read.csv(filename0); names(mydata)[1]<-"p_ID"; row.names(mydata) = mydata$p_ID; head(mydata[1:3, 1:5])
#     p_ID label           group   batch  f1
# 1   P1     1   Hospital I 1.5T     1    5644


# ---------------- Batch sort--------------------
mydata = data.frame(mydata[order(mydata$batch), ])
# ---------------- Parameter preparation --------------------
Hospital = unique(mydata$group); Hospital
Color = c("#0bb6c3", "#b394e5", "#efd15a", "#dc5c27", "#fcbea5", "#f17f4f","#d9af09" )



# ============================= 2.Data preprocessing =====================================================
# ------------- Preprocessing: deletes features that duplicate 60% ----
ref.batch = 3; cc = ifelse(cc0>0, cc0+2, cc); V = c()
ref.data = mydata[which(mydata$batch==ref.batch), ]
for (v in c((cc+1):ncol(ref.data))) {Num = length(which(duplicated(ref.data[, v])))
if (Num > 0.6*nrow(ref.data) || length(which(is.na(mydata[, v])))>0) {V = c(V, v)}}; if (is.null(V) == FALSE) {mydata = mydata[, -V]}



# # ----------- Preprocessing: Improved - Remove header and tail 2.5%~97.5% (3/24 fixed Low/High string problem)--------------------
Low = as.numeric(); High = as.numeric()
ref.data = mydata[which(mydata$batch==ref.batch), ]
for (j in c((cc+1):ncol(mydata))) {
  Low[j-cc] = quantile(ref.data[,j], 0.025)
  index05 = which(ref.data[,j] <= min(ref.data[ref.data[,j]>Low[j-cc], j]))
  ref.data[index05, j]= min(ref.data[ref.data[,j]>Low[j-cc], j])
  
  High[j-cc] = quantile(ref.data[,j],0.975)
  index95 = which(ref.data[,j] >= max(ref.data[ref.data[,j]<High[j-cc],j]))
  ref.data[index95, j]= max(ref.data[ref.data[,j]<High[j-cc],j])}
mydata[row.names(ref.data), ] = ref.data



# ------------- Data standardization--------------------
Mean = as.numeric(); SD = as.numeric()
ref.data = mydata[which(mydata$batch==ref.batch), ]
for (v in c((cc+1):ncol(mydata))) {
  # Z-score standardization
  Mean[v-cc] = mean(ref.data[,v]); SD[v-cc] = sd(mydata[,v]);
  # Max[v-cc] = max(ref.data[,v]); Min[v-cc] = min(mydata[,v])
  for (j in c(1:length(Groupname))) {
    Index = which(mydata$group==Groupname[j])
    mydata[Index, v] = (mydata[Index, v]-Mean[v-cc])/SD[v-cc]
    # mydata[Index, v] = (mydata[Index, v]-Min[v-cc])/(Max[v-cc]-Min[v-cc])
  }
}


# ================================ 3.ComBat ========================================
total_x = mydata[, -c(1:cc)]; combat_x = ComBat(dat=t(total_x), batch=mydata$batch, mod=NULL, ref.batch=ref.batch, par.prior=TRUE, prior.plots=FALSE); total_x <- t(combat_x)

# -------------- Result save --------------------------
ComBat.data = cbind(mydata[, c(1:cc)], total_x)
if (length(filename0)==1){write.csv(ComBat.data, paste0('ComBat_', filename0), row.names = F)}
if (length(filename0)>=2){write.csv(ComBat.data, paste0('ComBat_Mixdata.csv'), row.names = F)}




# ================================= 4.PCA =============================================
library(ggplot2); library(factoextra); library(FactoMineR)
# ---- PCA calculation: before correction -----------
Bef.Data = mydata
Bef.pca = PCA(Bef.Data[,-c(1:cc)], graph=FALSE)

# ---- PCA calculation: after correction-----------
Af.Data  = ComBat.data
Af.pca  = PCA( Af.Data[,-c(1:cc)], graph=FALSE)

lim = 100
# ---------------- PCA drawing -------------------
pdf('Before_ComBat.pdf', width=4, height = 4)
fviz_pca_ind(Bef.pca, geom.ind="point", pointsize=3, pointshape=21, 
             fill.ind = as.factor(Bef.Data$group),
             palette = Color, addEllipses=TRUE, title = 'Before batch calibration',
             legend.title="")+ 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  theme(legend.position='none', plot.title = element_text(hjust = 0.5))+
  xlim(-lim, lim) + ylim(-lim, lim) #+
dev.off()

pdf('After_ComBat.pdf', width=5, height = 4)
fviz_pca_ind(Af.pca, geom.ind="point", pointsize=3, pointshape=21, 
             fill.ind = as.factor(Af.Data$group),
             palette = Color, addEllipses=TRUE, title = 'After batch calibration',legend.title="Group")+ 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(-lim, lim) + ylim(-lim, lim)
dev.off()

print('finish')
