library(xgboost)
library(Matrix)

del = commandArgs(trailingOnly=TRUE)[1]

del_test = commandArgs(trailingOnly=TRUE)[2]

out_del = commandArgs(trailingOnly=TRUE)[3]

param = list(booster = "gbtree", objective = "multi:softmax", eval_metric = "mlogloss", num_class = 2, nthread = 6, eta = 0.2)


df <- read.table(del,header=T,sep="\t")
df2 <- df[,c(1,3,4,5,6,7,8,9)]

options(na.action='na.pass')
train_dmat <- xgb.DMatrix(sparse.model.matrix(TF ~ ., data = df2), label = df2$TF)
options(na.action='na.omit')

model_del <- xgb.train(data = train_dmat, nrounds = 100, params = param)

test_del <- read.table(del_test,header=T,sep="\t")
test_del2 <- test_del[,c(1,3,4,5,6,7,8,9)]
test_del_dmat <- xgb.DMatrix(sparse.model.matrix(TF ~ ., data = test_del2), label = test_del2$TF)
pred_del <- predict(model_del, test_del_dmat)

write.table(pred_del, out_del, quote=F, col.names=F)
