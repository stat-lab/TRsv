library(xgboost)
library(Matrix)

ins = commandArgs(trailingOnly=TRUE)[1]

ins_test = commandArgs(trailingOnly=TRUE)[2]

out_ins = commandArgs(trailingOnly=TRUE)[3]

param = list(booster = "gbtree", objective = "multi:softmax", eval_metric = "mlogloss", num_class = 2, nthread = 6, eta = 0.2)


df <- read.table(ins,header=T,sep="\t")
df2 <- df[,c(1,3,4,5,6,7,8,9)]

options(na.action='na.pass')
train_dmat <- xgb.DMatrix(sparse.model.matrix(TF ~ ., data = df2), label = df2$TF)
options(na.action='na.omit')

model_ins <- xgb.train(data = train_dmat, nrounds = 30, params = param)

test_ins <- read.table(ins_test,header=T,sep="\t")
test_ins2 <- test_ins[,c(1,3,4,5,6,7,8,9)]
test_ins_dmat <- xgb.DMatrix(sparse.model.matrix(TF ~ ., data = test_ins2), label = test_ins2$TF)
pred_ins <- predict(model_ins, test_ins_dmat)

write.table(pred_ins, out_ins, quote=F, col.names=F)
