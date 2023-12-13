#===========================
#Author:Chen Shan
#First Creation Date:2023-11-25
#Last Modify Date:2023-12-4
#===========================
rm(list = ls())

sample.fun=function(DT,GT,trait,idx,sampling=TRUE,transform=FALSE,percentage,cycle.number){
  
  pheno.f=DT
  #Phenotype file processing, extracting material columns and trait value columns
  
  pheno.f=pheno.f[!is.na(pheno.f[2]),]
  names(pheno.f)[1]="Genotype"#表型材料列名为Genotype，和后续代码一致
  names(pheno.f)[2]=trait
  #trait=names(pheno.f)[2]
  pheno.f <- pheno.f[!duplicated(pheno.f$Genotype),]#表型删除重复的材料
  
  geno.f=GT
  names(geno.f)[1]="Genotype"
  geno.f.colnames <- geno.f$Genotype#取基因型的材料名
  dim(geno.f)
  geno.f.colnames <- geno.f.colnames[!duplicated(geno.f.colnames)]# delete repeat (materials names)
  #geno.f <- geno.f[geno.f.colnames,]
  
  
  #Intersection of genotype material names and phenotype material names
  all_pheno_geno <- intersect(pheno.f$Genotype,geno.f.colnames)#表型的材料名和基因型的材料名取交集
  index <- geno.f.colnames %in% all_pheno_geno#基因型材料名只取共有材料名
  
  geno.f <- geno.f[index,]#基因型数据
  
  index <- pheno.f$Genotype %in% all_pheno_geno#表型材料名只取共有材料名
  pheno.f <- pheno.f[index,]
  
  #保存对齐的基因型和表型文件
  write.table(pheno.f,paste0(idx,"_",trait,"_pheno_final.txt"),sep="\t",quote=F,row.names=F)
  write.table(geno.f,paste0(idx,"_",trait,"_geno_final.txt"),sep="\t",quote=F,row.names=F,col.names=T)
  
  
  if(sampling==TRUE){
    #sampling
    set.seed(123)
    percentage <- as.numeric(percentage)
    cycle.number <- as.numeric(cycle.number)
    
    pheno.f.nrow <- nrow(pheno.f)#表型的行数
    index.coefficient <- round(pheno.f.nrow*percentage)#训练群体的个数
    
    if (!dir.exists(paste0("percentage_",percentage))){
      dir.create(paste0("percentage_",percentage))
    }#新建percentage_文件夹
    
    for (i in 1:cycle.number){
      if (!dir.exists(paste0("percentage_",percentage,"/",i))){
        dir.create(paste0("percentage_",percentage,"/",i))
      }#新建每个循环的文件夹
      pheno.f.new <- pheno.f[sample(nrow(pheno.f),index.coefficient),]#抽样训练群体
      write.table(pheno.f.new,paste0("percentage_",percentage,"/",i,"/",idx,"_",trait,"_",percentage,"_",i,".txt"),sep="\t",quote=F,row.names=F)
      #保存每一次抽样的结果
    }
  }
  
  if(transform==TRUE){
    # ###对基因型文件进行处理#####
    
    Gx<-as.data.frame(geno.f[,1])
    
    Gy<-geno.f[1,-1]
    
    Gx<-as.data.frame(geno.f[,1])
    
    colnames(Gx)="Genotype"
    
    GD1<-geno.f[,-1]#去第一列材料名
    
    GD1 <- GD1-1 #基因型转为-1，0，1
    
    GD1 <- as.data.frame(GD1)
    
    GD2 <- cbind(Gx,GD1)#
    
    write.table(GD2,paste0(idx,"_",trait,"_GD_final.txt"),sep="\t",quote=F,row.names=F,col.names=T)
  }
}

library(tidyverse)

setwd("/Output") #

#read files
cat("Please choose a Phenotype file.\n")

pheno.fa <- read.csv("../datafiles/DT.csv",header=TRUE)#read phenotype file

cat("Please choose a genotype file. \n")

#Read the genotype file with column names are SNPs and row names are materials
geno.fa <- read.table("../datafiles/GT.txt",header=TRUE,check.names = F)

#表型文件处理，取单个性状
dim(pheno.fa)

unique(pheno.fa$Year)

pheno.fa=pheno.fa %>% filter(Year %in% "2020")#取不同年份的数据

idx="2020"
  
traits=unique(pheno.fa$Trait) #

#基因型文件处理
GT=geno.fa[,-c(1,3:6)]

for(trait in traits){
  
  DT=pheno.fa %>% filter(Trait %in% trait)
  
  DT=DT[,c("Name","observation_value")]
  
  sample.fun(DT=DT,GT=GT,trait,idx,sampling=TRUE,transform=FALSE,percentage=0.8,cycle.number=1)

}
##########预测########

setwd(paste0("percentage_",percentage))#设置工作路径为抽样得到的文件夹中

library(rrBLUP)
library(sommer)
library(randomForest)

cat("Such as: C:/Users/Dell/DH2_percentage_0.8")

models=c("rrBLUP_rrBLUP","GBLUP","randomForest")

accuracy.all=NULL

for (trait in traits) {
  all_folder <- sort(as.numeric(dir(pattern="[0-9]")))
  gf <- dir(path="../",pattern = paste0(idx,"_",trait,"_geno_final.txt"))   # 在上级目录（"../"）读取基因型文件(.NumericGenoForGS.csv)
  
  #gf <- dir(path="../",pattern = paste0(trait,"_GD_final.txt")   # 在上级目录（"../"）读取基因型文件(.NumericGenoForGS.csv)
  
  pf <- dir(path="../",pattern =paste0(idx,"_",trait,"_pheno_final.txt"))   # 读取表型文件
  myGeno <- read.delim(paste0("../",gf),header=T,check.names = F)
  myPheno <- read.delim(paste0("../",pf),header=T)
  
  assistGenoCol1 <- as.vector(myGeno[,1])   # 获取基因型名称
  assistPhenoCol1 <- as.vector(myPheno[,1])   # 获取表型名称
  
  
  #idx=substr(myPheno$Genotype[1],1,5)
  
  PRE_PHE.all=NULL
  accuracy = matrix(nrow=cycle.number, ncol=length(models)+2)
  colnames(accuracy) <- c("Trait", models,"geno")
  
  
  for (ix in 1:max(all_folder)){
    
    plf <- dir(path=as.character(ix), pattern = paste0(idx,"_",trait,"_",percentage))   # 读取列表文件
    
    myList <- read.delim(paste0(as.character(ix),"/",plf),header=T) #每个文件夹下训练群体表型文件
    
    resultNames <- intersect(assistGenoCol1,assistPhenoCol1)   # 基因型名和表型名取交集
    length(assistGenoCol1)   # 基因型
    length(assistPhenoCol1)   # 表型
    length(myList$Genotype)   # 列表
    length(resultNames)   # 交集
    pop <-2
    assistCol <- cbind(resultNames,pop)   # 获得辅助列
    colnames(assistCol) <- c("Genotype","pop")
    index <- assistCol[,1] %in% myList$Genotype #找到训练群体材料名
    assistCol[index,2] <- "1"   # 训练群体的pop改为1
    
    #write.csv(assistCol,paste0(as.character(ix),"/finalTaxa.csv"),row.names = F)   # 生成文件，手动编号
    #assistCol <- read.csv(paste0(as.character(ix),"/finalTaxa.csv"),header=T)   # 读取编号后的文件
    
    colnames(assistCol)[1] <- names(myGeno)[1]   # 统一材料列名
    names(myPheno)[1] <- names(myGeno)[1]   # 统一材料列名
    
    bigFile <- merge(assistCol,myPheno)   # 合并数据框 表型
    bigFile <- merge(bigFile,myGeno)   # 合并数据框 基因型和表型
    rownames(bigFile)=bigFile$Genotype
    
    pop1 <- which(bigFile$pop==1)   
    pop2 <- which(bigFile$pop==2)   
    
    
    traits.tmp <- dim(myPheno)[2]-1 #表型文件包括第一列材料名，第二列及以后性状，看有多少个性状
    Markers_impute <- bigFile[,-c(1:(2+traits.tmp))]#取基因型文件
    Markers_impute<-as.matrix(Markers_impute)
    
    #impute = A.mat(Markers,impute.method="EM",return.imputed=T)   #### A.mat函数
    #该函数主要用来过滤和填充基因型，返回加性效应关系矩阵（即Kinship）
    
    #str(impute)
    # Markers_impute <- impute$imputed
    # Markers_impute[1:6,1:6]
    # dim(Markers_impute)
    # write.table(Markers_impute, paste(as.character(ix),"/","Geno", ".Imputed.txt",sep=""),
    #             quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
    #Markers_impute=Markers
    
    Pheno <-bigFile[,c(1,3:(2+traits.tmp))]  # 第4列，表型列，多表型预留位置
    #colnames(Pheno) <- names(bigFile)[3:(2+traits)]   # 列名
    head(Pheno)
    dim(Pheno)
    #nS <- dim(Markers_impute)[1]   # 混合群体数量
    
    train <- pop1   # 建模群体
    test <- pop2   # 验证群体
    #accuracy = matrix(nrow=1, 2) #accuracy的大小，行 列
    
    PRE_PHE=NULL
    PRE_PHE=Pheno[test,]  
    # PRE_PHE$
    # output_file <- paste0(as.character(ix), "/", gf, "_combined_predictions.csv")
    
    #write.csv(PRE_PHE, output_file, row.names = FALSE)
    
    #######-----------rrBLUP-RRBLUP-------------------------------------- 
    if("rrBLUP_rrBLUP" %in% models ){
      model1="rrBLUP_rrBLUP"
      cat (paste("processing",trait,model1,ix,"......\n",sep=" "))
      
      
      Pheno_train=as.matrix(Pheno[train,2]) #建模群体的表型
      m_train=Markers_impute[train,] #建模群体的基因型
      Pheno_valid=as.matrix(Pheno[test,2]) #预测群体的表型
      m_valid=Markers_impute[test,] #预测群体的基因型
      yield <- as.matrix(as.numeric(Pheno_train[,1]))
      yield_answer<-mixed.solve(yield, #观测值
                                Z=m_train, #随机效应矩阵
                                K=NULL, #随机效应协变量矩阵
                                #X=NULL,#固定效应矩阵
                                #method="REML",  #最大似然估计方法，ML/REML
                                SE = FALSE, #是否计算标准误
                                return.Hinv=FALSE)#是否H取逆，一般在GWAS中用，忽略之
      
      
      
      # u$mean=colMeans(Markers_impute)/2#列求平均，计算基因频率，基因型编码为012，
      # range(u$mean)
      # 
      # u$Va=2*u$mean*(1-u$mean)*u$V1^2
      # u=u[order(-u$Va),]
      # plot(u$Va)
      # u$cummu=cumsum(u$Va)
      # plot(u$cummu)
      #write.csv(u,"marker_effect.csv")
      u = as.matrix(yield_answer$u)
      pred_yield_valid =  m_valid %*% u #预测群体的基因型*标记效应=估计育种值
      
      pred_yield_rrblup = as.data.frame(as.numeric(pred_yield_valid[, 1]) + as.vector(yield_answer$beta)) #预测预测群体的表型值,加了群体平均
      names(pred_yield_rrblup)="pred_rrblup"
      PRE_PHE=cbind(PRE_PHE,pred_yield_rrblup)
      
      yield_valid = as.matrix(as.numeric(Pheno_valid[,1])) #实际预测群体的表型值
      accuracy[ix,model1] <-cor(pred_yield_valid, yield_valid, use="complete" ) #相关系数
    }
    #######-----------GBLUP_sommer--------------------------------------
    if("GBLUP" %in% models ){
      model2="GBLUP"
      cat (paste("processing",trait,model2,ix,"......\n",sep=" "))
      
      # GBLUP pedigree-based approach
      #pheno有两列，第一列为材料名，第二列为性状值
      y.trn=Pheno
      colnames(y.trn)=c("id","X1")
      y.trn$id=as.factor(y.trn$id)
      y.trn[test,"X1"] <- NA#预测群体的表型值设为NA
      
      if(ix==1){
        A <- sommer::A.mat(Markers_impute)
        saveRDS(A,"../Amat.rds")
      } else{
        A <- readRDS("../Amat.rds")
      }
      
      #test first trait X1
      system.time(
        ans <- mmer(X1~1, 
                    random=~vsr(id,Gu=A),
                    rcov=~units,
                    data=y.trn,verbose = FALSE) # kinship based
      )
      #summary(ans)$varcomp  #查看方差组分
      #vpredict(ans,h2~V1/(V1+V2))  #计算遗传力
      ansU <- as.data.frame(ans$U$`u:id`$X1)#所有材料的预测值
      rownames(ansU) <- gsub("id","",rownames(ansU))
      
      pred_yield_gblup = as.data.frame(c(ansU[, 1]) + as.vector(ans$Beta$Estimate))#加群体平均，预测的表型值
      
      PRE_PHE$pred_gblup <- pred_yield_gblup[test,]
      
      accuracy[ix,model2]<-cor(ansU[test,],Pheno[test,2], use="complete")
    }
    #######-----------rrBLUP-sommer-------------------------------------- 
    if("rrBLUP_sommer" %in% models){
      cat (paste("processing",trait,model3,ix,"......\n",sep=" "))
      
      model3="rrBLUP_sommer"
      y.trn=Pheno
      colnames(y.trn)=c("id","X1")
      y.trn$id=as.factor(y.trn$id)
      y.trn[test,"X1"] <- NA#预测群体的表型值设为NA
      
      system.time(
        ans2 <- mmer(X1~1,
                     random=~vsr(list(Markers_impute)),
                     rcov=~units,
                     data=y.trn,verbose = FALSE) # kinship based
      )
      u <- Markers_impute %*% as.matrix(ans2$U$`u:Markers_impute`$X1) # BLUPs for individuals
      rownames(u) <- rownames(Markers_impute)
      
      accuracy[ix,model3]<-cor(u[test,],Pheno[test,2], use="complete")# same correlation
      
      pred_yield_rrblup_sommer = as.data.frame(c(u[, 1]) + as.vector(ans2$Beta$Estimate))#加群体平均，预测的表型值
      names(pred_yield_rrblup_sommer)="pred_rrblup_sommer"
      PRE_PHE$pred_rrblup_sommer <- pred_yield_rrblup_sommer[test,]
      
      
    }
    #######-----------randomForest-------------------------------------- 
    
    if("randomForest" %in% models ){
      model4="randomForest"
      
      bigFile=bigFile[,c(-1,-2)]
      names(bigFile)[1]="trait"
      
      rf_ntree<- randomForest(trait ~ ., data=bigFile[train,],  
                              ntree=200,important=TRUE,proximity=TRUE)#除了trait以外的所有列都作为特征，data为训练数据集。
      plot(rf_ntree)#画出错误率与决策树的数量的关系曲线
      rf_pred<-predict(rf_ntree, newdata=bigFile[test,])
      
      PRE_PHE$pred_rf <- rf_pred
      
      accuracy[ix,model4]<-cor(rf_pred,Pheno[test,2], use="complete")
      
      #varImpPlot(rf_ntree,n.var=20)#看变量的重要性
    }
    #预测精度
    accuracy[ix,"Trait"]=paste0(idx,"_",trait)
    accuracy[ix,"geno"]=gf
    rownames(accuracy)<- paste0("Cyc_",1:cycle.number)
    accuracy=as.data.frame(accuracy)
    accuracy$cyc=rownames(accuracy)
    
    
    #####结果处理######
    #预测的表型值输出
    #names(PRE_PHE)=c("Genotype","valid","pred_rrblup","pred_gblup")
    PRE_PHE$trait=trait
    PRE_PHE$geno=gf
    PRE_PHE$cycle=ix
    
    PRE_PHE.all=rbind(PRE_PHE.all,PRE_PHE)
  }
  
  write.csv(PRE_PHE.all,paste0(gf,"_prediction",".csv",sep=''),row.names=F)
  
  accuracy.all=rbind(accuracy.all,accuracy)
}

#预测精度输出

# idx=paste0(cat,"_",trait)
# column_names <- c(idx,"rrBLUP-RRBLUP","GBLUP")
write.table(accuracy.all,paste0(idx,"_acc_resultList.csv"),row.names=F,sep=",")


