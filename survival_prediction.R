
#Script one
#Survival function prediction for the common DEGs

# 1. Set your work folder
setwd(".../workfolder")


# 2. Load libraries for script one
library(dplyr)
library(plyr)
library(survival)

#3. Working with clinical data
dato<-read.csv("data_clinical_patient_CRC.csv",header = TRUE,stringsAsFactors = FALSE)

View(dato)

datos<-select(dato,Patient_ID,AGE,SEX,Cancer_STAGE,ETHNICITY,RACE,Censor_Status,Cancer_type,Tumor_site,OS_DAYS)
View(datos)

#ranme the column
datos<-rename(datos,Rectime=OS_DAYS,Cancer_stage=Cancer_STAGE,Race=RACE,Censor.Status=Censor_Status,Ethnicity=ETHNICITY)

datos$race=as.character(lapply(datos$race,function(x){gsub("^$","Others",x)}))
datos$tumour_site=as.character(lapply(datos$tumour_site,function(x){gsub("^$",NA,x)}))
datos$cancer_stage=as.character(lapply(datos$cancer_stage,function(x){gsub("^$",NA,x)}))
datos$anatomic_site=as.character(lapply(datos$anatomic_site,function(x){gsub("^$",NA,x)}))
datos$histologic_grade=as.character(lapply(datos$histologic_grade,function(x){gsub("^$",NA,x)}))

#4. View the distribution of the clinical features.

count(datos,'Cancer_stage')
count(datos,'Race')
count(datos,'SEX')
count(datos,'Tumor_site')
count(datos,'Ethnicity')
count(datos,'Cancer_type')


#5. Working with the  mRNA data

expr_mr<-read.csv("Genes_Zscore_CRC.csv",header=TRUE,stringsAsFactors = FALSE)
View(expr_mr)

#Replace missing values (NA) with 0
expr_mr[is.na(expr_mr)] <- 0

#function for labelling each expression value
altered_test<-function(x){
  
  if(typeof(x)=="character"){
    d=x
  }
  else{
    
    if (abs(x)>=1.5){
      d="Altered"
      
    }
    else{
      d="Normal"
      
    }
  }
  d
}

applyfunc<-function(df,f){
  ds<-matrix(0,nrow = nrow(df),ncol=ncol(df))
  colnames(ds)<-colnames(df)
  for (i in seq(1:ncol(df))){
    ds[,i]<-(sapply(expr_mr[,i],f))
  }
  ds<-as.data.frame(ds)
}
gene_status<-applyfunc(expr_mr,altered_test)

#remove the 01 from patient iD

remove_01<-function(x){
  x<-unlist(strsplit(x,split=""))
  x<-paste(x[0:(length(x)-3)],collapse = "")
  x
}
gene_status$Patient_ID<-as.character(gene_status$Patient_ID)
gene_status$Patient_ID=unlist(lapply(gene_status$Patient_ID,remove_01))


# 6. Merge the tables
gene_status$Patient_ID=as.character(gene_status$Patient_ID)
combined<-datos%>%inner_join(gene_status)

#relevel the genes as normal as reference factor
applyrevel<-function(combined){
  
  col_names<-colnames(combined)[11:ncol(combined)]
  for(i in col_names){
    combined[,i]<-as.factor(combined[,i])
    combined[,i]<-relevel(combined[,i],ref="Normal")
  }
  combined
}
combined<-applyrevel(combined)

#7. Univariate analysis

kmsurvo<-Surv(combined$Rectime,combined$Censor.Status)
applycox<-function(combined){
  models<-list()
  col_names<-colnames(combined)[11:ncol(combined)]
  for(i in col_names){
    
    fit<-coxph(kmsurvo~factor(combined[,i]),data=combined)
    tss<-summary(fit)
    coefs<-c(tss$coefficients[c(1,2,5)])
    models[[i]]=coefs
  }
  
  final_mode<-as.data.frame(models) 
  final_model=t(final_mode)
  colnames(final_model)<-c("coef","exp.coef","p")
  
  as.data.frame(final_model)
}
fs<-applycox(combined)

View(fs)
fs<-fs%>%mutate(gene=rownames(.))
write.csv(fs, file="R_unique_univariate_CRC.csv")


#8. Multivariate Analysis

fitt<-coxph(kmsurvo~.,data=combined[,11:ncol(combined)])
fitt

#9. Combined analysis for Clinical and rna Expression survival analysis

combined_reor<-combined[,c(1,7,10,2:6,8,9,11:ncol(combined))]
View(combined_reor)
combined_reor$race=factor(combined_reor$race)
fitt_grand<-coxph(kmsurvo~.-1,data=combined_reor[,4:ncol(combined_reor)])

fitt_grand
