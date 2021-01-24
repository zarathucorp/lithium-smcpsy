library(data.table);library(magrittr);library(DT);library(jstable);library(dplyr)

#setwd("~/ShinyApps/jihyunbaek/lithium")
lithium <- readRDS("lithium.RDS")


# CKD-EPI----------------------------------------
CKDEPI<-function(scr,age,sex){
  if(sex=="F"){ k<-0.7; alpha<-(-0.329); const<-1.018 }
  else{ k<-0.9; alpha<-(-0.411); const<-1 }
  GFR <- 141 * min(scr/k,1)^alpha * max(scr/k,1)^(-1.209) * 0.993^age * const
  return(GFR)
}

df <- lithium$MEDI[, c("NO","처방일","처방명","일수","횟수")] 
names(df) <- c("NO","date","drug","day","times")
df[, drug := ifelse(drug == "Lithium carbonate 300mg", "Lithium", "Valproate")]
df <- unique(df)
df <- df[order(-day), .(maxday = max(day, na.rm = T), maxnotqd = day[which(times != 1)[1]]), by=c("NO","date","drug")]
df <- df[, .(maxday, qd = ifelse(is.na(maxnotqd), maxday , maxday - maxnotqd)), by=c("NO","date","drug")]
df <- df[, .(totDay = sum(maxday, na.rm = T), qd = sum(qd, na.rm = T)/sum(maxday, na.rm = T)),by = c("NO", "drug")]

df.long <- dcast(df, NO ~ drug, value.var = c("totDay", "qd"))

#left_join 위해서 NO의 class를 맞춰주기----------------------------------------
lithium$`clinical data`$NO <- lithium$`clinical data`$NO %>% as.numeric() %>% as.character()
df.long$NO <- df.long$NO %>% as.character()
lithium$`clinical data` <- merge(lithium$`clinical data`, df.long, by = "NO")
## Dx group
lithium$`clinical data`[, group_bipolar_schizoaffective_other := factor(ifelse(grepl("Bipolar|bipolar", lithium$`clinical data`$주상병명), "Bipolar disorder",
                                                                               ifelse(grepl("Schizoaffective|schizoaffective", lithium$`clinical data`$주상병명), "Schizoaffective disorder", "vOthers")))]


# Data inclusion - 총처방일수 180일 초과 : tot = 4360----------------------------------------
a <- lithium$`clinical data`[xor(is.na(totDay_Lithium),is.na(totDay_Valproate)) & (totDay_Lithium>180 | totDay_Valproate>180),
                             .(NO,성별,생년월일,totDay_Lithium,totDay_Valproate,qd_Lithium,qd_Valproate, 
                               HTN = factor(as.integer(!is.na(`고혈압 여부`))), DM = factor(as.integer(!is.na(`당뇨 여부`))), group_bipolar_schizoaffective_other)]


N_profile<-cbind("총처방일수 180일 초과",NA,nrow(a))

## Date age----------------------------------------
df <- lithium$MEDI[, NO := as.character(NO)][,.SD,]
setnames(df,c("처방일","함량단위투여량","일수"),c("date","dose","day"))

data.main <- a %>% 
  merge(df[,.(firstPrescriptionDay=min(date, na.rm = T)), by = "NO"], by = "NO",all.x = T) %>% 
  merge(df[,.(lastPrescriptionDay=max(date, na.rm = T)), by = "NO"], by = "NO",all.x = T) %>% 
  merge(df[, .(avgDose_1day = sum(dose * day)/sum(day)), by = "NO"], by = "NO",all.x = T)

data.main[, Age := floor(as.numeric(as.Date(firstPrescriptionDay) - as.Date(`생년월일`))/365.25)]

data.main <- data.main[Age>=18,,]
N_profile<-rbind(N_profile,cbind("첫처방일기준 만 18세 이상",NA,data.main[,.N,]))

# LithiumToxicity----------------------------------------

df <- lithium$`renal function & TDM`[, NO := as.character(NO)][] %>% 
  merge(data.main[, .(NO, firstPrescriptionDay, lastPrescriptionDay)], by = "NO", all.x = T)
setnames(df,c("세부검사명","결과","시행일시"),c("test","result","testdate"))

data.main <- data.main %>% 
  merge(df[test=="Lithium" & as.numeric(result) > 1.0 & (testdate - firstPrescriptionDay >= 0) & (lastPrescriptionDay - testdate  >= 0),  .(LithiumToxicity1.0 = .N), by="NO"], by="NO", all.x = T) %>% 
  merge(df[test=="Lithium" & as.numeric(result) > 0.8 & (testdate - firstPrescriptionDay >= 0) & (lastPrescriptionDay - testdate  >= 0), .(LithiumToxicity0.8 = .N), by="NO"], by="NO", all.x = T) %>% 
  merge(df[test=="Lithium" & as.numeric(result) > 1.2 & (testdate - firstPrescriptionDay >= 0) & (lastPrescriptionDay - testdate  >= 0), .(LithiumToxicity1.2 = .N), by="NO"], by="NO", all.x = T) %>% 
  merge(df[test=="Lithium" & (testdate - firstPrescriptionDay >= 0), .(avgTDM_Lithium = mean(as.numeric(result), na.rm = T)), by="NO"], by="NO", all.x = T) %>% 
  merge(df[test=="Valproic Acid" & (testdate - firstPrescriptionDay >= 0), .(avgTDM_Valproate = mean(as.numeric(result), na.rm = T)), by="NO"], by="NO", all.x = T)

for (v in c("LithiumToxicity1.0", "LithiumToxicity0.8", "LithiumToxicity1.2")){
  data.main[[v]] <- ifelse(is.na(data.main[[v]]), 0, data.main[[v]])
}


df<-lithium$`renal function & TDM`
df$NO <- as.character(df$NO)
df <- merge(df, lithium$`clinical data`[,.(NO,`성별`,`생년월일`),], by="NO", mult=all)
df$`결과`<- as.numeric(df$`결과`) 
df$시행일시 <- as.Date(df$시행일시); df$생년월일 <- as.Date(df$생년월일)
setnames(df,c("세부검사명","시행일시","생년월일","결과","성별"),c("test", "testDate", "birthDate", "result", "sex"))
df[,age:=as.numeric(testDate-birthDate)/365.25]
df[,eGFR:=ifelse(test=="Creatinine",CKDEPI(result,age,sex),NA),by=seq_len(nrow(df))]


## data for figure 1----------------------------------------
data.f1 <- df[!is.na(eGFR),.(NO,testDate,eGFR)]
setnames(data.f1, "testDate", "date")


## Main data----------------------------------------

data.main <- merge(data.main, df[eGFR < 60, .(eGFRbelow60Date = min(testDate)), by = "NO"], all.x = TRUE) %>% 
  merge(df[test == "Creatinine", .(testNum = .N), by="NO"], by="NO", all.x=TRUE) %>% 
  .[!is.na(testNum) & (is.na(eGFRbelow60Date) | as.Date(firstPrescriptionDay) < as.Date(eGFRbelow60Date))]

data.main<-merge(data.main,data.f1[,.(lastTestDate=max(date)),by="NO"])
data.main<-data.main[testNum>=2 & (as.Date(lastTestDate)-as.Date(firstPrescriptionDay))/365.25>=0.5,,]
N_profile<-rbind(N_profile,cbind("최소 2개 이상의 eGFR data\n(baseline & 최소 6개월 이상의 post-baseline data)",NA,data.main[,.N,]))

data.main[, duration := ifelse(is.na(eGFRbelow60Date),as.Date(lastPrescriptionDay) - as.Date(firstPrescriptionDay), as.Date(eGFRbelow60Date) - as.Date(firstPrescriptionDay))]
# duration Full
data.main[, year_FU_full := as.numeric(as.Date(lastPrescriptionDay) - as.Date(firstPrescriptionDay))/365.25]
data.main[, drug := factor(ifelse(is.na(totDay_Lithium), 0, 1))]
data.main[, eGFRbelow60 := factor(as.integer(!is.na(eGFRbelow60Date)))]
data.main[, `:=`(year_FU= duration/365.25, totYear_Lithium = totDay_Lithium/365.25, totYear_Valproate = totDay_Valproate/365.25)]
setnames(data.main, "성별", "Sex")
data.main[, Sex := factor(Sex)]

data.main <- data.main[, .SD, .SDcols = -c("생년월일", "firstPrescriptionDay", "lastPrescriptionDay", "eGFRbelow60Date", "duration", "totDay_Valproate", "totDay_Lithium","testNum")]

## Figure 1 data----------------------------------------

# 처방 정보

df <- lithium$MEDI[, c("NO","처방일","처방명","일수","횟수")]
names(df) <- c("NO","date","drug","day","times")
df[, drug := factor(ifelse(drug == "Lithium carbonate 300mg", 1, 0))]
df <- unique(df)[, `:=`(NO = as.character(NO), date = as.Date(date))][]
df <- df[, .(maxday = max(day, na.rm = T)), by=c("NO","date","drug")]

# 일단 합친 후 cumsum
data.f1 <- rbindlist(list(data.f1, df),use.names = TRUE, fill=TRUE)[order(NO,date)][, maxday:=ifelse(is.na(maxday),0,maxday)][]
data.f1[, cumulativePrescriptionDay := cumsum(maxday),by=.(NO)]

data.f1 <- data.f1[!is.na(eGFR), !c("maxday","drug")]
data.f1 <- merge(data.f1, data.main[,.(NO,drug),], by="NO")

## Base eGFR
#data.f1[, .(base_eGFR = ifelse(any(as.integer(cumulativePrescriptionDay)) == 0, eGFR[max(which(cumulativePrescriptionDay == 0))], eGFR[1])), by ="NO"]

## GFR change----------------------------------------

data.GFRchange <- data.f1[cumulativePrescriptionDay<365.25,.(year0GFR=mean(eGFR,na.rm=T)),by="NO"] %>% 
  merge(., data.f1[, .(base_eGFR = ifelse(any(as.integer(cumulativePrescriptionDay)) == 0, eGFR[max(which(cumulativePrescriptionDay == 0))], eGFR[1])), by ="NO"], by= "NO", all = T) %>% ## base eGFR
  merge(.,data.f1[(365.25*3)<cumulativePrescriptionDay & cumulativePrescriptionDay<(365.25*4),.(year3GFR=mean(eGFR,na.rm=T)),by="NO"],all=T) %>% 
  merge(.,data.f1[(365.25*5)<cumulativePrescriptionDay & cumulativePrescriptionDay<(365.25*6),.(year5GFR=mean(eGFR,na.rm=T)),by="NO"],all=T) %>% 
  merge(.,data.f1[(365.25*7)<cumulativePrescriptionDay & cumulativePrescriptionDay<(365.25*8),.(year7GFR=mean(eGFR,na.rm=T)),by="NO"],all=T) %>% 
  merge(.,data.f1[(365.25*10)<cumulativePrescriptionDay & cumulativePrescriptionDay<(365.25*11),.(year10GFR=mean(eGFR,na.rm=T)),by="NO"],all=T) %>% 
  merge(.,data.f1[(365.25*12)<cumulativePrescriptionDay & cumulativePrescriptionDay<(365.25*13),.(year12GFR=mean(eGFR,na.rm=T)),by="NO"],all=T) %>% 
  merge(.,data.f1[(365.25*15)<cumulativePrescriptionDay & cumulativePrescriptionDay<(365.25*16),.(year15GFR=mean(eGFR,na.rm=T)),by="NO"],all=T) %>% 
  merge(.,data.f1[(365.25*20)<cumulativePrescriptionDay & cumulativePrescriptionDay<(365.25*21),.(year20GFR=mean(eGFR,na.rm=T)),by="NO"],all=T)

data.main<-merge(data.main,data.GFRchange,all=T)

## ----------------------------------------

ICD_data <- readRDS("ICD_data.RDS")
setnames(ICD_data,c("개인정보동의여부","정렬순서"),c("Privacy Consent","NO"))
ICD_data$NO<-ICD_data$NO %>% as.character()

ICD_data<-merge(ICD_data,data.main[,.(NO,base_eGFR),],by="NO")

ICD_data<-ICD_data[`Privacy Consent`=="Y" & !is.na(base_eGFR),,]
N_profile<-rbind(N_profile,cbind("개인정보사용미동의",as.integer(N_profile[nrow(N_profile),3])-ICD_data[,.N,],ICD_data[,.N,]))
                 
ICD_data<-ICD_data[base_eGFR>=30,,]
N_profile<-rbind(N_profile,cbind("baseline eGFR<30",as.integer(N_profile[nrow(N_profile),3])-ICD_data[,.N,],ICD_data[,.N,]))


ICD_data<-ICD_data[,alldiagnosis:=Reduce(paste,.SD),.SDcols=grep("진단코드",colnames(ICD_data))][,c("NO","base_eGFR","alldiagnosis"),]
ICD_data<-ICD_data[!(alldiagnosis %like% "N0|N1" & !(alldiagnosis %like% "N09")),.SD,]
N_profile<-rbind(N_profile,cbind("ICD N00-N08 or N10-N19",as.integer(N_profile[nrow(N_profile),3])-ICD_data[,.N,],ICD_data[,.N,]))

ICD_data<-ICD_data[!(alldiagnosis %like% "T86.1"),.SD,]
N_profile<-rbind(N_profile,cbind("ICD T86.1",as.integer(N_profile[nrow(N_profile),3])-ICD_data[,.N,],ICD_data[,.N,]))

ICD_data<-ICD_data[!(alldiagnosis %like% "Z94.0"),.SD,]
N_profile<-rbind(N_profile,cbind("ICD Z94.0",as.integer(N_profile[nrow(N_profile),3])-ICD_data[,.N,],ICD_data[,.N,]))
colnames(N_profile)<-c("조건","제외","N")

data.main <- merge(data.main,ICD_data[,.(NO),],by="NO")


## 복용년수별 n수 ----------------------------------------

Year_N<-data.frame(Year=0:26,
           Lithium_N=sapply(0:26,function(x) data.main[totYear_Lithium>x,.N,]),
           Valproate_N=sapply(0:26,function(x) data.main[totYear_Valproate>x,.N,]))


## 5년 뒤, 10년 뒤 eGFR<60의 비율 (n수) ----------------------------------------

eGFRbelow60ratio<-rbind(c(5,data.main[!is.na(year5GFR),.N,],data.main[year5GFR<60,.N,],round(data.main[year5GFR<60,.N,]/data.main[!is.na(year5GFR),.N,]*100,3)),
                        c(10,data.main[!is.na(year10GFR),.N,],data.main[year10GFR<60,.N,],round(data.main[year10GFR<60,.N,]/data.main[!is.na(year10GFR),.N,]*100,3)))

colnames(eGFRbelow60ratio)<-c("Year","eGFR 측정값 있는 N","eGFR<60","Percentage")

## ----------------------------------------

data.main <- data.main[, -c("NO")]  ## NO 제외

label.main <- jstable::mk.lev(data.main)

label.main[variable == "eGFRbelow60", `:=`(var_label = "eGFR < 60", val_label = c("No", "Yes"))]
label.main[variable == "drug", `:=`(var_label = "Drug", val_label = c("Valproate", "Lithium"))]
label.main[variable == "DM", `:=`(var_label = "DM", val_label = c("No", "Yes"))]
label.main[variable == "HTN", `:=`(var_label = "HTN", val_label = c("No", "Yes"))]
label.main[variable == "LithiumToxicity1.0", `:=`(var_label = "Lithium > 1.0 횟수")]
label.main[variable == "LithiumToxicity1.2", `:=`(var_label = "Lithium > 1.2 횟수")]
label.main[variable == "LithiumToxicity0.8", `:=`(var_label = "Lithium > 0.8 횟수")]
label.main[variable == "avgDose_1day", `:=`(var_label = "Average 1day dose")]
label.main[variable == "totYear_Lithium", `:=`(var_label = "Cumulative Lithium year")]
label.main[variable == "totYear_Valproate", `:=`(var_label = "Cumulative Valproate year")]
label.main[variable == "qd_Lithium", `:=`(var_label = "Lithium QD proportion")]
label.main[variable == "qd_Valproate", `:=`(var_label = "Valproate QD proportion")]

label.main[variable == "year0GFR", `:=`(var_label = "복용 1년 이내 GFR")]
label.main[variable == "year3GFR", `:=`(var_label = "복용 3년차 GFR")]
label.main[variable == "year5GFR", `:=`(var_label = "복용 5년차 GFR")]
label.main[variable == "year7GFR", `:=`(var_label = "복용 7년차 GFR")]
label.main[variable == "year10GFR", `:=`(var_label = "복용 10년차 GFR")]
label.main[variable == "year12GFR", `:=`(var_label = "복용 12년차 GFR")]
label.main[variable == "year15GFR", `:=`(var_label = "복용 15년차 GFR")]
label.main[variable == "year20GFR", `:=`(var_label = "복용 20년차 GFR")]

## variable order: 미리 만들어놓은 KM, cox 모듈용
varlist_kmcox <- list(variable = c("eGFRbelow60", "year_FU", "drug", setdiff(names(data.main), c("eGFRbelow60", "year_FU", "drug" ))))