library(arules)
library(arulesSequences)
library(ggplot2)
library(dplyr)

# Loading diab_trans.data file
diabtrans_df=read.csv("diab_trans.txt",header=T,stringsAsFactors=F)
diabtrans_df$Code=as.numeric(substring(diabtrans_df$Code,4,5))
head(diabtrans_df)

# Loading diab_description.data file 
diabdes_df=read.csv("diab_description.txt",header=T)
head(diabdes_df)

# Inner join operation
df_input=inner_join(diabtrans_df,diabdes_df)
head(df_input)

# Grouping the data by Code and selecting only distinct value(i.e attribute). Filter is applied on the rows 
# which have count greater than 1
to_chunkdf=df_input %>%
  group_by(Code) %>% distinct(value) %>% summarise(count=n()) %>% filter(count>1)
head(to_chunkdf)

# Adding new column called measure 
to_chunkdf$measure=T
head(to_chunkdf)

# left join 
base_df=left_join(df_input,to_chunkdf)
head(base_df)

# Created a variable called type which holds values based on certain condition with reference of code variable. 
dose=c(33,34,35)
base_df=base_df %>% mutate(type=ifelse(measure==T,ifelse(Code %in% dose,'dose','measurement'),'event'))
head(base_df)

#to_chunk_temp.df=inner_join(diabdes_df,to_chunkdf)
#head(to_chunk_temp.df)

# Histogram plotting of different types(i.e dose,measurement,event)
#to_hist_measurement.df=base_df %>% na.omit() %>% filter(measure=='T',type=='measurement')
#ggplot(to_hist_measurement.df,aes(x=value,y=Code))+geom_histogram()

#to_hist_dose.df=base_df %>% na.omit() %>% filter(measure=='T',type=='dose')
#ggplot(to_hist_dose.df, aes(x=value,y=Code))+geom_histogram(binwidth = 1)+scale_x_continuous(limits = c(2,40))

to_hist.df <- base_df %>% na.omit() %>% filter(measure==TRUE, Code==62)
head(to_hist.df)

# Histogram plot of values for Code=62
q <- quantile(to_hist.df$value, c(0.20, 0.80))
ggplot(to_hist.df, aes(x=value))+stat_bin(bins=30)+geom_histogram()+geom_vline(xintercept = c(q[1], q[2]))

# Omitting null values and selecting rows based on measure and extract distinct codes
to_divide.df=base_df %>% na.omit() %>% filter(measure=TRUE)
distinct_df=to_divide.df %>% distinct(Code)

qa=c()
qb=c()
code=c()
for(i in distinct_df$Code){
  to_temp_divide.df=base_df %>%na.omit() %>% filter(measure==T,Code==i)
  qa=c(qa,quantile(to_temp_divide.df$value,c(0.20)))
  qb=c(qb,quantile(to_temp_divide.df$value,c(0.80)))
  code=c(code,i)
}
print(qa)
print(qb)
print(code)

df=data.frame(qa,qb,code)

cluster=function(code_id,value){
  if(is.na(value))
    return(NA)
  if(is.null(value))
    return(NA)
  if(value<qa){
    return('low')
  }
  else if(value>qb){
    return('high')
  }
  else if(value>=qa & value<=qb){
    return('normal')
  }
  else{
    return(NA)}
}
base_df$value.level=mapply(cluster,base_df$Code,base_df$value)
head(base_df)

df <- base_df %>% mutate(description=if_else(!is.na(value.level),paste(Description,"level",value.level),as.character(Description))) %>%
  select(patient_id,time_sek,description)
head(df)

new_id=df %>% arrange(description) %>% distinct(description)
new_id$id=seq.int(nrow(new_id))

df = inner_join(df, new_id)
head(df)

write.table(df[,c(1,2,3)],file="clean_data.csv",row.names = F,col.names = F,sep=",")

dseq <- read_baskets(con = "clean_data.csv", sep =",", info = c("sequenceID","eventID"))
summary(dseq)

frameS =   as(dseq,"data.frame")
View(frameS)

nitems(dseq)

# getting frequency of items
freqItem = itemFrequency(dseq)
freqItem = sort(freqItem, decreasing = TRUE )
freqItem
head(freqItem,30)

# C-Spade algorithm parameters 
# max-size: 1
# max-gap: 28,800
# support=0.05
# max_len: 4
parameter = new ("SPparameter",support = 0.05, maxsize = 1, maxgap = 28800, maxlen=4)
CSpade= cspade(dseq, parameter, control = list(verbose = TRUE, tidLists = FALSE, summary= TRUE))

regularly=ruleInduction(CSpade,confidence = 0.8)
length(regularly)
inspect(head(regularly,30))

hypo=subset(regularly,!(lhs(regularly) %in% c("\"Hypoglycemic symptoms level low\"")) & rhs(regularly) %in% c("\"Hypoglycemic symptoms level low\""))
hypo = sort(hypo, by = "lift", decreasing=TRUE)
inspect(hypo)
