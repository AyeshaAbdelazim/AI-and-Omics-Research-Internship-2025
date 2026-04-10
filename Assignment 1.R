dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")

data = read.csv(file.choose())
View(data)
str(data)

data$diagnosis = as.factor(data$diagnosis)
class(data$diagnosis)
levels(data$diagnosis)
data$diagnosis = factor(data$diagnosis , levels = c("Normal" , "Cancer"))

data$gender <- as.factor(data$gender)
data$gender_num <- ifelse(data$gender == "Female", 1, 0)

str(data)
data$gender_num <- as.factor(data$gender_num)
str(data)
data$smoker_binary = data$smoker
data$smoker_binary <- ifelse(data$smoker_binary == "Yes", 1, 0)


write.csv(data , file="clean_data/patient_info_clean.csv")
