
####1. Check Cholesterol level (using if)#### 

cholesterol <- 230
if(cholesterol > 240)
{print("High Cholesterol")
}


####2. Blood Pressure Status (using if...else)#### 

Systolic_bp <- 130

if(Systolic_bp < 120) {
  print ("Blood Pressure is normal")} else {
    print ("Blood Pressure is high")
  }

####3. Automating Data Type Conversion with for loop####

#patient_info
patient_info <- read.csv(file.choose())
patient.copy <- patient_info
str(patient.copy)
patient.factor_cols <- c("gender", "diagnosis", "smoker")

for (p in patient.factor_cols) { 
  patient.copy[[p]] <- as.factor(patient.copy[[p]])
  
}


#metadata
metadata <- read.csv(file.choose())
meta.copy <- metadata
str(meta.copy)
meta.factor_cols <- c("height", "gender")

for (m in meta.factor_cols) { 
  meta.copy[[m]] <- as.factor(meta.copy[[m]])
}


####4. Converting Factors to Numeric Codes####

#patient_info
binary.patient_cols <- c("smoker")

for (s in binary.patient_cols) {
  patient.copy[[s]] <- ifelse(patient.copy[[s]] == "Yes", 1, 0)
}


#metadata
binary.meta_cols <- c("gender")

for (g in binary.meta_cols) {
  meta.copy[[g]] <- ifelse(meta.copy[[g]] == "Female", 1, 0)
}

####Verification:####
str(patient_info)
str(patient.copy)

str(metadata)
str(meta.copy)

write.csv(patient.copy, "patient_copy.csv", row.names = FALSE)
write.csv(meta.copy, "meta.copy.csv", row.names = FALSE)



