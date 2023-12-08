# ==============================================================================
# ================================ Description =================================
# ==============================================================================
# Author: Amirreza Sahebi
# Last modified: 12.01.2023
# Input: Comprehensive data set detailing ICU and ED encounters, code status, 
#        and ventilator use across several hospitals.
# Output: Statistical analysis results, visualizations, and summaries of 
#         code status transitions and trends, particularly related to the impact
#         of COVID-19 on patient care choices.
# Purpose: To analyze and understand the trends in in-hospital code status 
#          updates over time, with a focus on the impact of the COVID-19 
#          on patient code status choices in a Level I Trauma hospital setting.


# ============================================================================ #
# ============================= Required libraries =========================== #
# ============================================================================ #
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(car)
library(lubridate)
library(emmeans)
library(semTools)
library(gridExtra)
library(scales)
library(grid)
library(MASS)
library(boot)


# ============================================================================ #
# ================================== Functions =============================== #
# ============================================================================ #
# Function factory for secondary axis transforms
train_sec <- function(primary, secondary, na.rm = TRUE) {
  # Thanks Henry Holm for including the na.rm argument!
  from <- range(secondary, na.rm = na.rm)
  to   <- range(primary, na.rm = na.rm)
  # Forward transform for the data
  forward <- function(x) {
    rescale(x, from = from, to = to)
  }
  # Reverse transform for the secondary axis
  reverse <- function(x) {
    rescale(x, from = to, to = from)
  }
  list(fwd = forward, rev = reverse)
}


# ============================================================================ #
# ================================== Read Data =============================== #
# ============================================================================ #
# data path -- please modify paths accordingly
pathread <- "G:/My Drive/North Carolina State University/Project - Code Status 0/data/" 
pathwrite <- "G:/My Drive/North Carolina State University/Project - Code Status 0/result/"
# Comprehensive dataset 'care-planning-capacity-data-v5.csv' detailing ICU and ED encounters,
# Code status, and ventilator use across seven hospitals from 3/1/2019 to 12/31/2022. 
# We are not allowed to share this dataset.
df <- fread(paste0(pathread,"care-planning-capacity-data-v5.csv"))


# ============================================================================ #
# ================================ Data Cleaning ============================= #
# ============================================================================ #
df <- as.data.frame(df)
names(df)[names(df) %in% c("date","hospital_code")] <- c("Date","Hospital_Code")
df$Date <- as.Date(df$Date)
# Second Hospital data has problems and removed
df <- subset(df, Hospital_Code != 2)
df$Hospital_Code <- ifelse(df$Hospital_Code == 1,df$Hospital_Code, df$Hospital_Code - 1)

# Filters out rows in df with discrepancies between old and new census.
OldEDCols <- names(df)[grepl("ED_CENSUS", names(df)) & !grepl("CURRENT", names(df))]
NewEDCols <- names(df)[grepl("ED_CENSUS", names(df)) & grepl("CURRENT", names(df))]
OldICUCols <- names(df)[grepl("ICU_CENSUS", names(df)) & !grepl("CURRENT", names(df))]
NewICUCols <- names(df)[grepl("ICU_CENSUS", names(df)) & grepl("CURRENT", names(df))]
df$diffED <- rowSums(df[,OldEDCols]) - rowSums(df[,NewEDCols])
df$diffICU <- rowSums(df[,OldICUCols]) - rowSums(df[,NewICUCols])
# 14 records after Dec 18, 2022 has 1 extra record in ED old census
df <- subset(df, diffICU == 0 & diffED == 0)
df <- subset(df, select = -c(diffICU, diffED))


# ============================================================================ #
# ============================= Data Pre-Processing ========================== #
# ============================================================================ #
# Defining COVID time
COVIDend <- as.Date("2022-05-01")     # end of o-micron wave
COVIDstart <- as.Date("2020-04-01")     # first case in the state of the study


# Separate old and new census data
dfOld <- df[,c("Date","Hospital_Code","N_vent_pts","ventilator_pt_census",OldEDCols,OldICUCols)]
colnames(dfOld) <- c("Date","Hospital_Code","Vent1","Vent2",
                       "ED_Lim","ED_DNR","ED_Full","ED_Non",
                       "ICU_Lim","ICU_DNR","ICU_Full","ICU_Non")

dfNew <- df[,c("Date","Hospital_Code","N_vent_pts","ventilator_pt_census",NewEDCols,NewICUCols)]
colnames(dfOld) <- c("Date","Hospital_Code","Vent1","Vent2",
                       "ED_Lim","ED_DNR","ED_Full","ED_Non",
                       "ICU_Lim","ICU_DNR","ICU_Full","ICU_Non")


# Separate code status columns
dfSub <- df[,c("Date","Hospital_Code",OldEDCols,OldICUCols,NewEDCols,NewICUCols)]
colnames(dfSub) <- c("Date","Hospital_Code",
                     "ED_LC_Old","ED_DNR_Old","ED_FC_Old","ED_NS_Old",
                     "ICU_LC_Old","ICU_DNR_Old","ICU_FC_Old","ICU_NS_Old",
                     "ED_LC_New","ED_DNR_New","ED_FC_New","ED_NS_New",
                     "ICU_LC_New","ICU_DNR_New","ICU_FC_New","ICU_NS_New")
dfSub[,"ED_Sum"] <- rowSums(dfSub[, c("ED_LC_Old","ED_DNR_Old","ED_FC_Old","ED_NS_Old")])
dfSub[,"ICU_Sum"] <- rowSums(dfSub[, c("ICU_LC_Old","ICU_DNR_Old","ICU_FC_Old","ICU_NS_Old")])
dfSub[,"ICU_LC_Diff"] <- dfSub[,"ICU_LC_New"] - dfSub[,"ICU_LC_Old"]
dfSub[,"ICU_FC_Diff"] <- dfSub[,"ICU_FC_New"] - dfSub[,"ICU_FC_Old"]
dfSub[,"ICU_DNR_Diff"] <- dfSub[,"ICU_DNR_New"] - dfSub[,"ICU_DNR_Old"]
dfSub[,"ICU_NS_Diff"] <- dfSub[,"ICU_NS_New"] - dfSub[,"ICU_NS_Old"]
dfSub[,"ED_LC_Diff"] <- dfSub[,"ED_LC_New"] - dfSub[,"ED_LC_Old"]
dfSub[,"ED_FC_Diff"] <- dfSub[,"ED_FC_New"] - dfSub[,"ED_FC_Old"]
dfSub[,"ED_DNR_Diff"] <- dfSub[,"ED_DNR_New"] - dfSub[,"ED_DNR_Old"]
dfSub[,"ED_NS_Diff"] <- dfSub[,"ED_NS_New"] - dfSub[,"ED_NS_Old"]
dfSub[,"ICU_LC_ChngPct"] <- dfSub[,"ICU_LC_Diff"] / dfSub[,"ICU_Sum"] * 100
dfSub[,"ICU_FC_ChngPct"] <- dfSub[,"ICU_FC_Diff"] / dfSub[,"ICU_Sum"] * 100
dfSub[,"ICU_DNR_ChngPct"] <- dfSub[,"ICU_DNR_Diff"] / dfSub[,"ICU_Sum"] * 100
dfSub[,"ICU_NS_ChngPct"] <- dfSub[,"ICU_NS_Diff"] / dfSub[,"ICU_Sum"] * 100
dfSub[,"ED_LC_ChngPct"] <- dfSub[,"ED_LC_Diff"] / dfSub[,"ED_Sum"] * 100
dfSub[,"ED_FC_ChngPct"] <- dfSub[,"ED_FC_Diff"] / dfSub[,"ED_Sum"] * 100
dfSub[,"ED_DNR_ChngPct"] <- dfSub[,"ED_DNR_Diff"] / dfSub[,"ED_Sum"] * 100
dfSub[,"ED_NS_ChngPct"] <- dfSub[,"ED_NS_Diff"] / dfSub[,"ED_Sum"] * 100


# Convert wide to long
dfSubLong <- dfSub[,c(1:20)] %>% 
  dplyr::select(-c("Hospital_Code")) %>%
  group_by(Date) %>%
  summarize_at(vars(matches("^(ED|ICU)_.*")), sum) %>%
  group_by(Date) %>%
  mutate(across(matches("ED|ICU"), ~ ifelse(cur_column() %>% str_detect("ICU"), round(. / ICU_Sum * 100,2), round(. / ED_Sum * 100, 2)))) %>%
  dplyr::select(-c("ICU_Sum","ED_Sum")) %>%
  pivot_longer(cols = matches("Old|New"), names_to = c("Section", "CodeType", "Census"), names_sep = "_", values_to = "Value")
dfSubLong$Census <- factor(dfSubLong$Census, levels = c("Old", "New"))


# Hospital 6 is the largest one and analyzed separately
dfSub6Long <- dfSub[dfSub$Hospital_Code == 6,c(1:20)] %>% 
  dplyr::select(-c("Hospital_Code")) %>%
  group_by(Date) %>%
  summarize_at(vars(matches("^(ED|ICU)_.*")), sum) %>%
  group_by(Date) %>%
  mutate(across(matches("ED|ICU"), ~ ifelse(cur_column() %>% str_detect("ICU"), round(. / ICU_Sum * 100,2), round(. / ED_Sum * 100, 2)))) %>%
  dplyr::select(-c("ICU_Sum","ED_Sum")) %>%
  pivot_longer(cols = matches("Old|New"), names_to = c("Section", "CodeType", "Census"), names_sep = "_", values_to = "Value")
dfSub6Long$Census <- factor(dfSub6Long$Census, levels = c("Old", "New"))


# Looking at code change rate after updating the code status documentation system
ColList <- c("Date","Hospital_Code",names(dfSub)[grep("_ChngPct",names(dfSub))])
dfCodeTransition <- dfSub[,ColList]
dfCodeTransition <- gather(dfCodeTransition, condition, measurement, 3:10, factor_key=TRUE)
dfCodeTransition$Section <- word(dfCodeTransition$condition, 1, sep = "_")
dfCodeTransition$Code <- word(dfCodeTransition$condition, 2, sep = "_")
names(dfCodeTransition)[which(names(dfCodeTransition) == "measurement")] = "ChangePct"
names(dfCodeTransition)[which(names(dfCodeTransition) == "Hospital_Code")] = "HospitalCode"
dfCodeTransition = subset(dfCodeTransition, select = -c(condition))
dfCodeTransition$quarter <- paste(year(dfCodeTransition$Date), "-Q",
                                  quarter(dfCodeTransition$Date), sep = "")
dfCodeTransition$CovidStatus <- "in-Covid"
dfCodeTransition$CovidStatus[which(dfCodeTransition$Date < COVIDstart)] <- "pre-Covid"
dfCodeTransition$CovidStatus[which(dfCodeTransition$Date > COVIDend)] <- "post-Covid" # after Omicron
dfCodeTransition$CovidStatus <- factor(dfCodeTransition$CovidStatus, levels = c("pre-Covid", "in-Covid", "post-Covid"))
dfCodeTransition$SectionCode <- factor(paste0(dfCodeTransition$Section, "-", dfCodeTransition$Code),
                                       levels = c("ICU-FC","ED-FC","ICU-DNR","ED-DNR","ICU-LC","ED-LC","ICU-NS","ED-NS"))
dfCodeTransition$HospitalCode <- factor(paste0("Hospital ",dfCodeTransition$HospitalCode),
                                        levels = c("Hospital 1","Hospital 2","Hospital 3","Hospital 4","Hospital 5","Hospital 6"))


# Separating only NS code transitions
dfNSTransitionICU <- subset(dfSub, abs(ICU_NS_Diff) == (ICU_FC_Diff + ICU_LC_Diff + ICU_DNR_Diff) 
                              & ICU_NS_Diff != 0 & ICU_NS_Old > 0 & ICU_FC_Diff >= 0 & ICU_LC_Diff >= 0 & ICU_DNR_Diff >= 0)
dfNSTransitionED <- subset(dfSub, abs(ED_NS_Diff) == (ED_FC_Diff + ED_LC_Diff + ED_DNR_Diff) &
                               ED_NS_Diff != 0 & ED_NS_Old > 0 & ED_FC_Diff >= 0 & ED_LC_Diff >= 0 & ED_DNR_Diff >= 0)
dfNSTransitionICU <- dfNSTransitionICU %>%  
  mutate(ICU_NS2LCPct = round(ICU_LC_Diff / ICU_NS_Old * 100,2),
         ICU_NS2DNRPct = round(ICU_DNR_Diff / ICU_NS_Old * 100,2),
         ICU_NS2FCPct = round(ICU_FC_Diff / ICU_NS_Old * 100,2),
         ICU_NS2NSPct = round(ICU_NS_New / ICU_NS_Old * 100,2))
dfNSTransitionICU$CovidStatus <- "in-Covid"
dfNSTransitionICU$CovidStatus[which(dfNSTransitionICU$Date < COVIDstart)] <- "pre-Covid"
dfNSTransitionICU$CovidStatus[which(dfNSTransitionICU$Date > COVIDend)] <- "post-Covid" 
dfNSTransitionICU$CovidStatus <- factor(dfNSTransitionICU$CovidStatus,levels = c("pre-Covid","in-Covid","post-Covid"))
dfNSTransitionED <- dfNSTransitionED %>%  
  mutate(ED_NS2LCPct = round(ED_LC_Diff / ED_NS_Old * 100,2),
         ED_NS2DNRPct = round(ED_DNR_Diff / ED_NS_Old * 100,2),
         ED_NS2FCPct = round(ED_FC_Diff / ED_NS_Old * 100,2),
         ED_NS2NSPct = round(ED_NS_New / ED_NS_Old * 100,2))
dfNSTransitionED$CovidStatus <- "in-Covid"
dfNSTransitionED$CovidStatus[which(dfNSTransitionED$Date < COVIDstart)] <- "pre-Covid"
dfNSTransitionED$CovidStatus[which(dfNSTransitionED$Date > COVIDend)] <- "post-Covid"
dfNSTransitionED$CovidStatus <- factor(dfNSTransitionED$CovidStatus,levels = c("pre-Covid","in-Covid","post-Covid"))


# Regression data. Only ventilator columns and new census data
ColList <- c("Date","Hospital_Code",names(df)[grep("vent",names(df))],names(df)[grep("CURRENT",names(df))])
dfReg <- df[,ColList]
colnames(dfReg) <- c("Date","Hospital_Code","Vent1","Vent2",
                     "ED_LC","ED_DNR","ED_FC","ED_NS",
                     "ICU_LC","ICU_DNR","ICU_FC","ICU_NS")
dfReg[,"ICU_DNR_Rt"] <- dfReg[,"ICU_DNR"] / rowSums(dfReg[, c("ICU_LC","ICU_DNR","ICU_FC","ICU_NS")])
dfReg[,"ICU_FC_Rt"] <- dfReg[,"ICU_FC"] / rowSums(dfReg[, c("ICU_LC","ICU_DNR","ICU_FC","ICU_NS")])
dfReg[,"ICU_LC_Rt"] <- dfReg[,"ICU_LC"] / rowSums(dfReg[, c("ICU_LC","ICU_DNR","ICU_FC","ICU_NS")])
dfReg[,"ICU_NS_Rt"] <- dfReg[,"ICU_NS"] / rowSums(dfReg[, c("ICU_LC","ICU_DNR","ICU_FC","ICU_NS")])
dfReg[,"ED_DNR_Rt"] <- dfReg[,"ED_DNR"] / rowSums(dfReg[, c("ED_LC","ED_DNR","ED_FC","ED_NS")])
dfReg[,"ED_FC_Rt"] <- dfReg[,"ED_FC"] / rowSums(dfReg[, c("ED_LC","ED_DNR","ED_FC","ED_NS")])
dfReg[,"ED_LC_Rt"] <- dfReg[,"ED_LC"] / rowSums(dfReg[, c("ED_LC","ED_DNR","ED_FC","ED_NS")])
dfReg[,"ED_NS_Rt"] <- dfReg[,"ED_NS"] / rowSums(dfReg[, c("ED_LC","ED_DNR","ED_FC","ED_NS")])
dfReg$CovidStatus <- "in-Covid"
dfReg$CovidStatus[which(dfReg$Date < COVIDstart)] <- "pre-Covid"
dfReg$CovidStatus[which(dfReg$Date > COVIDend)] <- "post-Covid"
dfReg$CovidStatus <- factor(dfReg$CovidStatus,levels = c("pre-Covid","in-Covid","post-Covid"))

# ============================================================================ #
# ================================= Analysis ================================= #
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# -------------------------- Census Code Comparison -------------------------- #
# ---------------------------------------------------------------------------- #


# ---------------------- #
# All Hospitals together #
# ---------------------- #
p <- ggplot(dfSubLong, aes(x = CodeType, y = Value, fill = Census)) +
  geom_boxplot() +
  labs(x='Code Status', y = "Percentage") + 
  facet_wrap(~Section, ncol = 2) +
  scale_fill_manual(values = c("red", "blue")) +
  labs(colour = "Census") + 
  ggtitle("Comparison of Code Status Percentage in Old and New Censuses") +
  theme(strip.text.x = element_text(size = 28, color = "black", face = "bold",family = "mono"),
        strip.text.y = element_text(size = 28, color = "black", face = "bold",family = "mono"),
        legend.title = element_text(family = "mono",size = 24),
        legend.text = element_text(family = "mono",size = 24),
        legend.position="right",
        axis.title.y  = element_text(family = "mono",size = 24, face = "bold"),
        axis.title.x  = element_text(family = "mono",size = 24, face = "bold"),
        axis.text.x   = element_text(family = "mono",size = 24, face = "bold"),
        axis.text.y   = element_text(family = "mono",size = 24, face = "bold"),
        plot.title = element_text(family = "mono",size = 26,face = "bold", hjust = 0.5))
ggsave(paste0(pathwrite,"CodeStatusBoxPlot_AllHosp.jpg"), plot = p, width = 14, height = 8)

# Conduct statistical tests for comparing performance of new and old censuses
Tresult <- data.frame()
for (section in c("ICU","ED")){
  for (code in c("FC","LC","DNR","NS")){
    result <- list()
    g1 <- subset(dfSubLong, Section == section & CodeType == code & Census == "Old")$Value
    g2 <- subset(dfSubLong, Section == section & CodeType == code & Census == "New")$Value
    vartest <- var.test(g1,g2)
    ttest <- t.test(g2, g1, var.equal = vartest$p.value > 0.05)
    if (ttest$p.value > 0.1){
      result["Ttest"] <- "No Significant Difference"
    } else if (ttest$p.value > 0.05){
      result["Ttest"] <- ifelse(code == "NS", "Decreased (.)", "Increased (.)")
    } else if (ttest$p.value > 0.01){
      result["Ttest"] <- ifelse(code == "NS", "Decreased (*)", "Increased (*)")
    } else if (ttest$p.value > 0.001){
      result["Ttest"] <- ifelse(code == "NS", "Decreased (**)", "Increased (**)")
    } else {
      result["Ttest"] <- ifelse(code == "NS", "Decreased (***)", "Increased (***)")
    }
    result["EqualVariance"] <- ifelse(vartest$p.value <= 0.05,"No","Yes")
    result$CodeType <- code
    result$Section = section
    Tresult <- rbind(Tresult,result)
  }
}
tmp <- dfSubLong %>%
  group_by(Section,CodeType,Census) %>%
  summarize(mean = round(mean(Value, na.rm = T),2),
            lower_CI = round(quantile(Value, 0.025, na.rm = T),2),
            upper_CI = round(quantile(Value, 0.975, na.rm = T),2)) %>%
  mutate(CodePct = paste0(mean," [",lower_CI, ", ",upper_CI,"]")) %>%
  ungroup() %>%
  select(Section,CodeType,Census,CodePct) %>%
  spread(Census , CodePct) %>%
  merge(Tresult, by = c("CodeType","Section")) %>%
  arrange(factor(Section, levels = c("ICU", "ED")), match(CodeType, c("FC", "DNR", "LC", "NS")))
fwrite(tmp, paste0(pathwrite,"CodeStat.csv"))
knitr::kable(tmp)
#   |CodeType |Section |Old                  |New                 |Ttest           |EqualVariance |
#   |:--------|:-------|:--------------------|:-------------------|:---------------|:-------------|
#   |FC       |ICU     |61.81 [53.89, 70.29] |87.92 [82.5, 92.66] |Increased (***) |No            |
#   |DNR      |ICU     |2.97 [0.53, 6.63]    |7.86 [3.87, 12.44]  |Increased (***) |No            |
#   |LC       |ICU     |1.59 [0, 3.87]       |3.69 [1.09, 6.97]   |Increased (***) |No            |
#   |NS       |ICU     |33.63 [24.73, 42.15] |0.53 [0, 1.79]      |Decreased (***) |No            |
#   |FC       |ED      |32.25 [25.5, 39.72]  |47.1 [38.07, 55.27] |Increased (***) |No            |
#   |DNR      |ED      |2.89 [0.9, 5.15]     |3.85 [1.63, 6.33]   |Increased (***) |No            |
#   |LC       |ED      |0.28 [0, 1.1]        |0.34 [0, 1.21]      |Increased (***) |No            |
#   |NS       |ED      |64.58 [56.64, 71.84] |48.7 [40.56, 58.23] |Decreased (***) |No            |


# --------------- #
# Only Hospital 6 #
# --------------- #
p <- ggplot(dfSub6Long, aes(x = CodeType, y = Value, fill = Census)) +
  geom_boxplot() +
  labs(x='Code Status', y = "Percentage") + 
  facet_wrap(~Section, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("red", "blue")) +
  labs(colour = "Census") + 
  ggtitle("Comparison of Code Status Percentage in Old and New Censuses") +
  theme(strip.text.x = element_text(size = 28, color = "black", face = "bold",family = "mono"),
        strip.text.y = element_text(size = 28, color = "black", face = "bold",family = "mono"),
        legend.title = element_text(family = "mono",size = 24),
        legend.text = element_text(family = "mono",size = 24),
        legend.position="right",
        axis.title.y  = element_text(family = "mono",size = 24, face = "bold"),
        axis.title.x  = element_text(family = "mono",size = 24, face = "bold"),
        axis.text.x   = element_text(family = "mono",size = 24, face = "bold"),
        axis.text.y   = element_text(family = "mono",size = 24, face = "bold"),
        plot.title = element_text(family = "mono",size = 26,face = "bold", hjust = 0.5))
ggsave(paste0(pathwrite,"CodeStatusBoxPlot_Hosp6.jpg"), plot = p, width = 14, height = 8)

# Conduct statistical tests for comparing performance of new and old censuses
Tresult <- data.frame()
for (section in c("ICU","ED")){
  for (code in c("FC","LC","DNR","NS")){
    result <- list()
    g1 <- subset(dfSub6Long, Section == section & CodeType == code & Census == "Old")$Value
    g2 <- subset(dfSub6Long, Section == section & CodeType == code & Census == "New")$Value
    vartest <- var.test(g1,g2)
    ttest <- t.test(g2, g1, var.equal = vartest$p.value > 0.05)
    if (ttest$p.value > 0.1){
      result["Ttest"] <- "No Significant Difference"
    } else if (ttest$p.value > 0.05){
      result["Ttest"] <- ifelse(code == "NS", "Decreased (.)", "Increased (.)")
    } else if (ttest$p.value > 0.01){
      result["Ttest"] <- ifelse(code == "NS", "Decreased (*)", "Increased (*)")
    } else if (ttest$p.value > 0.001){
      result["Ttest"] <- ifelse(code == "NS", "Decreased (**)", "Increased (**)")
    } else {
      result["Ttest"] <- ifelse(code == "NS", "Decreased (***)", "Increased (***)")
    }
    result["EqualVariance"] <- ifelse(vartest$p.value <= 0.05,"No","Yes")
    result$CodeType <- code
    result$Section = section
    Tresult <- rbind(Tresult,result)
  }
}
tmp <- dfSub6Long %>%
  group_by(Section,CodeType,Census) %>%
  summarize(mean = round(mean(Value, na.rm = T),2),
            lower_CI = round(quantile(Value, 0.025, na.rm = T),2),
            upper_CI = round(quantile(Value, 0.975, na.rm = T),2)) %>%
  mutate(CodePct = paste0(mean," [",lower_CI, ", ",upper_CI,"]")) %>%
  ungroup() %>%
  select(Section,CodeType,Census,CodePct) %>%
  spread(Census , CodePct) %>%
  merge(Tresult, by = c("CodeType","Section")) %>%
  arrange(factor(Section, levels = c("ICU", "ED")), match(CodeType, c("FC", "DNR", "LC", "NS")))
fwrite(tmp, paste0(pathwrite,"CodeStatHosp6.csv"))
knitr::kable(tmp)
#   |CodeType |Section |Old                  |New                  |Ttest           |EqualVariance |
#   |:--------|:-------|:--------------------|:--------------------|:---------------|:-------------|
#   |FC       |ICU     |68.41 [58.42, 77.87] |90.11 [83.51, 95.83] |Increased (***) |No            |
#   |DNR      |ICU     |2.61 [0, 6.93]       |6.15 [1.96, 11.54]   |Increased (***) |No            |
#   |LC       |ICU     |1.56 [0, 4.49]       |3.16 [0, 7.69]       |Increased (***) |No            |
#   |NS       |ICU     |27.41 [17.35, 37.8]  |0.57 [0, 2.27]       |Decreased (***) |No            |
#   |FC       |ED      |33.87 [22.37, 45.45] |48.78 [36.36, 60.71] |Increased (***) |No            |
#   |DNR      |ED      |2.98 [0, 7.58]       |4.11 [0, 9.28]       |Increased (***) |No            |
#   |LC       |ED      |0.3 [0, 1.89]        |0.37 [0, 2.13]       |Increased (**)  |No            |
#   |NS       |ED      |62.85 [50.68, 75]    |46.74 [34.49, 59.37] |Decreased (***) |Yes           |


# ---------------------------------------------------------------------------- #
# -------------------------- All Code Transitions ---------------------------- #
# ---------------------------------------------------------------------------- #


# ------------------------ #
# Each hospital separately #
# ------------------------ #
tmp <- dfCodeTransition %>%
  group_by(Section,Code,HospitalCode) %>%
  summarize(mean = round(mean(ChangePct, na.rm = T),2),
            lower_CI = round(quantile(ChangePct, 0.025, na.rm = T),2),
            upper_CI = round(quantile(ChangePct, 0.975, na.rm = T),2))
tmp["95% CI"] <- paste0("(",tmp$lower_CI, ", ",tmp$upper_CI,")")
tmp <- tmp[,-c(5,6)]
fwrite(tmp, paste0(pathwrite,"CodeTransistionStat_AllHosp.csv"))
knitr::kable(tmp)
#   |Section |Code |HospitalCode |   mean|95% CI           |
#   |:-------|:----|:------------|------:|:----------------|
#   |ED      |DNR  |Hospital 1   |   0.78|(0, 7.69)        |
#   |ED      |DNR  |Hospital 2   |   0.73|(0, 5.41)        |
#   |ED      |DNR  |Hospital 3   |   0.45|(-1.82, 3.57)    |
#   |ED      |DNR  |Hospital 4   |   1.72|(0, 9.09)        |
#   |ED      |DNR  |Hospital 5   |   1.07|(-1.52, 5)       |
#   |ED      |DNR  |Hospital 6   |   1.13|(0, 4.65)        |
#   |ED      |FC   |Hospital 1   |  15.33|(0, 37.5)        |
#   |ED      |FC   |Hospital 2   |  17.49|(0, 34.78)       |
#   |ED      |FC   |Hospital 3   |  11.55|(-2.5, 27.58)    |
#   |ED      |FC   |Hospital 4   |  15.41|(0, 35)          |
#   |ED      |FC   |Hospital 5   |  16.15|(3.43, 29.72)    |
#   |ED      |FC   |Hospital 6   |  14.91|(4.45, 25.78)    |
#   |ED      |LC   |Hospital 1   |   0.07|(0, 0)           |
#   |ED      |LC   |Hospital 2   |   0.03|(0, 0)           |
#   |ED      |LC   |Hospital 3   |   0.04|(0, 1.52)        |
#   |ED      |LC   |Hospital 4   |   0.18|(0, 4)           |
#   |ED      |LC   |Hospital 5   |   0.08|(0, 1.72)        |
#   |ED      |LC   |Hospital 6   |   0.07|(0, 1.2)         |
#   |ED      |NS   |Hospital 1   | -16.18|(-39.13, 0)      |
#   |ED      |NS   |Hospital 2   | -18.24|(-35.71, -2.86)  |
#   |ED      |NS   |Hospital 3   | -12.04|(-28.26, 2.33)   |
#   |ED      |NS   |Hospital 4   | -17.31|(-37.5, 0)       |
#   |ED      |NS   |Hospital 5   | -17.29|(-30.95, -3.92)  |
#   |ED      |NS   |Hospital 6   | -16.11|(-26.98, -5.41)  |
#   |ICU     |DNR  |Hospital 1   |  10.56|(-8.33, 33.33)   |
#   |ICU     |DNR  |Hospital 2   |   6.56|(-20, 40)        |
#   |ICU     |DNR  |Hospital 3   |   4.67|(-16.67, 33.33)  |
#   |ICU     |DNR  |Hospital 4   |  11.91|(-12.5, 50)      |
#   |ICU     |DNR  |Hospital 5   |   5.32|(-1.75, 14.01)   |
#   |ICU     |DNR  |Hospital 6   |   3.54|(-0.98, 8.69)    |
#   |ICU     |FC   |Hospital 1   |  29.20|(-18.18, 72.73)  |
#   |ICU     |FC   |Hospital 2   |  27.31|(-20, 75.14)     |
#   |ICU     |FC   |Hospital 3   |  29.49|(-25, 80.33)     |
#   |ICU     |FC   |Hospital 4   |  21.70|(-25, 66.67)     |
#   |ICU     |FC   |Hospital 5   |  33.25|(16.36, 49.13)   |
#   |ICU     |FC   |Hospital 6   |  21.70|(10.75, 33.01)   |
#   |ICU     |LC   |Hospital 1   |   2.41|(-11.11, 20)     |
#   |ICU     |LC   |Hospital 2   |   1.55|(-6.67, 20)      |
#   |ICU     |LC   |Hospital 3   |   3.83|(-16.67, 33.33)  |
#   |ICU     |LC   |Hospital 4   |   5.58|(-12.5, 33.33)   |
#   |ICU     |LC   |Hospital 5   |   2.46|(-2.08, 8.06)    |
#   |ICU     |LC   |Hospital 6   |   1.60|(-2.06, 5.94)    |
#   |ICU     |NS   |Hospital 1   | -42.18|(-80, 0)         |
#   |ICU     |NS   |Hospital 2   | -35.42|(-77.84, 0)      |
#   |ICU     |NS   |Hospital 3   | -37.99|(-100, 0)        |
#   |ICU     |NS   |Hospital 4   | -39.19|(-80, 0)         |
#   |ICU     |NS   |Hospital 5   | -41.03|(-55.49, -26.58) |
#   |ICU     |NS   |Hospital 6   | -26.84|(-37.37, -16.67) |

# ---------------------- #
# All hospitals together #
# ---------------------- #
tmp <- dfCodeTransition %>%
  group_by(Section,Code) %>%
  summarize(mean = round(mean(ChangePct, na.rm = T),2),
            lower_CI = round(quantile(ChangePct, 0.025, na.rm = T),2),
            upper_CI = round(quantile(ChangePct, 0.975, na.rm = T),2))
tmp["95% CI"] <- paste0("(",tmp$lower_CI, ", ",tmp$upper_CI,")")
tmp <- tmp[,-c(4,5)]
fwrite(tmp, paste0(pathwrite,"CodeTransistionStat_AllHosp.csv"))
knitr::kable(tmp)
#   |Section |Code |   mean|95% CI        |
#   |:-------|:----|------:|:-------------|
#   |ED      |DNR  |   0.98|(-1.39, 6.25) |
#   |ED      |FC   |  15.14|(0, 33.33)    |
#   |ED      |LC   |   0.08|(0, 1.64)     |
#   |ED      |NS   | -16.19|(-34.1, 0)    |
#   |ICU     |DNR  |   7.10|(-10, 33.33)  |
#   |ICU     |FC   |  27.10|(-20, 75)     |
#   |ICU     |LC   |   2.91|(-6.81, 25)   |
#   |ICU     |NS   | -37.10|(-77.78, 0)   |

# -------------------------------- #
# Plot of each hospital separately #
# -------------------------------- #
p <- ggplot(dfCodeTransition, aes(x = HospitalCode, y = ChangePct, fill = Section)) +
  geom_boxplot() +
  labs(x='Code Status', y = "Change Percentage") + 
  facet_wrap(~Code, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("red", "blue")) +
  labs(colour = "Section") + 
  ggtitle("Percentage of code status transitions after updating the documentation system") +
  theme(strip.text.x = element_text(size = 20, color = "black", face = "bold",family = "mono"),
        strip.text.y = element_text(size = 20, color = "black", face = "bold",family = "mono"),
        legend.title = element_text(family = "mono",size = 18),
        legend.text = element_text(family = "mono",size = 18),
        legend.position="bottom",
        axis.title.y  = element_text(family = "mono",size = 18, face = "bold"),
        axis.title.x  = element_text(family = "mono",size = 18, face = "bold"),
        axis.text.x   = element_text(family = "mono",size = 16, face = "bold", angle = 90),
        axis.text.y   = element_text(family = "mono",size = 16, face = "bold"),
        plot.title = element_text(family = "mono",size = 20,face = "bold", hjust = 0.5))
ggsave(paste0(pathwrite,"CodeStatusTransitionBoxPlot_AllHosp.jpg"), plot = p, width = 14, height = 10)


# -------------------------- #
# T test for code transition #
# -------------------------- #
Tresult <- data.frame()
for (section in c("ICU","ED")){
  for (code in c("FC","LC","DNR","NS")){
    result <- list()
    nhosp <- 6
    for (ii in 1:nhosp){
      g1 <- subset(dfCodeTransition, HospitalCode == paste0("Hospital ",ii) & Section == section & Code == code)$ChangePct
      tmp <- t.test(g1, mu = 0, alternative = ifelse(code == "NS", "greater", "less"))
      result[paste0("Hospital ",ii)] <- tmp$p.value
    }
    tmp <- names(result)
    result <- as.data.frame(unlist(result))
    names(result) <- "Pvalue"
    result$group <- tmp
    result["Pvalue"] <- round(result["Pvalue"], 4)
    rownames(result) <- NULL
    result$Code <- code
    result$Section = section
    Tresult <- rbind(Tresult,result)
  }
}
fwrite(Tresult,paste0(pathwrite,"ttest_codetransitions.csv"))

# ------------------------------- #
# T test for comparing ICU and ED #
# ------------------------------- #
Tresult <- data.frame()
nhosp <- 6
for (ii in 1:nhosp){
  for (code in c("FC","LC","DNR","NS")){
    result <- list()
    g1 <- subset(dfCodeTransition, HospitalCode == paste0("Hospital ",ii) & Section == "ICU" & Code == code)$ChangePct
    g2 <- subset(dfCodeTransition, HospitalCode == paste0("Hospital ",ii) & Section == "ED" & Code == code)$ChangePct
    tmp <- t.test(g1, g2, alternative = ifelse(code == "NS", "greater", "less"))
    result[paste0("Hospital ",ii)] <- tmp$p.value
    tmp <- names(result)
    result <- as.data.frame(unlist(result))
    names(result) <- "Pvalue"
    result$group <- tmp
    result["Pvalue"] <- round(result["Pvalue"], 4)
    rownames(result) <- NULL
    result$Code <- code
    Tresult <- rbind(Tresult,result)
  }
}
fwrite(Tresult,paste0(pathwrite,"ttest_sectioncomparison.csv"))

# ----------------------------- #
# ANOVA for comparing hospitals #
# ----------------------------- #
print(leveneTest(ChangePct ~ HospitalCode, data = subset(dfCodeTransition, Code == "DNR"& Section == "ICU"), center = mean))
print(oneway.test(ChangePct ~ HospitalCode, data = subset(dfCodeTransition, Code == "DNR"& Section == "ICU"), var.equal = F))
TukeyResult <- data.frame()
KSresult <- data.frame()
for (section in c("ICU","ED")){
  for (code in c("FC","LC","DNR","NS")){
    result <- list()
    nhosp <- 6
    for (ii in 1:nhosp){
      if (ii < 6){
        for (jj in (ii+1):nhosp){
          g1 <- subset(dfCodeTransition, HospitalCode == paste0("Hospital ",ii) & Section == section & Code == code)$ChangePct
          g2 <- subset(dfCodeTransition, HospitalCode == paste0("Hospital ",jj) & Section == section & Code == code)$ChangePct
          tmp <- ks.test(g1, g2)
          result[paste0(ii,"-",jj)] <- tmp$p.value
        }
      }
    }
    tmp <- names(result)
    result <- as.data.frame(unlist(result))
    names(result) <- "Padj"
    result$group <- tmp
    result["Padj"] <- round(result["Padj"], 4)
    rownames(result) <- NULL
    result$Code <- code
    result$Section = section
    KSresult <- rbind(KSresult,result)
  }
}
KSresultSig <- KSresult[which(KSresult["Padj"] > 0.05 & KSresult["Padj"] < 1),]


# ---------------------------------------------------------------------------- #
# --------------------------- Only NS Transitions ---------------------------- #
# ---------------------------------------------------------------------------- #

# ---------------- #
# Days without NS  #
# ---------------- #
# Looking at days without "Not Specified" code in hospitals under two coding systems
knitr::kable(dfSub %>% group_by(Hospital_Code) %>% 
               summarise(ICU_NS_Old_Count = round(sum(ifelse(ICU_NS_Old == 0,0,1)) / 1424*100,2),
                         ICU_NS_New_Count = round(sum(ifelse(ICU_NS_New == 0,0,1)) / 1424*100,2),
                         ED_NS_Old_Count = round(sum(ifelse(ED_NS_Old == 0,0,1)) / 1424*100,2),
                         ED_NS_New_Count = round(sum(ifelse(ED_NS_New == 0,0,1)) / 1424*100,2)))
# | Hospital_Code| ICU_NS_Old_Count| ICU_NS_New_Count| ED_NS_Old_Count| ED_NS_New_Count|
# |-------------:|----------------:|----------------:|---------------:|---------------:|
# |             1|            95.65|            10.39|           98.46|           98.46|
# |             2|            88.34|             8.71|           98.46|           98.46|
# |             3|            83.92|             1.19|           98.46|           98.46|
# |             4|            92.77|             1.69|           98.31|           98.38|
# |             5|            97.47|            15.45|           97.47|           97.47|
# |             6|            98.46|            41.64|           98.46|           98.46|


# ------------------- #
# ICU NS  transitions #
# ------------------- #
tmpICU <- dfNSTransitionICU %>% group_by(Hospital_Code,CovidStatus) %>% 
  summarise(ICU_NS2FCPct = paste0(round(mean(ICU_NS2FCPct, na.rm = T),2),
                                  " [",round(sd(ICU_NS2FCPct, na.rm = T),2),"]"),
            ICU_NS2DNRPct = paste0(round(mean(ICU_NS2DNRPct, na.rm = T),2),
                                  " [",round(sd(ICU_NS2DNRPct, na.rm = T),2),"]"),
            ICU_NS2LCPct = paste0(round(mean(ICU_NS2LCPct, na.rm = T),2),
                                  " [",round(sd(ICU_NS2LCPct, na.rm = T),2),"]"),
            ICU_NS2NSPct = paste0(round(mean(ICU_NS2NSPct, na.rm = T),2),
                                  " [",round(sd(ICU_NS2NSPct, na.rm = T),2),"]"))
tmpICU_AllTime <- dfNSTransitionICU %>%
  group_by(Hospital_Code) %>%
  summarise(
    ICU_NS2FCPct = paste0(round(mean(ICU_NS2FCPct, na.rm = TRUE), 2),
                          " [", round(sd(ICU_NS2FCPct, na.rm = TRUE), 2), "]"),
    ICU_NS2DNRPct = paste0(round(mean(ICU_NS2DNRPct, na.rm = TRUE), 2),
                           " [", round(sd(ICU_NS2DNRPct, na.rm = TRUE), 2), "]"),
    ICU_NS2LCPct = paste0(round(mean(ICU_NS2LCPct, na.rm = TRUE), 2),
                          " [", round(sd(ICU_NS2LCPct, na.rm = TRUE), 2), "]"),
    ICU_NS2NSPct = paste0(round(mean(ICU_NS2NSPct, na.rm = TRUE), 2),
                          " [", round(sd(ICU_NS2NSPct, na.rm = TRUE), 2), "]"),
    CovidStatus = "All Time"
  )
tmpICU <- rbind(tmpICU, tmpICU_AllTime) %>%
  arrange(Hospital_Code, CovidStatus)
tmpICU_AllTime_AllHosp <- dfNSTransitionICU %>%
  summarise(
    ICU_NS2FCPct = paste0(round(mean(ICU_NS2FCPct, na.rm = TRUE), 2),
                          " [", round(sd(ICU_NS2FCPct, na.rm = TRUE), 2), "]"),
    ICU_NS2DNRPct = paste0(round(mean(ICU_NS2DNRPct, na.rm = TRUE), 2),
                           " [", round(sd(ICU_NS2DNRPct, na.rm = TRUE), 2), "]"),
    ICU_NS2LCPct = paste0(round(mean(ICU_NS2LCPct, na.rm = TRUE), 2),
                          " [", round(sd(ICU_NS2LCPct, na.rm = TRUE), 2), "]"),
    ICU_NS2NSPct = paste0(round(mean(ICU_NS2NSPct, na.rm = TRUE), 2),
                          " [", round(sd(ICU_NS2NSPct, na.rm = TRUE), 2), "]"),
    CovidStatus = "All Time",
    Hospital_Code = "All"
  )
tmpICU$Hospital_Code <- as.character(tmpICU$Hospital_Code)
tmpICU <- rbind(tmpICU, tmpICU_AllTime_AllHosp) %>%
  arrange(Hospital_Code, CovidStatus)
fwrite(tmpICU, paste0(pathwrite,"NStransitionStat_ICU.csv"))
knitr::kable(tmpICU)
#   |Hospital_Code |CovidStatus |ICU_NS2FCPct  |ICU_NS2DNRPct |ICU_NS2LCPct  |ICU_NS2NSPct |
#   |:-------------|:-----------|:-------------|:-------------|:-------------|:------------|
#   |1             |All Time    |65.79 [30.62] |25.31 [27.95] |6.54 [14.71]  |2.36 [7.43]  |
#   |1             |in-Covid    |70.17 [26.85] |20.61 [23.25] |6.7 [13.81]   |2.52 [7.37]  |
#   |1             |post-Covid  |57.44 [35.39] |34.7 [35.11]  |6.2 [14.51]   |1.67 [6.87]  |
#   |1             |pre-Covid   |61.58 [33.33] |29.58 [30.1]  |6.4 [16.58]   |2.44 [7.85]  |
#   |2             |All Time    |70.1 [32.96]  |23.33 [31.39] |4.63 [14.23]  |1.94 [7.15]  |
#   |2             |in-Covid    |73.03 [29.03] |19.8 [27.4]   |5.07 [14.69]  |2.1 [6.35]   |
#   |2             |post-Covid  |65.21 [35.23] |29.57 [34.31] |2.67 [9.04]   |2.55 [9.12]  |
#   |2             |pre-Covid   |66.29 [39.45] |27.66 [37.1]  |5.01 [16.14]  |1.05 [7.31]  |
#   |3             |All Time    |77.85 [34.16] |12.09 [27.44] |9.49 [24.23]  |0.57 [4.98]  |
#   |3             |in-Covid    |77.02 [34.01] |12.69 [27.76] |9.79 [23.71]  |0.5 [4.69]   |
#   |3             |post-Covid  |81.06 [32.69] |10.56 [24.81] |8.38 [24.59]  |0 [0]        |
#   |3             |pre-Covid   |77.51 [35.33] |11.83 [28.38] |9.58 [25.09]  |1.08 [6.73]  |
#   |4             |All Time    |56.1 [35.01]  |28.24 [31.22] |15.12 [25.01] |0.54 [4.35]  |
#   |4             |in-Covid    |54.05 [33.95] |26.52 [29.17] |19.1 [25.79]  |0.34 [2.95]  |
#   |4             |post-Covid  |49.64 [32.57] |35.23 [31.11] |14.86 [24.36] |0.27 [3.68]  |
#   |4             |pre-Covid   |63.74 [37.11] |27.42 [34.42] |7.74 [22.08]  |1.1 [6.4]    |
#   |5             |All Time    |78.29 [11.75] |13.72 [9.28]  |7.21 [6.63]   |0.77 [2.07]  |
#   |5             |in-Covid    |78.9 [10.75]  |12.64 [8.5]   |7.97 [7.27]   |0.5 [1.48]   |
#   |5             |post-Covid  |74.78 [15.35] |16.92 [10.93] |7 [6.37]      |1.3 [3.3]    |
#   |5             |pre-Covid   |79.25 [10.64] |13.79 [9.17]  |5.99 [5.24]   |0.96 [1.97]  |
#   |6             |All Time    |76.78 [10.74] |13.53 [8.52]  |7.57 [6.41]   |2.11 [2.83]  |
#   |6             |in-Covid    |77.61 [10.75] |12.86 [8.14]  |7.59 [6.6]    |1.94 [2.69]  |
#   |6             |post-Covid  |76.69 [11.81] |13.47 [9.36]  |7.89 [6.21]   |1.95 [3.01]  |
#   |6             |pre-Covid   |75.17 [9.92]  |14.92 [8.66]  |7.38 [6.13]   |2.54 [2.97]  |
#   |All           |All Time    |70.95 [28.65] |19.28 [25.19] |8.37 [16.95]  |1.4 [5.23]   |
  
  
# ------------------ #
# ED NS  transitions #
# ------------------ #
tmpED <- dfNSTransitionED %>% group_by(Hospital_Code,CovidStatus) %>% 
  summarise(ED_NS2FCPct = paste0(round(mean(ED_NS2FCPct, na.rm = T),2),
                                  " [",round(sd(ED_NS2FCPct, na.rm = T),2),"]"),
            ED_NS2DNRPct = paste0(round(mean(ED_NS2DNRPct, na.rm = T),2),
                                   " [",round(sd(ED_NS2DNRPct, na.rm = T),2),"]"),
            ED_NS2LCPct = paste0(round(mean(ED_NS2LCPct, na.rm = T),2),
                                  " [",round(sd(ED_NS2LCPct, na.rm = T),2),"]"),
            ED_NS2NSPct = paste0(round(mean(ED_NS2NSPct, na.rm = T),2),
                                  " [",round(sd(ED_NS2NSPct, na.rm = T),2),"]"))
tmpED_AllTime <- dfNSTransitionED %>%
  group_by(Hospital_Code) %>%
  summarise(
    ED_NS2FCPct = paste0(round(mean(ED_NS2FCPct, na.rm = TRUE), 2),
                          " [", round(sd(ED_NS2FCPct, na.rm = TRUE), 2), "]"),
    ED_NS2DNRPct = paste0(round(mean(ED_NS2DNRPct, na.rm = TRUE), 2),
                           " [", round(sd(ED_NS2DNRPct, na.rm = TRUE), 2), "]"),
    ED_NS2LCPct = paste0(round(mean(ED_NS2LCPct, na.rm = TRUE), 2),
                          " [", round(sd(ED_NS2LCPct, na.rm = TRUE), 2), "]"),
    ED_NS2NSPct = paste0(round(mean(ED_NS2NSPct, na.rm = TRUE), 2),
                          " [", round(sd(ED_NS2NSPct, na.rm = TRUE), 2), "]"),
    CovidStatus = "All Time"
  )
tmpED <- rbind(tmpED, tmpED_AllTime) %>%
  arrange(Hospital_Code, CovidStatus)
tmpED_AllTime_AllHosp <- dfNSTransitionED %>%
  summarise(
    ED_NS2FCPct = paste0(round(mean(ED_NS2FCPct, na.rm = TRUE), 2),
                          " [", round(sd(ED_NS2FCPct, na.rm = TRUE), 2), "]"),
    ED_NS2DNRPct = paste0(round(mean(ED_NS2DNRPct, na.rm = TRUE), 2),
                           " [", round(sd(ED_NS2DNRPct, na.rm = TRUE), 2), "]"),
    ED_NS2LCPct = paste0(round(mean(ED_NS2LCPct, na.rm = TRUE), 2),
                          " [", round(sd(ED_NS2LCPct, na.rm = TRUE), 2), "]"),
    ED_NS2NSPct = paste0(round(mean(ED_NS2NSPct, na.rm = TRUE), 2),
                          " [", round(sd(ED_NS2NSPct, na.rm = TRUE), 2), "]"),
    CovidStatus = "All Time",
    Hospital_Code = "All"
  )
tmpED$Hospital_Code <- as.character(tmpED$Hospital_Code)
tmpED <- rbind(tmpED, tmpED_AllTime_AllHosp) %>%
  arrange(Hospital_Code, CovidStatus)
fwrite(tmpED, paste0(pathwrite,"NStransitionStat_ED.csv"))
knitr::kable(tmpED)
#   |Hospital_Code |CovidStatus |ED_NS2FCPct   |ED_NS2DNRPct |ED_NS2LCPct |ED_NS2NSPct   |
#   |:-------------|:-----------|:-------------|:------------|:-----------|:-------------|
#   |1             |All Time    |23.15 [12.05] |1.31 [3.42]  |0.13 [1.08] |75.4 [12.16]  |
#   |1             |in-Covid    |24.44 [11.63] |1.33 [3.22]  |0.16 [1.17] |74.07 [11.78] |
#   |1             |post-Covid  |27.58 [13.19] |1.78 [4.39]  |0.17 [1.21] |70.47 [12.97] |
#   |1             |pre-Covid   |16.78 [9.42]  |0.88 [2.98]  |0.06 [0.69] |82.28 [9.21]  |
#   |2             |All Time    |27.19 [11.91] |1.3 [2.66]   |0.06 [0.57] |71.45 [12.04] |
#   |2             |in-Covid    |29.19 [11.2]  |1.51 [2.87]  |0.05 [0.5]  |69.25 [11.22] |
#   |2             |post-Covid  |33.3 [11.08]  |1.63 [2.61]  |0.11 [0.71] |64.96 [10.91] |
#   |2             |pre-Covid   |19.39 [9.73]  |0.66 [2.1]   |0.07 [0.6]  |79.88 [9.58]  |
#   |3             |All Time    |21.21 [11.03] |1.05 [1.93]  |0.13 [0.63] |77.61 [11.26] |
#   |3             |in-Covid    |22.17 [8.38]  |1.11 [2.04]  |0.13 [0.61] |76.6 [8.42]   |
#   |3             |post-Covid  |30.53 [8.82]  |1.33 [1.92]  |0.16 [0.71] |67.98 [8.92]  |
#   |3             |pre-Covid   |5.72 [4.29]   |0.5 [1.28]   |0.08 [0.56] |93.7 [4.2]    |
#   |4             |All Time    |23.31 [11.57] |2.79 [4.2]   |0.35 [1.61] |73.55 [11.58] |
#   |4             |in-Covid    |25.03 [11.3]  |2.82 [4.13]  |0.33 [1.67] |71.81 [11.2]  |
#   |4             |post-Covid  |26.17 [11.46] |3.06 [4.17]  |0.65 [1.9]  |70.13 [11.53] |
#   |4             |pre-Covid   |17.74 [10.32] |2.52 [4.37]  |0.2 [1.16]  |79.53 [10.19] |
#   |5             |All Time    |25.45 [10.1]  |1.85 [2.47]  |0.15 [0.71] |72.55 [10.41] |
#   |5             |in-Covid    |28.17 [9.06]  |2.05 [2.53]  |0.18 [0.78] |69.6 [9.11]   |
#   |5             |post-Covid  |31.89 [7.9]   |2.05 [2.56]  |0.17 [0.67] |65.89 [7.87]  |
#   |5             |pre-Covid   |16.83 [7.15]  |1.34 [2.24]  |0.1 [0.59]  |81.73 [7.58]  |
#   |6             |All Time    |23.51 [8.05]  |1.88 [2.15]  |0.13 [0.52] |74.48 [8.23]  |
#   |6             |in-Covid    |25.46 [7.09]  |2.07 [2.27]  |0.16 [0.58] |72.32 [6.91]  |
#   |6             |post-Covid  |28.88 [6.55]  |2.12 [1.86]  |0.16 [0.53] |68.84 [6.45]  |
#   |6             |pre-Covid   |16.66 [5.9]   |1.38 [1.96]  |0.05 [0.35] |81.91 [6.24]  |
#   |All           |All Time    |24.07 [11]    |1.71 [2.97]  |0.16 [0.94] |74.06 [11.17] |


# ------------------------------ #
# Bar plot of ICU NS transitions #
# ------------------------------ #
ColList <- c("Date","Hospital_Code",names(dfNSTransitionICU)[grep("ICU_NS2",names(dfNSTransitionICU))])
tmp <- dfNSTransitionICU[,ColList]
tmp$Date <- floor_date(as.Date(tmp$Date), unit = "week")
names(tmp)[which(names(tmp) == "Hospital_Code")] = "HospitalCode"
tmp <- tmp %>% group_by(Date,HospitalCode) %>% summarize_all(mean)
tmp <- gather(tmp, TrnasitionType, ChangePct, 3:6, factor_key=TRUE)
tmp$TrnasitionType <- word(tmp$TrnasitionType, 2, sep = "_")

p <- ggplot(tmp, aes(x = Date, y = 100, fill = TrnasitionType)) + 
  geom_bar(aes(y = (ChangePct)), stat = "identity", width = 7) +
  facet_wrap(~HospitalCode, labeller = labeller(HospitalCode = function(x) {paste("Hospital", x)}))+
  geom_vline(xintercept = as.Date("2020-04-01",format = "%Y-%m-%d"), color = "red", size = 1) +
  geom_vline(xintercept = as.Date("2022-01-01",format = "%Y-%m-%d"), color = "blue", size = 1) +
  scale_fill_manual(values = c("NS2FCPct" = "#FFB6C1", "NS2DNRPct" = "#ADD8E6", "NS2LCPct" = "#90EE90", "NS2NSPct" = "black"),
                    labels = c("NS2FCPct" ="FC", "NS2DNRPct"="DNR", "NS2LCPct"="LC", "NS2NSPct"="NS")) +
  labs(x = "Date", y = "Percentage (%)", fill = "Transitions") +
  theme(strip.text.x = element_text(size = 20, color = "black", face = "bold",family = "mono"),
        strip.text.y = element_text(size = 20, color = "black", face = "bold",family = "mono"),
        legend.title = element_text(family = "mono",size = 16),
        legend.text = element_text(family = "mono",size = 16),
        legend.position="bottom",
        axis.title.y  = element_text(family = "mono",size = 16, face = "bold"),
        axis.title.x  = element_text(family = "mono",size = 16, face = "bold"),
        axis.text.x   = element_text(family = "mono",size = 14, face = "bold", angle = 45),
        axis.text.y   = element_text(family = "mono",size = 14, face = "bold"),
        plot.title = element_text(family = "mono",size = 20,face = "bold", hjust = 0.5))
ggsave(paste0(pathwrite,"NSTransitionBarPlot_AllHosp.jpg"), plot = p, width = 14, height = 8)


# --------------------------------------------- #
# test of COVID effect on NS transitions in ICU #
# --------------------------------------------- #
ColList <- c("Date","Hospital_Code",names(dfNSTransitionICU)[grep("ICU_NS2",names(dfNSTransitionICU))])
tmp <- dfNSTransitionICU[,ColList]
names(tmp)[which(names(tmp) == "Hospital_Code")] = "HospitalCode"
tmp <- gather(tmp, TrnasitionType, ChangePct, 3:6, factor_key=TRUE)
tmp$TrnasitionType <- word(tmp$TrnasitionType, 2, sep = "_")
tmp$CovidStatus <- "in-Covid"
tmp$CovidStatus[which(tmp$Date < as.Date("2020-04-01"))] <- "pre-Covid"
tmp$CovidStatus[which(tmp$Date > as.Date("2022-05-01"))] <- "post-Covid"
tmp$CovidStatus <- factor(tmp$CovidStatus,levels = c("pre-Covid","in-Covid","post-Covid"))
TukeyResult <- data.frame()
for (type in c("FC","DNR","LC","NS")){
  result <- list()
  for (ii in 1:6){
    res <- oneway.test(ChangePct ~ CovidStatus, data = subset(tmp, TrnasitionType == paste0("NS2",type,"Pct") &
                                                                   HospitalCode == ii), var.equal = F)
    if (is.na(res$p.value)){
      result[paste0(ii,"-","pre_in")] <- NaN
      result[paste0(ii,"-","in_post")] <- NaN
      result[paste0(ii,"-","pre_post")] <- NaN
    }else{
        g1 <- subset(tmp, HospitalCode == ii & TrnasitionType  == paste0("NS2",type,"Pct") & CovidStatus == "pre-Covid")$ChangePct
        g2 <- subset(tmp, HospitalCode == ii & TrnasitionType  == paste0("NS2",type,"Pct") & CovidStatus == "in-Covid")$ChangePct
        res <- tukeySEM(mean(g1,na.rm=T), mean(g2,na.rm=T), var(g1,na.rm=T), var(g2,na.rm=T), length(g1), length(g2), 3)["p"]
        result[paste0(ii,"-","pre_in")] <- paste0(round(res,3), "_",mean(g2)-mean(g1))
        
        g1 <- subset(tmp, HospitalCode == ii & TrnasitionType  == paste0("NS2",type,"Pct") & CovidStatus == "in-Covid")$ChangePct
        g2 <- subset(tmp, HospitalCode == ii & TrnasitionType  == paste0("NS2",type,"Pct") & CovidStatus == "post-Covid")$ChangePct
        res <- tukeySEM(mean(g1,na.rm=T), mean(g2,na.rm=T), var(g1,na.rm=T), var(g2,na.rm=T), length(g1), length(g2), 3)["p"]
        result[paste0(ii,"-","in_post")] <- paste0(round(res,3), "_",mean(g2)-mean(g1))

        g1 <- subset(tmp, HospitalCode == ii & TrnasitionType  == paste0("NS2",type,"Pct") & CovidStatus == "pre-Covid")$ChangePct
        g2 <- subset(tmp, HospitalCode == ii & TrnasitionType  == paste0("NS2",type,"Pct") & CovidStatus == "post-Covid")$ChangePct
        res <- tukeySEM(mean(g1,na.rm=T), mean(g2,na.rm=T), var(g1,na.rm=T), var(g2,na.rm=T), length(g1), length(g2), 3)["p"]
        result[paste0(ii,"-","pre_post")] <- paste0(round(res,3), "_",mean(g2)-mean(g1))
    }
  }
  result <- as.data.frame(unlist(result))
  result$HospitalCode <- word(rownames(result), 1, sep = "-")
  result$Comparison <- word(rownames(result), 2, sep = "-")
  names(result)[1] <- "Value"
  result$Pvalue <- ifelse(result$Value == "NaN",NaN,word(result$Value,1,sep="_"))
  result$Change <- ifelse(result$Value == "NaN",NaN,word(result$Value,2,sep="_"))
  result <- select(result,-c("Value"))
  result <- reshape(result, idvar = "HospitalCode", timevar = "Comparison", direction = "wide")
  result$`Chosen Code` <- type
  rownames(result) <- NULL
  TukeyResult <- rbind(TukeyResult, result)
}
fwrite(TukeyResult,paste0(pathwrite,"ttest_NStransitions_ICU_MayCutoff.csv"))
knitr::kable(TukeyResult, digits = 2)
#   |HospitalCode |Pvalue.pre_in |Change.pre_in      |Pvalue.in_post |Change.in_post      |Pvalue.pre_post |Change.pre_post    |Chosen Code |
#   |:------------|:-------------|:------------------|:--------------|:-------------------|:---------------|:------------------|:-----------|
#   |1            |0             |8.59496143739855   |0              |-12.7366483065656   |0.386           |-4.14168686916702  |FC          |
#   |2            |0.039         |6.74403121689211   |0.014          |-7.82259014648805   |0.95            |-1.07855892959594  |FC          |
#   |3            |0.979         |-0.488604927557063 |0.311          |4.04096852838852    |0.495           |3.55236360083146   |FC          |
#   |4            |0             |-9.69402180851621  |0.249          |-4.41036200926384   |0               |-14.10438381778    |FC          |
#   |5            |0.864         |-0.358414098194814 |0.001          |-4.11873712813636   |0.001           |-4.47715122633117  |FC          |
#   |6            |0.001         |2.43780752026218   |0.614          |-0.917184800527522  |0.307           |1.52062271973466   |FC          |
#   |1            |0             |-8.96874756741729  |0              |14.0884667922751    |0.21            |5.11971922485779   |DNR         |
#   |2            |0.007         |-7.85472891566265  |0.001          |9.76974489795919    |0.839           |1.91501598229653   |DNR         |
#   |3            |0.903         |0.857492108943566  |0.58           |-2.12885920969507   |0.861           |-1.27136710075151  |DNR         |
#   |4            |0.916         |-0.903044297240449 |0.002          |8.71325745379323    |0.025           |7.81021315655278   |DNR         |
#   |5            |0.118         |-1.15240482764911  |0              |4.2826372587077     |0.001           |3.13023243105859   |DNR         |
#   |6            |0.001         |-2.05378860888436  |0.704          |0.611685789647215   |0.203           |-1.44210281923715  |DNR         |
#   |1            |0.957         |0.298881953599325  |0.904          |-0.500907115615456  |0.988           |-0.20202516201613  |LC          |
#   |2            |0.998         |0.0626584867075666 |0.016          |-2.3993761737824    |0.131           |-2.33671768707483  |LC          |
#   |3            |0.992         |0.211579206306212  |0.768          |-1.41287470970842   |0.861           |-1.20129550340221  |LC          |
#   |4            |0             |11.3567456231229   |0.104          |-4.23745019561991   |0.003           |7.11929542750297   |LC          |
#   |5            |0             |1.97574214118096   |0.149          |-0.970224159226453  |0.124           |1.0055179819545    |LC          |
#   |6            |0.869         |0.212742061207317  |0.838          |0.299184141114408   |0.643           |0.511926202321725  |LC          |
#   |1            |0.988         |0.0752929412503249 |0.294          |-0.851168462157932  |0.466           |-0.775875520907607 |NS          |
#   |2            |0.115         |1.04831891985118   |0.794          |0.451664892951045   |0.147           |1.49998381280223   |NS          |
#   |3            |NaN           |NaN                |NaN            |NaN                 |NaN             |NaN                |NS          |
#   |4            |0.112         |-0.759805146558671 |0.974          |-0.0652704951207447 |0.159           |-0.825075641679415 |NS          |
#   |5            |0             |-0.464668399711436 |0.002          |0.806291350531108   |0.353           |0.341622950819672  |NS          |
#   |6            |0.006         |-0.596210062447407 |1              |0.00570342894823606 |0.084           |-0.590506633499171 |NS          |
  
  
# -------------------------------------------- #
# test of COVID effect on NS transitions in ED #
# -------------------------------------------- #
ColList <- c("Date","Hospital_Code",names(dfNSTransitionED)[grep("ED_NS2",names(dfNSTransitionED))])
tmp <- dfNSTransitionED[,ColList]
names(tmp)[which(names(tmp) == "Hospital_Code")] = "HospitalCode"
tmp <- gather(tmp, TrnasitionType, ChangePct, 3:6, factor_key=TRUE)
tmp$TrnasitionType <- word(tmp$TrnasitionType, 2, sep = "_")
tmp$CovidStatus <- "in-Covid"
tmp$CovidStatus[which(tmp$Date < as.Date("2020-04-01"))] <- "pre-Covid"
tmp$CovidStatus[which(tmp$Date > as.Date("2022-05-01"))] <- "post-Covid"
tmp$CovidStatus <- factor(tmp$CovidStatus,levels = c("pre-Covid","in-Covid","post-Covid"))
TukeyResult <- data.frame()
for (type in c("FC","DNR","LC","NS")){
  result <- list()
  for (ii in 1:6){
    res <- oneway.test(ChangePct ~ CovidStatus, data = subset(tmp, TrnasitionType == paste0("NS2",type,"Pct") &
                                                                HospitalCode == ii), var.equal = F)
    if (is.na(res$p.value)){
      result[paste0(ii,"-","pre_in")] <- NaN
      result[paste0(ii,"-","in_post")] <- NaN
      result[paste0(ii,"-","pre_post")] <- NaN
    }else{
        g1 <- subset(tmp, HospitalCode == ii & TrnasitionType  == paste0("NS2",type,"Pct") & CovidStatus == "pre-Covid")$ChangePct
        g2 <- subset(tmp, HospitalCode == ii & TrnasitionType  == paste0("NS2",type,"Pct") & CovidStatus == "in-Covid")$ChangePct
        res <- tukeySEM(mean(g1,na.rm=T), mean(g2,na.rm=T), var(g1,na.rm=T), var(g2,na.rm=T), length(g1), length(g2), 3)["p"]
        result[paste0(ii,"-","pre_in")] <- paste0(round(res,3), "_",mean(g2)-mean(g1))
        
        g1 <- subset(tmp, HospitalCode == ii & TrnasitionType  == paste0("NS2",type,"Pct") & CovidStatus == "in-Covid")$ChangePct
        g2 <- subset(tmp, HospitalCode == ii & TrnasitionType  == paste0("NS2",type,"Pct") & CovidStatus == "post-Covid")$ChangePct
        res <- tukeySEM(mean(g1,na.rm=T), mean(g2,na.rm=T), var(g1,na.rm=T), var(g2,na.rm=T), length(g1), length(g2), 3)["p"]
        result[paste0(ii,"-","in_post")] <- paste0(round(res,3), "_",mean(g2)-mean(g1))

        g1 <- subset(tmp, HospitalCode == ii & TrnasitionType  == paste0("NS2",type,"Pct") & CovidStatus == "pre-Covid")$ChangePct
        g2 <- subset(tmp, HospitalCode == ii & TrnasitionType  == paste0("NS2",type,"Pct") & CovidStatus == "post-Covid")$ChangePct
        res <- tukeySEM(mean(g1,na.rm=T), mean(g2,na.rm=T), var(g1,na.rm=T), var(g2,na.rm=T), length(g1), length(g2), 3)["p"]
        result[paste0(ii,"-","pre_post")] <- paste0(round(res,3), "_",mean(g2)-mean(g1))
    }
  }
  result <- as.data.frame(unlist(result))
  result$HospitalCode <- word(rownames(result), 1, sep = "-")
  result$Comparison <- word(rownames(result), 2, sep = "-")
  names(result)[1] <- "Value"
  result$Pvalue <- ifelse(result$Value == "NaN",NaN,word(result$Value,1,sep="_"))
  result$Change <- ifelse(result$Value == "NaN",NaN,word(result$Value,2,sep="_"))
  result <- select(result,-c("Value"))
  result <- reshape(result, idvar = "HospitalCode", timevar = "Comparison", direction = "wide")
  result$`Chosen Code` <- type
  rownames(result) <- NULL
  TukeyResult <- rbind(TukeyResult, result)
}
fwrite(TukeyResult,paste0(pathwrite,"ttest_NStransitions_ED_MayCutoff.csv"))
knitr::kable(TukeyResult)
#   |HospitalCode |Pvalue.pre_in |Change.pre_in       |Pvalue.in_post |Change.in_post       |Pvalue.pre_post |Change.pre_post    |Chosen Code |
#   |:------------|:-------------|:-------------------|:--------------|:--------------------|:---------------|:------------------|:-----------|
#   |1            |0             |7.66139931501736    |0.004          |3.14170714427762     |0               |10.803106459295    |FC          |
#   |2            |0             |9.79491036901641    |0              |4.11186312451915     |0               |13.9067734935356   |FC          |
#   |3            |0             |16.4482914930051    |0              |8.36631528095601     |0               |24.8146067739611   |FC          |
#   |4            |0             |7.2895332160216     |0.396          |1.13375901009885     |0               |8.42329222612045   |FC          |
#   |5            |0             |11.3330855229355    |0              |3.72741199886641     |0               |15.0604975218019   |FC          |
#   |6            |0             |8.79931732789913    |0              |3.42353051539565     |0               |12.2228478432948   |FC          |
#   |1            |0.08          |0.451964612340886   |0.341          |0.443409068310629    |0.023           |0.895373680651515  |DNR         |
#   |2            |0             |0.849933040928084   |0.821          |0.120032581073835    |0               |0.969965622001919  |DNR         |
#   |3            |0             |0.612081758263061   |0.314          |0.216265827612509    |0               |0.82834758587557   |DNR         |
#   |4            |0.547         |0.299433016990271   |0.746          |0.231858892582654    |0.315           |0.531291909572925  |DNR         |
#   |5            |0             |0.719156109920094   |1              |-0.00511194152809891 |0.002           |0.714044168391995  |DNR         |
#   |6            |0             |0.684212712362111   |0.925          |0.0567528725061006   |0               |0.740965584868211  |DNR         |
#   |1            |0.202         |0.101331423345403   |0.985          |0.015100676145968    |0.398           |0.116432099491371  |LC          |
#   |2            |0.81          |-0.0224859389656566 |0.424          |0.0626321311169655   |0.754           |0.0401461921513089 |LC          |
#   |3            |0.597         |0.0464259636199154  |0.758          |0.0369944063056191   |0.377           |0.0834203699255345 |LC          |
#   |4            |0.347         |0.124363526869229   |0.059          |0.322597485027682    |0.005           |0.446961011896911  |LC          |
#   |5            |0.164         |0.0757810102585     |0.993          |-0.00602105041158366 |0.416           |0.0697599598469164 |LC          |
#   |6            |0             |0.112938398182937   |0.998          |-0.00256943678034063 |0.017           |0.110368961402597  |LC          |
#   |1            |0             |-8.21459055022585   |0.001          |-3.60014455105896    |0               |-11.8147351012848  |NS          |
#   |2            |0             |-10.6222490981288   |0              |-4.29447060126361    |0               |-14.9167196993924  |NS          |
#   |3            |0             |-17.1068284546541   |0              |-8.61967670989067    |0               |-25.7265051645448  |NS          |
#   |4            |0             |-7.713178649833     |0.132          |-1.68823090990216    |0               |-9.40140955973516  |NS          |
#   |5            |0             |-12.1278843748339   |0              |-3.71656998240307    |0               |-15.844454357237   |NS          |
#   |6            |0             |-9.59646182258167   |0              |-3.47685906632822    |0               |-13.0733208889099  |NS          |

# ---------------------------------------------------------------------------- #
# --------------------------------- Regression ------------------------------- #
# ---------------------------------------------------------------------------- #

# Simple correlation plots
for (depvar in c("ICU_FC","ICU_DNR","ICU_LC","ICU_NS","ED_FC","ED_DNR","ED_LC","ED_NS")){
  color1 <- "red"
  color2 <- "blue"
  counter = 1
  fig_list <- list()
  cor_list <- list()
  for (grp in unique(dfReg$Hospital_Code)){
    tmp <- subset(dfReg, Hospital_Code == grp)
    tmp$Date <- floor_date(as.Date(tmp$Date), unit = "quarter")
    tmp <- tmp %>% group_by(Date) %>% summarize_all(mean)
    tmp <- tmp[,c("Date","Vent2",depvar)]
    names(tmp) <- c("Date","Vent","Value")
    res <- cor.test(tmp$Vent,tmp$Value, use = "complete.obs")
    if (res$p.value > 0.1){
      res <- paste0(round(res$estimate,2))
    } else if (res$p.value > 0.05){
      res <- paste0(round(res$estimate,2),"(.)")
    } else if (res$p.value > 0.01){
      res <- paste0(round(res$estimate,2),"(*)")
    } else if (res$p.value > 0.001){
      res <- paste0(round(res$estimate,2),"(**)")
    } else{
      res <- paste0(round(res$estimate,2),"(***)")
    }
    cor_list[paste0(grp,"_",depvar)] <- res 
    textdata <- as.data.frame(tmp %>% summarize(ylevel = quantile(Value,0.90), yvalue = res))
    sec <- with(tmp, train_sec(Vent, Value))
    fig <- ggplot(tmp, aes(x = Date)) +
      geom_line(aes(y = Vent), color = color1, size=0.5, linetype = "dashed") +
      geom_line(aes(y = sec$fwd(Value)), color = color2, size=0.5, linetype = "dashed") +
      geom_label(data = textdata, aes(label = paste0("cor = ",yvalue),
                                      y = ylevel, x = as.Date("2021-01-01")),
                 parse = FALSE, hjust = 0, size = 4.5, fontface = "bold",
                 family = 'mono', color = 'brown3') +
      scale_y_continuous(name = "Ventilated Patients", sec.axis = sec_axis(~sec$rev(.), name = paste0(gsub("_", " ", depvar), " Count"))) +
      theme_bw() +
      scale_fill_manual()+
      labs(x = "Date") +
      theme(strip.text.x = element_text(size = 18, color = "black", face = "bold",family = "mono"),
            strip.text.y = element_text(size = 18, color = "black", face = "bold",family = "mono"),
            legend.title = element_text(family = "mono",size = 14),
            legend.position="bottom",
            legend.text = element_text(family = "mono",size = 14),
            axis.line.x = element_line(color="black", size = 0.8),
            axis.line.y = element_line(color="black", size = 0.8),
            axis.title.y  = element_text(family = "mono",size = 14, face = "bold", color = color1),
            axis.title.y.right = element_text(family = "mono",size = 14, face = "bold", color = color2),
            axis.title.x  = element_text(family = "mono",size = 14, face = "bold", hjust = 0.5),
            axis.text.x   = element_text(family = "mono",size = 14, face = "bold", angle = 90),
            axis.text.y   = element_text(family = "mono",size = 14, face = "bold"),
            plot.title = element_text(family = "mono",size = 18,face = "bold", hjust = 0.5)) +
      ggtitle(paste0("Hospital ", grp))
    fig_list[[counter]] <- fig
    counter = counter + 1
  }
  fig <- grid.arrange(arrangeGrob(fig_list[[1]]+ theme(legend.position="none"),
                                  fig_list[[2]]+ theme(legend.position="none"),
                                  fig_list[[3]]+ theme(legend.position="none"),
                                  fig_list[[4]]+ theme(legend.position="none"),
                                  fig_list[[5]]+ theme(legend.position="none"),
                                  fig_list[[6]]+ theme(legend.position="none"),
                                  nrow = 2), ncol = 1,
                      top =textGrob(paste0("Correlation between ventilated patient census and ",gsub("_", " ", depvar)),gp=gpar(fontsize=20, font=4)))
  ggsave(paste0(pathwrite,"Vent_vs_",depvar,".jpg"), plot = fig, width = 14, height = 10)
}

# Regression models
AIClist <- list()
modellist <- list()
model1 <- list()
model2 <- list()
model3 <- list()
model <- list()
for (ii in 1:6){
  tmp <- subset(dfReg, Hospital_Code == ii)
  # model 1, only ICU FC
  model1[[ii]] <- glm.nb(Vent2 ~  ICU_FC , data = tmp)
  # model 2, ED FC added
  model2[[ii]] <- glm.nb(Vent2 ~  ICU_FC + ED_FC, data = tmp)
  # model 3, best lag of ventilator (y) is added
  Error_List = list()
  set.seed(3456)
  for (d in 1:54){
    form <- as.formula(Vent2 ~  ED_FC + ICU_FC + c(rep(0,d),diff(Vent2,lag=d)))
    model <- glm.nb(form,data = tmp)
    Error_List[d] <- cv.glm(tmp,model,K=10)$delta[2]
  }
  Error_List <- unlist(Error_List)
  dopt = which(Error_List == min(Error_List))
  model3[[ii]] <- glm.nb(Vent2 ~  ICU_FC + ED_FC + c(rep(0,dopt),diff(Vent2,lag=dopt)), data = tmp)
  AIClist[[ii]] <- c(round(AIC(model1[[ii]]),2), round(AIC(model2[[ii]]),2), round(AIC(model3[[ii]]),2))
  idx <- which(AIClist[[ii]] == min(AIClist[[ii]]))
  if(idx == 1) {
    best_model <- model1[[ii]]
  } else if (idx == 2) {
    best_model <- model2[[ii]]
  } else{
    best_model <- model3[[ii]]
  }
  modellist[[ii]] <- best_model
}
