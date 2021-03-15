library(tidyverse)

Annotate_results_plate <- function(Results_plate, Plate_plan, Primer_efficiency) {
  Annotated_Results_plate <- Results_plate %>%
    inner_join(Plate_plan) %>%
    inner_join(Primer_efficiency)
}

Primer_efficiency <- read_csv("Primer_Efficiency_Japonicus_TE.csv")

Plate_plan_1 <- read_csv("210720_2742_3054_qRT_1_12_act1_fba1_his3_crm1_Plan.csv")

Results_plate_1 <- read_csv("210720_2742_3054_qRT_1_12_act1_fba1_his3_crm1_Result.csv")

Plate_plan_2 <- read_csv("210720_2742_3054_qRT_1_12_Tj1L_Tj1O_Tj2L_Tj2O_Plan.csv")

Results_plate_2 <- read_csv("210720_2742_3054_qRT_1_12_Tj1L_Tj1O_Tj2L_Tj2O_Result.csv")

Plate_plan_3 <- read_csv("210720_2742_3054_qRT_1_12_Tj3L_Tj3O_Tj4L_Tj4O_Plan.csv")

Results_plate_3 <- read_csv("210720_2742_3054_qRT_1_12_Tj3L_Tj3O_Tj4L_Tj4O_Result.csv")

Plate_plan_4 <- read_csv("210720_2742_3054_qRT_1_12_Tj5L_Tj5O_Tj6L_Tj6O_Plan.csv")

Results_plate_4 <- read_csv("210720_2742_3054_qRT_1_12_Tj5L_Tj5O_Tj6L_Tj6O_Result.csv")

Plate_plan_5 <- read_csv("210720_2742_3054_qRT_1_12_Tj7L_Tj7O_Tj8L_Tj8O_Plan.csv")

Results_plate_5 <- read_csv("210720_2742_3054_qRT_1_12_Tj7L_Tj7O_Tj8L_Tj8O_Result.csv")

Plate_plan_6 <- read_csv("210720_2742_3054_qRT_1_12_Tj9L_Tj9O_Tj10L_Tj10O_Plan.csv")

Results_plate_6 <- read_csv("210720_2742_3054_qRT_1_12_Tj9L_Tj9O_Tj10L_Tj10O_Result.csv")


Annotated_results_plate_1 <- Annotate_results_plate(Results_plate_1, Plate_plan_1, Primer_efficiency)

Annotated_results_plate_2 <- Annotate_results_plate(Results_plate_2, Plate_plan_2, Primer_efficiency)

Annotated_results_plate_3 <- Annotate_results_plate(Results_plate_3, Plate_plan_3, Primer_efficiency)

Annotated_results_plate_4 <- Annotate_results_plate(Results_plate_4, Plate_plan_4, Primer_efficiency)

Annotated_results_plate_5 <- Annotate_results_plate(Results_plate_5, Plate_plan_5, Primer_efficiency)

Annotated_results_plate_6 <- Annotate_results_plate(Results_plate_6, Plate_plan_6, Primer_efficiency)

All_annotated_results <- rbind(Annotated_results_plate_1, Annotated_results_plate_2, Annotated_results_plate_3, Annotated_results_plate_4, Annotated_results_plate_5,
                               Annotated_results_plate_6) %>%
  transmute(Strain = Strain, Replicate = Replicate, Cq = na_if(Cq, "-"), RT = RT, Primer = Primer, Efficiency = Efficiency) %>%
  mutate_all(type.convert) %>%
  mutate_if(is.factor, as.character)





Calculate_qRT_Relative_Enrichment <- function(Annotated_results_plate, Control_Gene_Name) {
  Results_plate_calc <- Annotated_results_plate %>% 
    mutate_all(type.convert) %>%
    mutate_if(is.factor, as.character) %>%
    mutate(Corrected_val = 1/(Efficiency^Cq)) %>%
    group_by(Strain, Replicate, Primer, RT) %>%
    summarise(mean_Corrected_val = mean(Corrected_val))
  
  Control_Primers <- Results_plate_calc %>%
    subset(Primer == Control_Gene_Name)
  
  Gene_Primers <- Results_plate_calc %>%
    subset(Primer != Control_Gene_Name)
  
  Mean_Gene_Control <- inner_join(Gene_Primers, Control_Primers, by = c("Strain", "Replicate", "RT"), suffix = c("_gene", "_control")) %>%
    transmute(Strain = Strain, Replicate = Replicate, Primer_gene = Primer_gene, Primer_control = Primer_control, RT = RT, 
              Primer_over_control = mean_Corrected_val_gene/mean_Corrected_val_control) %>%
    group_by(Strain, Primer_gene, RT) %>%
    summarise(mean_Primer_over_control = mean(Primer_over_control), sd_Primer_over_control = sd(Primer_over_control))
}


Primers <- c("Act1", "Fba1", "His3", "Crm1", "Tj1 LTR", "Tj1 ORF", "Tj2 LTR", "Tj2 ORF", "Tj3 LTR", "Tj3 ORF", "Tj4 LTR", 
             "Tj4 ORF", "Tj5 LTR", "Tj5 ORF", "Tj6 LTR", "Tj6 ORF", "Tj7 LTR", "Tj7 ORF", 
             "Tj8 LTR", "Tj8 ORF", "Tj9 LTR", "Tj9 ORF", "Tj10 LTR", "Tj10 ORF")

WT_Dcr1 <- c("wild_type", "dcr1")

qRT_calc <- Calculate_qRT_Relative_Enrichment(All_annotated_results, "His3") %>%
  mutate(Primer_gene = factor(Primer_gene, Primers), Strain = factor(Strain, WT_Dcr1))

qRT_calc_Crm1 <- Calculate_qRT_Relative_Enrichment(All_annotated_results, "Crm1") %>%
  mutate(Primer_gene = factor(Primer_gene, Primers), Strain = factor(Strain, WT_Dcr1))


qRT_calc_Tj1 <- subset(qRT_calc, RT == "Plus") %>% 
  subset(Primer_gene == "Tj1 LTR" | Primer_gene == "Tj1 ORF") %>%
  mutate(Strain = factor(Strain, WT_Dcr1))

qRT_calc_Tj2 <- subset(qRT_calc, RT == "Plus") %>% 
  subset(Primer_gene == "Tj2 LTR" | Primer_gene == "Tj2 ORF") %>%
  mutate(Strain = factor(Strain, WT_Dcr1))

qRT_calc_Tj3 <- subset(qRT_calc, RT == "Plus") %>% 
  subset(Primer_gene == "Tj3 LTR" | Primer_gene == "Tj3 ORF") %>%
  mutate(Strain = factor(Strain, WT_Dcr1))

qRT_calc_Tj4 <- subset(qRT_calc, RT == "Plus") %>% 
  subset(Primer_gene == "Tj4 LTR" | Primer_gene == "Tj4 ORF") %>%
  mutate(Strain = factor(Strain, WT_Dcr1))

qRT_calc_Tj5 <- subset(qRT_calc, RT == "Plus") %>% 
  subset(Primer_gene == "Tj5 LTR" | Primer_gene == "Tj5 ORF") %>%
  mutate(Strain = factor(Strain, WT_Dcr1))

qRT_calc_Tj6 <- subset(qRT_calc, RT == "Plus") %>% 
  subset(Primer_gene == "Tj6 LTR" | Primer_gene == "Tj6 ORF") %>%
  mutate(Strain = factor(Strain, WT_Dcr1))

qRT_calc_Tj7 <- subset(qRT_calc, RT == "Plus") %>% 
  subset(Primer_gene == "Tj7 LTR" | Primer_gene == "Tj7 ORF") %>%
  mutate(Strain = factor(Strain, WT_Dcr1))

qRT_calc_Tj8 <- subset(qRT_calc, RT == "Plus") %>% 
  subset(Primer_gene == "Tj8 LTR" | Primer_gene == "Tj8 ORF") %>%
  mutate(Strain = factor(Strain, WT_Dcr1))

qRT_calc_Tj9 <- subset(qRT_calc, RT == "Plus") %>% 
  subset(Primer_gene == "Tj9 LTR" | Primer_gene == "Tj9 ORF") %>%
  mutate(Strain = factor(Strain, WT_Dcr1))

qRT_calc_Tj10 <- subset(qRT_calc, RT == "Plus") %>% 
  subset(Primer_gene == "Tj10 LTR" | Primer_gene == "Tj10 ORF") %>%
  mutate(Strain = factor(Strain, WT_Dcr1))

qRT_calc_control <- subset(qRT_calc, RT == "Plus") %>% 
  subset(Primer_gene == "Act1" | Primer_gene == "Fba1"| Primer_gene == "His3"| Primer_gene == "Crm1") %>%
  mutate(Strain = factor(Strain, WT_Dcr1))

qRT_calc_control_Crm1 <- subset(qRT_calc_Crm1, RT == "Plus") %>% 
  subset(Primer_gene == "Act1" | Primer_gene == "Fba1"| Primer_gene == "His3"| Primer_gene == "Crm1") %>%
  mutate(Strain = factor(Strain, WT_Dcr1))

qRT_calc_TE <- subset(qRT_calc, RT == "Plus") %>% 
  subset(Primer_gene != "Act1" & Primer_gene != "Fba1"& Primer_gene != "His3"& Primer_gene != "Crm1") %>%
  mutate(Strain = factor(Strain, WT_Dcr1))



ggplot(qRT_calc_control_Crm1, aes(x = Primer_gene, y = mean_Primer_over_control, fill = Strain)) + 
  geom_bar(stat = "identity", width = 0.5, size = 0.75, colour = "black", position=position_dodge(0.6)) +
  geom_errorbar(aes(ymin=mean_Primer_over_control-sd_Primer_over_control, ymax=mean_Primer_over_control+sd_Primer_over_control), 
                width = 0.2, size = 0.75, position=position_dodge(0.6)) + 
  scale_fill_manual(values=c("black","darkgrey"), name = "Strain", labels = c("WT", "Dcr1\u0394")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,16), breaks=seq(0,16,2)) +
  labs(title=("qRT-PCR Dcr1 Control Gene Expression"), x="", y = "Ratio (TE/his3)") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size =16, face = "bold", colour = "black"),
    axis.text = element_text(size = 14, face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 45, hjust =1),
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
    panel.background = element_rect(fill = "white"), 
    axis.ticks = element_line(colour = "black", size = 1),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 16))

qPCR_plot_Tj10 + theme(legend.position = "none")

#To calculate significance
qRT_t_test <- function(Annotated_results_plate, Control_Gene_Name) {
  Results_plate_calc <- Annotated_results_plate %>% 
    mutate_all(type.convert) %>%
    mutate_if(is.factor, as.character) %>%
    mutate(Corrected_val = 1/(Efficiency^Cq)) %>%
    group_by(Strain, Replicate, Primer, RT) %>%
    summarise(mean_Corrected_val = mean(Corrected_val))
  
  Control_Primers <- Results_plate_calc %>%
    subset(Primer == Control_Gene_Name)
  
  Gene_Primers <- Results_plate_calc %>%
    subset(Primer != Control_Gene_Name)
  
  Mean_Gene_Control <- inner_join(Gene_Primers, Control_Primers, by = c("Strain", "Replicate", "RT"), suffix = c("_gene", "_control")) %>%
    transmute(Strain = Strain, Replicate = Replicate, Primer_gene = Primer_gene, Primer_control = Primer_control, RT = RT, 
              Primer_over_control = mean_Corrected_val_gene/mean_Corrected_val_control) %>%
    subset(RT == "Plus")
  
  names(Mean_Gene_Control)[names(Mean_Gene_Control) == "Primer_gene"] <- "Primer"
  
  t_test_calc <- Mean_Gene_Control %>%
    group_by(Strain, Primer) %>%
    summarise(Primer_over_control = list(Primer_over_control))%>%
    spread(Strain, Primer_over_control) %>%
    group_by(Primer) %>% 
    transmute(Primer = Primer, p_value = t.test(unlist(dcr1), unlist(wild_type), var.equal=TRUE)$p.value,
              t_value = t.test(unlist(dcr1), unlist(wild_type), var.equal=TRUE)$statistic)
  
}

t_test_qRT_His3 <- qRT_t_test(All_annotated_results, "His3")

t_test_qRT_Act1 <- qRT_t_test(All_annotated_results, "Act1")

# To calculate with WT set as 1

Calculate_qRT_Over_WT <- function(Annotated_results_plate, Control_Gene_Name) {
  Results_plate_calc <- Annotated_results_plate %>% 
    mutate_all(type.convert) %>%
    mutate_if(is.factor, as.character) %>%
    mutate(Corrected_val = 1/(Efficiency^Cq)) %>%
    group_by(Strain, Replicate, Primer, RT) %>%
    summarise(mean_Corrected_val = mean(Corrected_val))
  
  Control_Primers <- Results_plate_calc %>%
    subset(Primer == Control_Gene_Name)
  
  Gene_Primers <- Results_plate_calc %>%
    subset(Primer != Control_Gene_Name)
  
  Mean_Gene_Control <- inner_join(Gene_Primers, Control_Primers, by = c("Strain", "Replicate", "RT"), suffix = c("_gene", "_control")) %>%
    transmute(Strain = Strain, Replicate = Replicate, Primer_gene = Primer_gene, Primer_control = Primer_control, RT = RT, 
              Primer_over_control = mean_Corrected_val_gene/mean_Corrected_val_control)
  
  Wild_Type_Results <- Mean_Gene_Control %>%
    subset(RT == "Plus" & Strain == "wild_type") %>%
    group_by(Strain, Primer_gene, RT) %>%
    summarise(mean_Primer_over_control = mean(Primer_over_control))
  
  names(Wild_Type_Results)[names(Wild_Type_Results) == "Primer_gene"] <- "Primer"
  
  All_Results <- Mean_Gene_Control %>%
    subset(RT == "Plus")
  
  names(All_Results)[names(All_Results) == "Primer_gene"] <- "Primer"
  
  Over_WT_Results <- inner_join(All_Results, Wild_Type_Results, by = c("Primer", "RT")) %>%
    transmute(Strain.x = Strain.x, Replicate = Replicate, Primer = Primer, 
              Over_WT = Primer_over_control/mean_Primer_over_control) %>%
    group_by(Strain.x, Primer) %>%
    summarise(mean_Over_WT = mean(Over_WT), sd_Over_WT = sd(Over_WT))
  
}

Over_WT_Act1 <- Calculate_qRT_Over_WT(All_annotated_results, "Act1")

Over_WT_His3 <- Calculate_qRT_Over_WT(All_annotated_results, "His3")%>%
  mutate(Strain.x = factor(Strain.x, WT_Dcr1), Primer = factor(Primer, Primers))

Dcr1_only_over_WT <- Over_WT_His3 %>%
  subset(Strain.x == "dcr1") %>%
  mutate(Strain.x = factor(Strain.x, WT_Dcr1), Primer = factor(Primer, Primers))

#Over_WT subsetting

qRT_calc_Tj1_WT <- subset(Over_WT_His3) %>% 
  subset(Primer == "Tj1 LTR" | Primer == "Tj1 ORF") %>%
  mutate(Strain.x = factor(Strain.x, WT_Dcr1))

qRT_calc_Tj2_WT <- subset(Over_WT_His3) %>% 
  subset(Primer == "Tj2 LTR" | Primer == "Tj2 ORF") %>%
  mutate(Strain.x = factor(Strain.x, WT_Dcr1))

qRT_calc_Tj3_WT <- subset(Over_WT_His3) %>% 
  subset(Primer == "Tj3 LTR" | Primer == "Tj3 ORF") %>%
  mutate(Strain.x = factor(Strain.x, WT_Dcr1))

qRT_calc_Tj4_WT <- subset(Over_WT_His3) %>% 
  subset(Primer == "Tj4 LTR" | Primer == "Tj4 ORF") %>%
  mutate(Strain.x = factor(Strain.x, WT_Dcr1))

qRT_calc_Tj5_WT <- subset(Over_WT_His3) %>% 
  subset(Primer == "Tj5 LTR" | Primer == "Tj5 ORF") %>%
  mutate(Strain.x = factor(Strain.x, WT_Dcr1))

qRT_calc_Tj6_WT <- subset(Over_WT_His3) %>% 
  subset(Primer == "Tj6 LTR" | Primer == "Tj6 ORF") %>%
  mutate(Strain.x = factor(Strain.x, WT_Dcr1))

qRT_calc_Tj7_WT <- subset(Over_WT_His3) %>% 
  subset(Primer == "Tj7 LTR" | Primer == "Tj7 ORF") %>%
  mutate(Strain.x = factor(Strain.x, WT_Dcr1))

qRT_calc_Tj8_WT <- subset(Over_WT_His3) %>% 
  subset(Primer == "Tj8 LTR" | Primer == "Tj8 ORF") %>%
  mutate(Strain.x = factor(Strain.x, WT_Dcr1))

qRT_calc_Tj9_WT <- subset(Over_WT_His3) %>% 
  subset(Primer == "Tj9 LTR" | Primer == "Tj9 ORF") %>%
  mutate(Strain.x = factor(Strain.x, WT_Dcr1))

qRT_calc_Tj10_WT <- subset(Over_WT_His3) %>% 
  subset(Primer == "Tj10 LTR" | Primer == "Tj10 ORF") %>%
  mutate(Strain.x = factor(Strain.x, WT_Dcr1))

qRT_calc_control_WT <- subset(Over_WT_His3) %>% 
  subset(Primer == "Act1" | Primer == "Fba1") %>%
  mutate(Strain.x = factor(Strain.x, WT_Dcr1))

ggplot(Over_WT_His3, aes(x = Primer, y = mean_Over_WT, fill = Strain.x)) + 
  geom_bar(stat = "identity", width = 0.5, size = 0.75, colour = "black", position=position_dodge(0.6)) +
  geom_errorbar(aes(ymin=mean_Over_WT-sd_Over_WT, ymax=mean_Over_WT+sd_Over_WT), 
                width = 0.2, size = 0.75, position=position_dodge(0.6)) + 
  scale_fill_manual(values=c('black','darkgray'), name = "Strain", labels = c("WT", "Dcr1\u0394")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1100), breaks=seq(0,1100,100)) +
  labs(title=("Dcr1\u0394 qRT-PCR"), x="", y = "Ratio [(gene/his3)/WT]") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size =16, face = "bold", colour = "black"),
    axis.text = element_text(size = 14, face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 45, hjust =1),
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
    panel.background = element_rect(fill = "white"), 
    axis.ticks = element_line(colour = "black", size = 1),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 16))


qRT_t_test_Over_WT <- function(Annotated_results_plate, Control_Gene_Name) {
  Results_plate_calc <- Annotated_results_plate %>% 
    mutate_all(type.convert) %>%
    mutate_if(is.factor, as.character) %>%
    mutate(Corrected_val = 1/(Efficiency^Cq)) %>%
    group_by(Strain, Replicate, Primer, RT) %>%
    summarise(mean_Corrected_val = mean(Corrected_val))
  
  Control_Primers <- Results_plate_calc %>%
    subset(Primer == Control_Gene_Name)
  
  Gene_Primers <- Results_plate_calc %>%
    subset(Primer != Control_Gene_Name)
  
  Mean_Gene_Control <- inner_join(Gene_Primers, Control_Primers, by = c("Strain", "Replicate", "RT"), suffix = c("_gene", "_control")) %>%
    transmute(Strain = Strain, Replicate = Replicate, Primer_gene = Primer_gene, Primer_control = Primer_control, RT = RT, 
              Primer_over_control = mean_Corrected_val_gene/mean_Corrected_val_control)
  
  Wild_Type_Results <- Mean_Gene_Control %>%
    subset(RT == "Plus" & Strain == "wild_type") %>%
    group_by(Strain, Primer_gene, RT) %>%
    summarise(mean_Primer_over_control = mean(Primer_over_control))
  
  names(Wild_Type_Results)[names(Wild_Type_Results) == "Primer_gene"] <- "Primer"
  
  All_Results <- Mean_Gene_Control %>%
    subset(RT == "Plus")
  
  names(All_Results)[names(All_Results) == "Primer_gene"] <- "Primer"
  
  Over_WT_Results <- inner_join(All_Results, Wild_Type_Results, by = c("Primer", "RT")) %>%
    transmute(Strain.x = Strain.x, Replicate = Replicate, Primer = Primer, 
              Over_WT = Primer_over_control/mean_Primer_over_control) 
  
  names(Over_WT_Results)[names(Over_WT_Results) == "Strain.x"] <- "Strain"
  
  t_test_calc <- Over_WT_Results %>%
    group_by(Strain, Primer) %>%
    summarise(Over_WT = list(Over_WT))%>%
    spread(Strain, Over_WT) %>%
    group_by(Primer) %>% 
    transmute(Primer = Primer, p_value = t.test(unlist(dcr1), unlist(wild_type), var.equal=TRUE)$p.value,
              t_value = t.test(unlist(dcr1), unlist(wild_type), var.equal=TRUE)$statistic)
  
}


########## To look at transcripts in WT situation

qRT_calc_WT <- subset(qRT_calc, RT == "Plus") %>% 
  subset(Strain == "wild_type") %>%
  transmute(Strain = Strain, Primer_gene = Primer, mean_over_control = mean_Primer_over_control, sd_over_control = sd_Primer_over_control)%>%
mutate(Primer_gene = factor(Primer_gene, Primers))

qRT_calc_WT_TE <- qRT_calc_WT %>%
  subset(Primer_gene != "Act1" & Primer_gene != "Fba1" & Primer_gene != "Crm1")


ggplot(qRT_calc_WT_TE , aes(x = Primer_gene, y = mean_over_control)) + 
  geom_bar(stat = "identity", width = 0.5, size = 0.75, colour = "black", position=position_dodge(0.6)) +
  geom_errorbar(aes(ymin=mean_over_control-sd_over_control, ymax=mean_over_control+sd_over_control), 
                width = 0.2, size = 0.75, position=position_dodge(0.6)) + 
  scale_fill_manual(values=c('black','darkgray'), name = "Strain", labels = c("WT", "Dcr1\u0394")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.25), breaks=seq(0,0.25,0.05)) +
  labs(title=("TE WT qRT-PCR"), x="", y = "TE/his3") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size =16, face = "bold", colour = "black"),
    axis.text = element_text(size = 14, face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 45, hjust =1),
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
    panel.background = element_rect(fill = "white"), 
    axis.ticks = element_line(colour = "black", size = 1),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 16))

####### adding a comment for github




######## added comment on github
