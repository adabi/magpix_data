library(tidyverse)

file_path <- '../../Data/Compiled/ALL.csv'

df <- read.csv(file_path)

remove_last_two_chars <- function(x) {
  substr(x, 1, nchar(as.character(x)) - 2 )
}

df <- 
  df %>% 
  filter(
    !grepl("Background", Group) & 
      !grepl("Standard", Group) &
      !(Group %in% c("Control1", "Control2", "DEP Interference", "Ox66 Interference"))
  ) %>% 
  mutate(Group = remove_last_two_chars(Group)) %>%   
  separate(Group, into=c("Condition", "Compound"), sep = " ") %>% 
  mutate(Compound = str_replace(Compound, "Blank", "Control")) %>% 
  mutate(Compound = factor(Compound, levels = c('Control', "DEP", "Ox66", "Mix"))) %>% 
  mutate(Condition = factor(Condition, levels = c("Normoxia", "Hypoxia"))) %>%
  arrange(Condition, Compound) %>% 
  unite(Condition, cells, Condition, sep=" - ")

# Function to remove outliers using Tukey's Fences. Returns FALSE for outliers
isnt_out_tukey <- function(x, k = 1.5, na.rm = TRUE) {
  quar <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- diff(quar)
  
  !((quar[1] - k * iqr <= x) & (x <= quar[2] + k * iqr))
}

replace_outliers <- function(y){
  replace(y, isnt_out_tukey(y), NaN)
}


df_splt_lst <- 
  df %>% 
  split(.$Condition)

df_no_out <- 
  lapply(names(df_splt_lst), function(condit){
    df_splt_lst[[condit]] %>% 
      mutate_if(is.numeric, replace_outliers)
  }) %>% bind_rows()

df <- df_no_out

df_splt_sep <-
  df %>% 
  split(.$Condition)

df_scaled <-
  lapply(names(df_splt_sep), function(condit){
    means_df <-
      df %>% 
      filter(Condition == condit, Compound == "Control") %>% 
      group_by(Condition, Compound) %>% 
      summarise_all(function(x) mean(x, na.rm = TRUE)) %>% 
      ungroup
    
    df_splt_sep[[condit]] %>% 
      mutate_if(is.numeric, funs ({
        mean <- 
          means_df %>%
          select(quo_name(quo(.))) %>% 
          pull
        
        ((. - mean) / sd(., na.rm = TRUE))  + 1
      }))
  }) %>% bind_rows()

df_scaled <- 
  df_scaled %>% 
  select(Condition, Compound, everything())

dftest <- 
  df_scaled %>% 
  separate(Condition, into = c('Cells', 'Condition'), sep=" - ") %>% 
  gather(Analyte, Value, GRO.KC.CINC.1:TNFa) %>% 
  filter(!is.na(Value))

test_wilc <- function(x, y){
  wilcox.test(as.numeric(unlist(x)), as.numeric(unlist(y)))$p.value
}

df_sig <-
  dftest %>% 
  filter(Compound != "Control")

df_sig <-
  df_sig %>% 
  pivot_wider(names_from = Condition, values_from = Value)

df_sig <-
  df_sig %>% 
  group_by(Cells, Compound, Analyte) %>% 
  mutate(p_val = test_wilc(Hypoxia, Normoxia)) %>% 
  select(-Hypoxia, -Normoxia)

dftestmean <-
  dftest %>% 
  group_by(Cells, Condition, Compound, Analyte) %>% 
  summarise(std_err = sd(Value)/sqrt(length(Value)), Value = mean(Value)) %>% 
  filter(Compound != "Control") %>% 
  ungroup()

df_join <-
  right_join(dftestmean, df_sig)

df_join <-
  df_join %>% 
  mutate(p_val = case_when(
    p_val < 0.05 ~ "*",
    TRUE ~ ""
  ))

df_join <-
  df_join %>% 
  group_by(Cells, Compound, Analyte) %>% 
  mutate(pval_x = max(Value + std_err) + 0.2) %>% 
  ungroup()


df_human <-
  df_join %>% 
  filter(Cells %in% c("A549", "SVG")) %>% 
  mutate(Cells = str_replace(Cells, "A549", "Human Lung")) %>% 
  mutate(Cells = str_replace(Cells, "SVG", "Human Astroglial")) %>% 
  mutate(Condition = factor(Condition, levels=c("Normoxia", "Hypoxia"))) %>% 
  arrange(Condition)

ggplot(df_human) +
  geom_bar(aes(x=Analyte, y=Value, fill=Condition), stat='identity', position='dodge', width=0.8) +
  geom_errorbar(aes(x=Analyte, group=Condition, ymin=Value-std_err, ymax=Value+std_err), width=0.2, position=position_dodge(width = 0.8)) +
  geom_hline(yintercept = 1, linetype='dashed') +
  geom_text(aes(x=Analyte, y=pval_x, label=p_val), size = 6) + 
  labs(title = "Cytokines and Chemokins for Human Cells", y="Mean Normalized Fluorescence") +
  theme_bw() +
  facet_wrap(Cells~Compound, scales="free_x")

ggsave("../../Graphs/Human_Cells.png", device="png", type="cairo", width=8, height = 4.5)

df_rat <-
  df_join %>% 
  filter(Cells %in% c("DITNC1", "RLE6TN")) %>% 
  mutate(Cells = str_replace(Cells, "RLE6TN", "Rat Lung")) %>% 
  mutate(Cells = str_replace(Cells, "DITNC1", "Rat Astroglial")) %>% 
  mutate(Condition = factor(Condition, levels=c("Normoxia", "Hypoxia"))) %>% 
  arrange(Condition) %>% 
  mutate(Analyte = factor(Analyte, levels=c('IL.6', 'GRO.KC.CINC.1', 'TNFa'))) %>% 
  arrange(Analyte)

ggplot(df_rat) +
  geom_bar(aes(x=Analyte, y=Value, fill=Condition), stat='identity', position='dodge', width=0.8) +
  geom_errorbar(aes(x=Analyte, group=Condition, ymin=Value-std_err, ymax=Value+std_err), width=0.2, position=position_dodge(width = 0.8)) +
  geom_hline(yintercept = 1, linetype='dashed') +
  geom_text(aes(x=Analyte, y=pval_x, label=p_val), size = 6) + 
  labs(title = "Cytokines and Chemokins for Rat Cells", y="Mean Normalized Fluorescence") +
  theme_bw() +
  facet_wrap(Cells~Compound, scales="free_x")

ggsave("../../Graphs/Rat_Cells.png", device="png", type="cairo", width=8, height = 4.5)
