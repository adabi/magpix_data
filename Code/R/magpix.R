library(tidyverse)
cell_line <- 'A549'
file_path <- sprintf('../../Data/Compiled/%s.csv', cell_line)
df <- read.csv(file_path)

remove_last_two_chars <- function(x) {
  substr(x, 1, nchar(as.character(x)) - 2 )
}

df <- 
  df %>% 
  select(-X) %>% 
  filter(
    !grepl("Background", Group) & 
      !grepl("Standard", Group) &
      !(Group %in% c("Control1", "Control2", "DEP Interference", "Ox66 Interference"))
         ) %>% 
  mutate(Group = remove_last_two_chars(Group)) %>%   
  separate(Group, into=c("Condition", "Compound"), sep = " ") %>% 
  mutate(Compound = factor(Compound, levels = c('Control', "DEP", "Ox66", "Mix"))) %>% 
  mutate(Condition = factor(Condition, levels = c("Normoxia", "Hypoxia"))) %>%
  arrange(Condition, Compound)

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

std_err <- function(x){
  x <- na.omit(x)
  sd(x)/sqrt(length(x))
}


df_means <-
  df_scaled %>% 
  group_by(Condition, Compound) %>% 
  summarise_all(list(~mean(.,na.rm=TRUE), ~std_err(.)))

analytes <- names(df)[-(1:2)]

for (analyte in analytes){
  analyte_mean <- sprintf('%s_mean', analyte)
  analyte_err <- sprintf('%s_std_err', analyte)

  plt <-
    ggplot(data=df_means) +
    geom_bar(mapping = aes_string(x='Condition', y=analyte_mean, fill='Compound'), stat='identity', position = 'dodge', width=0.8) +
    geom_errorbar(aes(x=Condition, 
                      ymax = !!as.name(analyte_mean)+!!as.name(analyte_err),
                      ymin=!!as.name(analyte_mean)-!!as.name(analyte_err)), 
                  width = 0.2,
                  position = position_dodge(0.8)) +
    theme_bw(base_size = 12)
  
  print(plt)
}  

