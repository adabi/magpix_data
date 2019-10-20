library(tidyverse)
cell_line <- 'SVG'
file_path <- sprintf('../../Data/Compiled/%s.csv', cell_line)
df <- read.csv(file_path)
df <- 
  df %>% 
  filter(
    !grepl("Background", Group) & 
      !grepl("Standard", Group) &
      !(Group %in% c("Control1", "Control2", "DEP Interference", "Ox66 Interference"))
         ) %>% 
  mutate(Group = function(x) substr(x, 1, nchar(x) - 1))

view(count(df, Group))
   