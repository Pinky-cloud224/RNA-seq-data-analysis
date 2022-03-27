library(tidyverse)
go <- read.csv("~/Downloads/gProfiler_scerevisiae_1-29-2022_4-44-50-PM__intersections.csv")
go <- go %>% select("term_name", "adjusted_p_value", "term_size", "intersection_size")
go <- go %>% mutate(Gene_ratio=intersection_size/term_size)
# mutate function is used to create a new variable from the given data file
ggplot(go, aes(x=Gene_ratio, y=term_name, size=intersection_size, color=adjusted_p_value)) + geom_point()
# geom_point() is used to add layer of points in the plot