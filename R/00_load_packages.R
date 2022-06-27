# load libraries ----------

#.libPaths(new="/volumes/lab/users/tapsi/R/")
#.libPaths(new="/usr/share/R/library")
#.libPaths(new = c("/usr/lib64/R/library", "/usr/share/R/library", "/volumes/lab/users/tapsi/R/"))
library(pacman)
library(reticulate)
library(Seurat,SeuratWrappers)
library(dplyr)
library(Matrix)
library(pacman)
library(future)
p_load(params, tidyverse, glue, janitor, yarrr)
p_load(RColorBrewer)
library("wesanderson");library("yarrr");library("rcartocolor")


#options(future.globals.maxSize = Inf)
options(scipen=5)
set.seed(21)

# colors ------
colors_neon = c("aquamarine1", "bisque", "blue", "brown", "brown1", "cadetblue", "cyan", "chartreuse3", "chocolate", "coral1", "darkorange", 
                "cornflowerblue", "darkgoldenrod", "darkolivegreen", "darkmagenta", "darkolivegreen1", "deeppink1", "grey", "purple", "darkblue", "black", "pink1", "plum1", "yellow1", "olivedrab1",
                "coral4", "hotpink", "tan1")
colourCount = 14
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
colors_set1 = getPalette(colourCount)
colors_dark = c("blue", "red", "green", "gold2", "navy", "magenta", "darkgreen","deepskyblue",
                "orange", "deeppink", "springgreen4","darkturquoise", "brown3", "purple",
                "yellowgreen", "plum3", "darkgrey", "lavender", "cyan", "forestgreen", "indianred1", "lightblue", "slategray2", "cornsilk3", "black", colors_neon, colors_set1)

# ***Colors --------------------
# http://colorlisa.com/
# https://colorhunt.co/
yel_blue = c("#edc951","#eb6841","#cc2a36","#4f372d","#00a0b0")
botticelli = c("#7A989A", "#849271", "#CF9546", "#C67052", "#C1AE8D")
chagall = c("#3F6F76", "#C65840", "#F4CE4B", "#62496F", "#69B7CE")
delaunay = c("#A4B7E1", "#B8B87A", "#EFBD37", "#A85E5E", "#EFDE80")
ernst = c("#91323A", "#3A4960", "#6D7345", "#554540", "#D7C969")
escher = c("#C1395E", "#AEC17B", "#E07B42", "#89A7C2", "#F0CA50")
gauguin = c("#21344F", "#8AAD05", "#DF5D22", "#E17976", "#E2CE1B")
hopper = c("#67161C", "#3F6148", "#A4804C", "#4B5F80", "#DBD3A4")
jean = c("#51394E", "#C8AF8A", "#658385", "#B04838")
kandinsky = c("#d2981a", "#a53e1f", "#457277", "#8f657d", "#8dcee2")
michelangelo = c("#42819F", "#86AA7D", "#CBB396", "#4D280F", "#555234")
munch = c("#E69253", "#EDB931", "#E4502E", "#4378A0", "#272A2A")

colors_darjee = wes_palette(name = "Darjeeling1", n = 5)
colors_bottle2 = wes_palette(name = "BottleRocket2", n = 5)
colors_bottle2 = wes_palette(name = "BottleRocket2", n = 5)
colors_rcart = carto_pal(7, "Vivid");colors_rcart = colors_rcart[c(1:2,5:6)]
a = wes_palette(name = "Darjeeling1")
b = wes_palette(name = "GrandBudapest1")
c = wes_palette(name = "Royal2")
colors_darj <- c("#FF0000", "#00A08A", "#F2AD00", "plum4", "#F98400", "#5BBCD6", b[1:4], "steelblue4", c[1:5])

