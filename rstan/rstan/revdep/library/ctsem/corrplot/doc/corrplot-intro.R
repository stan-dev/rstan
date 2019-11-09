## ----setup, include=FALSE------------------------------------------------
set.seed(0) # we need reproducible results
knitr::opts_chunk$set(
  out.extra = 'style="display:block; margin: auto"',
  fig.align = "center",
  fig.path = "webimg/",
  dev = "png")

## ----methods-------------------------------------------------------------
library(corrplot)
M <- cor(mtcars)
corrplot(M, method = "circle")
corrplot(M, method = "square")
corrplot(M, method = "ellipse")
corrplot(M, method = "number") # Display the correlation coefficient
corrplot(M, method = "shade")
corrplot(M, method = "color")
corrplot(M, method = "pie")

## ----layout--------------------------------------------------------------
corrplot(M, type = "upper")
corrplot(M, type = "upper")

## ----mixed---------------------------------------------------------------
corrplot.mixed(M)
corrplot.mixed(M, lower.col = "black", number.cex = .7)
corrplot.mixed(M, lower = "ellipse", upper = "circle")
corrplot.mixed(M, lower = "square", upper = "circle", tl.col = "black")

## ----order---------------------------------------------------------------
corrplot(M, order = "AOE")
corrplot(M, order = "hclust")
corrplot(M, order = "FPC")
corrplot(M, order = "alphabet")

## ----rectangles----------------------------------------------------------
corrplot(M, order = "hclust", addrect = 2)
corrplot(M, order = "hclust", addrect = 3)

## ----hclust-lightblue----------------------------------------------------
# Change background color to lightblue
corrplot(M, type = "upper", order = "hclust",
         col = c("black", "white"), bg = "lightblue")

## ----color---------------------------------------------------------------
col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("red", "white", "blue"))	
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                           "cyan", "#007FFF", "blue", "#00007F"))
whiteblack <- c("white", "black")

## using these color spectra
corrplot(M, order = "hclust", addrect = 2, col = col1(100))
corrplot(M, order = "hclust", addrect = 2, col = col2(50))
corrplot(M, order = "hclust", addrect = 2, col = col3(20))
corrplot(M, order = "hclust", addrect = 2, col = col4(10))
corrplot(M, order = "hclust", addrect = 2, col = whiteblack, bg = "gold2")

## ----hclust-stdcolors----------------------------------------------------
corrplot(M, order = "hclust", addrect = 2, col = heat.colors(100))
corrplot(M, order = "hclust", addrect = 2, col = terrain.colors(100))
corrplot(M, order = "hclust", addrect = 2, col = cm.colors(100))
corrplot(M, order = "hclust", addrect = 2, col = gray.colors(100))

## ----hclust-rcolorbrewer-------------------------------------------------
library(RColorBrewer)

corrplot(M, type = "upper", order = "hclust",
         col = brewer.pal(n = 8, name = "RdBu"))
corrplot(M, type = "upper", order = "hclust",
         col = brewer.pal(n = 8, name = "RdYlBu"))
corrplot(M, type = "upper", order = "hclust",
         col = brewer.pal(n = 8, name = "PuOr"))

## ----color-label---------------------------------------------------------
## remove color legend and text legend 
corrplot(M, order = "AOE", cl.pos = "n", tl.pos = "n")  

## bottom  color legend, diagonal text legend, rotate text label
corrplot(M, order = "AOE", cl.pos = "b", tl.pos = "d", tl.srt = 60)

## a wider color legend with numbers right aligned
corrplot(M, order = "AOE", cl.ratio = 0.2, cl.align = "r")

## text labels rotated 45 degrees
corrplot(M, type = "lower", order = "hclust", tl.col = "black", tl.srt = 45)

## ----non-corr------------------------------------------------------------
corrplot(abs(M),order = "AOE", col = col3(200), cl.lim = c(0, 1))
## visualize a  matrix in [-100, 100]
ran <- round(matrix(runif(225, -100,100), 15))
corrplot(ran, is.corr = FALSE, method = "square")
## a beautiful color legend 
corrplot(ran, is.corr = FALSE, method = "ellipse", cl.lim = c(-100, 100))

## ----non-corr-asp--------------------------------------------------------
ran <- matrix(rnorm(70), ncol = 7)
corrplot(ran, is.corr = FALSE, win.asp = .7, method = "circle")

## ----NAs-----------------------------------------------------------------
M2 <- M
diag(M2) = NA
corrplot(M2)
corrplot(M2, na.label = "o")
corrplot(M2, na.label = "NA")

## ----plotmath------------------------------------------------------------
M2 <- M[1:5,1:5]
colnames(M2) <- c("alpha", "beta", ":alpha+beta", ":a[0]", "=a[beta]")
rownames(M2) <- c("alpha", "beta", NA, "$a[0]", "$ a[beta]")
corrplot(M2)

## ----test----------------------------------------------------------------
res1 <- cor.mtest(mtcars, conf.level = .95)
res2 <- cor.mtest(mtcars, conf.level = .99)

## specialized the insignificant value according to the significant level
corrplot(M, p.mat = res1$p, sig.level = .2)
corrplot(M, p.mat = res1$p, sig.level = .05)
corrplot(M, p.mat = res1$p, sig.level = .01)

## leave blank on no significant coefficient
corrplot(M, p.mat = res1$p, insig = "blank")

## add p-values on no significant coefficient
corrplot(M, p.mat = res1$p, insig = "p-value")

## add all p-values
corrplot(M, p.mat = res1$p, insig = "p-value", sig.level = -1)

## add cross on no significant coefficient 
corrplot(M, p.mat = res1$p, order = "hclust", insig = "pch", addrect = 3)

## ----ci------------------------------------------------------------------
corrplot(M, low = res1$lowCI, upp = res1$uppCI, order = "hclust",
         rect.col = "navy", plotC = "rect", cl.pos = "n")
corrplot(M, p.mat = res1$p, low = res1$lowCI, upp = res1$uppCI,
         order = "hclust", pch.col = "red", sig.level = 0.01,
         addrect = 3, rect.col = "navy", plotC = "rect", cl.pos = "n")

## ----ci_with_label-------------------------------------------------------
res1 <- cor.mtest(mtcars, conf.level = .95)

corrplot(M, p.mat = res1$p, insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white")
corrplot(M, p.mat = res1$p, method = "color",
         insig = "label_sig", pch.col = "white")
corrplot(M, p.mat = res1$p, method = "color", type = "upper",
         sig.level = c(.001, .01, .05), pch.cex = .9,
         insig = "label_sig", pch.col = "white", order = "AOE")
corrplot(M, p.mat = res1$p, insig = "label_sig", pch.col = "white",
         pch = "p<.05", pch.cex = .5, order = "AOE")

## ----pmat----------------------------------------------------------------
# matrix of the p-value of the correlation
p.mat <- cor.mtest(mtcars)$p
head(p.mat[, 1:5])

# Specialized the insignificant value according to the significant level
corrplot(M, type = "upper", order = "hclust", 
         p.mat = p.mat, sig.level = 0.01)

# Leave blank on no significant coefficient
corrplot(M, type = "upper", order = "hclust", 
         p.mat = p.mat, sig.level = 0.01, insig = "blank")

## ----customized----------------------------------------------------------
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method = "color", col = col(200),
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)

## ----large_matrix--------------------------------------------------------

# generating large feature matrix (cols=features, rows=samples)
num_features <- 60 # how many features
num_samples <- 300 # how many samples
DATASET <- matrix(runif(num_features * num_samples),
               nrow = num_samples, ncol = num_features)

# setting some dummy names for the features e.g. f23
colnames(DATASET) <- paste0("f", 1:ncol(DATASET))

# let's make 30% of all features to be correlated with feature "f1"
num_feat_corr <- num_features * .3
idx_correlated_features <- as.integer(seq(from = 1,
                                          to = num_features,
                                          length.out = num_feat_corr))[-1]
for (i in idx_correlated_features) {
  DATASET[,i] <- DATASET[,1] + runif(num_samples) # adding some noise
}

corrplot(cor(DATASET), diag = FALSE, order = "FPC",
         tl.pos = "td", tl.cex = 0.5, method = "color", type = "upper")

