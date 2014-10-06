data <- read.csv("multiple_regression.csv")
model <- 'y ~ x1 + x2 + x3 + x4 + x5'
R1 <- rbind(c(0,1,0,0,0,0), c(0,0,1,0,0,0), c(0,0,0,1,0,0),
            c(0,0,0,0,-1,0), c(0,0,0,0,0,-1))

## classical F test ##
fit.lm <- lm(model, data=data)
f <- summary(fit.lm)$fstatistic
# F statistic
f[1]
# p value
pf(f[1],f[2],f[3], lower.tail = FALSE)


## F-bar test ##
fit.csi <- csi.lm(model, data, ui = R1, mix.weights = "boot")
# F-bar statistics
fit.csi$T.obs
# p values
fit.csi$p.value

