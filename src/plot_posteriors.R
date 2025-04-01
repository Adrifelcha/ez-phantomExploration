plot_estimates <- function(data){

axisCol <- "black"
x_lim <- c(-1.8,1.8)
x_lab <- round(seq(-1.5,1.5,length.out=5),1)
binWidth = 0.1
binStarts <- c(0.9,1.1,1.9,2.1)
binMids <- binStarts + binWidth / 2
DOYrange <- range(data)
histList <- apply(data, 2, function(x, hCol) hist(x, plot = FALSE))
    
means <- apply(data, 2, mean)
CI <- apply(data, 2, quantile, prob=c(0.025,0.975))
shadeCI <- list()
for(i in 1:ncol(CI)){
  keep <- (data[,i] > CI[1,i]) & (data[,i] < CI[2,i])
  nBreaks <- sum(histList[[i]]$breaks >= CI[1,i] & histList[[i]]$breaks <= CI[2,i])
  partial <- hist(data[keep,i], breaks = nBreaks, plot = FALSE)
  shadeCI <- append(shadeCI, list(partial))
}
  
## Plotting
xlim <- c(0.5,2.5)
ylim <- c(-1.2,2)
xlabs <- c("Convexity", "Concavity")
fillCol <- c("#42B6EC", "#0D41D5", "#42B6EC", "#0D41D5")
CIcolor <- rep("#CAC7DA",4)

plot(c(0, 5), DOYrange, type = "n", xlim=xlim, ylim=DOYrange,
      ann = FALSE, axes = FALSE, xaxs = "i", yaxs = "i")
axis(1, 1.5, "Change type", cex.axis = 0.8, col = NA, line=0, f=2)
axis(1, 1:2, xlabs, cex.axis = 0.7, col = NA, line=-1, f=1)
axis(1, 1:2, c("",""), cex.axis = 0.7, tck=-0.04)
y.seq = format(round(seq(-1.2,2,length.out=7),digits = 1), nsmall = 1)
axis(2,ylim, c("",""), tck=-0, line=0)
axis(2,y.seq, rep("", length(y.seq)), tck=-0.04, line=0)
axis(2,y.seq, y.seq, cex.axis=0.6, lwd=0, line=-0.6, las=2)

mtext(side = 2, outer = F, line = 1.3, expression(paste(nu^pred)), cex = 0.8)
box(bty = "L", col = axisCol)

biggestDensity <- max(unlist(lapply(histList, function(h){max(h[[4]])})))
xscale <- binWidth * .9 / biggestDensity

## Plot the histograms
for (i in 1:4) {
  X <- binStarts[i]
  VerticalHist(x = X, xscale = xscale, xwidth = binWidth, 
                hist= histList[[i]], fillCol = "black")
  VerticalHist(x = X, xscale = xscale*0.8, xwidth = binWidth, 
                hist= shadeCI[[i]], fillCol = fillCol[i])
  points(X, means[i], pch=18, cex=0.5)
  lines(x=c(X-0.05, X+0.05), y=c(means[i],means[i]), lwd=1)
}
#lines(binStarts[c(1,3)], means[c(1,3)], lwd=1, lty=2)
#lines(binStarts[c(2,4)], means[c(2,4)], lwd=1, lty=2)
#text(1.1,1.85,"Change quality", cex=0.8, f=2)
#lines(c(0.59,0.64),c(1.55,1.55), lwd=3, col="#42B6EC")
#text(0.975,1.55,"Qualitative", cex=0.7, f=1)
#lines(c(0.59,0.64),c(1.25,1.25), lwd=3, col="#0D41D5")
#text(1.01,1.25,"Quantitative", cex=0.7, f=1)
}

plot_estimates(drift_pred)
plot_estimates(bound_pred)