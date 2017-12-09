cleared.txt is 3avalanches1001000000.txt without zeroes.

then, in R, perfrom
g<-c(..copypaste from cleared.txt)

h<-hist(g, breaks=100);
plot(h$breaks[-1], h$counts, log='xy', type='p',xlab='log(avalanche)',ylab='log(frequency)')
x<-100:10000; y=12000000*x^(-1.15)
points(x,y,type="l", col="red")
mtext('y=-1.15x+7',line=-5)

the picture is stored in logloghist.pdf


