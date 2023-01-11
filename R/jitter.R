
x <- 1:5

x

y<-exp(1:5) #지수 함수 e^y

y

plot(x,y)

#titles

plot(x, y, main="This is main", sub="This is sub",
     xlab="this is xlab", ylab="this is ylab")


# 3 column, 4 rows plots
par(mfrow=c(3,4))



# plot type(), lwd (선두께), col, bg(안쪽색), cex, pch(점모양)
plot (x,y, type = "p", lwd=2, col="darkseagreen", pch=24, bg="yellow", cex=0.5)
plot (x,y, type = "p", lwd=2, col="blue", pch=24, bg="red", cex=3)
plot (x,y, type = "l", lwd=3, col="blue", pch=24, bg="red", cex=1)
plot (x,y, type = "b", lwd=3, col="red", pch=24, bg="green", cex=2)
plot (x,y, type = "p", lwd=3, col="red", pch=1, bg="green", cex=2)
plot (x,y, type = "p", lwd=3, col="red", pch=21, bg="green", cex=2)


grDevices::colors()

par(mfrow=c(1,1))
#xlim, ylim, title 추가
plot(x,y, ann=FALSE, xlim=c(0, 6), ylim=c(0,200))
title(main="title", col.main="red", font.main=7)
title(xlab="xlab", col.lab="blue", font.main=5)
title(ylab="ylab", col.lab="green", font.main=5)

legend("topleft", legend=c("cell", "is", "tastey"), cex=0.8, 
       col=c("blue","red","green"), pch=21:23, lty=1:3)



#create data frame
games_started <- sample(1:10, 300, TRUE)

points_per_game <- 3*games_started + rnorm(300)

hist(rnorm(300))

data <- as.data.frame(cbind(games_started, points_per_game))


View(data)

head(data, 10)

plot(data$games_started, data$points_per_game)

#pch 21~25 제외 bg 안바뀜
plot(data$games_started, data$points_per_game, pch = 21, col = 'steelblue')
plot(data$games_started, data$points_per_game, pch = 21, col = 'steelblue',  bg="red")

boxplot(data$points_per_game~data$games_started)

plot(jitter(data$games_started, 2), data$points_per_game, pch = 16, col = 'steelblue')
plot(jitter(data$games_started, 20), data$points_per_game, pch = 16, col = 'steelblue')

boxplot(data$points_per_game~data$games_started, outline=FALSE)
stripchart(data$points_per_game~data$games_started, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

boxplot(data$points_per_game~data$games_started)

# install.packages("beeswarm")
library(beeswarm)
boxplot(data$points_per_game~data$games_started, outline=FALSE)
beeswarm(data$points_per_game~data$games_started, method="swarm", pch=16, col="blue", add=T)


library(ggplot2)


basegraph<-ggplot(data, aes(x=as.character(games_started), y=points_per_game) )

basegraph


basegraph+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()



basegraph_altered <- ggplot(data, aes(x=as.factor(games_started), y=points_per_game))

basegraph_altered +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()


basegraph_altered +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour = as.factor(games_started))) +
  theme_bw()


#error bar with barplot
birthwt<-MASS::birthwt

head(birthwt)

bwt<-birthwt$bwt

bwt.mean <- mean(bwt)
bwt.sd.upper <- mean(bwt) + sd(bwt)
bwt.sd.lower <- mean(bwt) - sd(bwt)

b <- barplot(bwt.mean, ylim = c(0, max (bwt)))
arrows(b, bwt.sd.upper, b, bwt.sd.lower, angle=90, code=3)


