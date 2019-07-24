# de martino
gain = rep(seq(20,50,length.out = 15),times=15)
loss = rep(seq(20,50,length.out = 15),each=15)

gm     = data.frame(gain,loss)
gm$ev  = 0.5*gain + 0.5*loss*(-1)
gm$var = (0.5*gain + 0.5*loss)^2
gm$sd  = sqrt((0.5*gain + 0.5*loss)^2)
plot(gm$ev,gm$var)
plot(gm$gain,gm$var)
plot(gm$ev,gm$sd)
cor(gm$gain,gm$var)
cor(gm$loss,gm$var)
cor(gm$ev,gm$var)
cor(gm$ev,gm$sd)

# our case
gain = rep(seq(14,36,length.out = 12),times=12)
loss = rep(seq(7,18,length.out = 12),each=12)
gmt     = data.frame(gain,loss)
gmt$ev  = 0.5*gain + 0.5*loss*(-1)
gmt$var = (0.5*gain + 0.5*loss)^2
plot(gmt$ev,gmt$var)
cor(gmt$ev,gmt$var)
cor(gmt$gain,gmt$var)
cor(gmt$loss,gmt$var)

# canessa
gain = rep(seq(1,99,length.out = 99),times=99)
loss = rep(seq(1,99,length.out = 99),each=99)

gm     = data.frame(gain,loss)
gm$ev  = 0.5*gain + 0.5*loss*(-1)
gm$var = (0.5*gain + 0.5*loss)^2
plot(gm$ev,gm$var)
plot(gm$gain,gm$var)
cor(gm$gain,gm$var)
cor(gm$ev,gm$var)

# charpentier 1
gain = rep(seq(2,10,length.out = 7),times=7)
loss = rep(seq(2,10,length.out = 7),each=7)

gm     = data.frame(gain,loss)
gm$ev  = 0.5*gain + 0.5*loss*(-1)
gm$var = (0.5*gain + 0.5*loss)^2
plot(gm$ev,gm$var)
plot(gm$gain,gm$var)
cor(gm$gain,gm$var)
cor(gm$ev,gm$var)

# charpentier 2
gain = rep(seq(11,19,length.out = 7),times=7)
loss = rep(seq(2,10,length.out = 7),each=7)

gm     = data.frame(gain,loss)
gm$ev  = 0.5*gain + 0.5*loss*(-1)
gm$var = (0.5*gain + 0.5*loss)^2
plot(gm$ev,gm$var)
plot(gm$gain,gm$var)
cor(gm$gain,gm$var)
cor(gm$ev,gm$var)
