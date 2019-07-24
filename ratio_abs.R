gain = seq(from = 14, to= 36, by =2)
loss = seq(from = -7, to=-18, by = -1)
losa = abs(loss)

# gamble matrix Gelskov
gains    = rep(gain,each=12)
losss    = rep(loss,times=12)
ratiog   = abs(gains/losss)
ratiog_2 = ratiog^2

# gamble matrix mine
gains   = rep(gain,each=12)
losas   = rep(losa,times=12)
ratioa  = gains/losas
ratioa  = ratioa - mean(ratioa)
ratioa_2 = ratioa^2
