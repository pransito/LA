## with cov
all_modc <- list()

# RiskMar
all_modc[[1]]      <- "glmer(accept.reject ~ RiskMar + (1+Age) + (RiskMar|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"
all_modc[[2]]      <- "glmer(accept.reject ~ (ed.abs + RiskMar) + (1+Age) + ((ed.abs + RiskMar)|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"

all_modc[[3]]      <- "glmer(accept.reject ~ RiskMar + (1+Bildungsjahre_ges) + (RiskMar|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"
all_modc[[4]]      <- "glmer(accept.reject ~ (ed.abs + RiskMar) + (1+Bildungsjahre_ges) + ((ed.abs + RiskMar)|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"

all_modc[[5]]      <- "glmer(accept.reject ~ RiskMar + (1+Age+Bildungsjahre_ges) + (RiskMar|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"
all_modc[[6]]      <- "glmer(accept.reject ~ (ed.abs + RiskMar) + (1+Age+Bildungsjahre_ges) + ((ed.abs + RiskMar)|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"


all_modc[[7]]      <- "glmer(accept.reject ~ RiskMar*group +(1+Age) + (RiskMar|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"
all_modc[[8]]      <- "glmer(accept.reject ~ (ed.abs + RiskMar)*group +(1+Age) + ((ed.abs + RiskMar)|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"

all_modc[[9]]      <- "glmer(accept.reject ~ RiskMar*group +(1+Bildungsjahre_ges) + (RiskMar|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"
all_modc[[10]]      <- "glmer(accept.reject ~ (ed.abs + RiskMar)*group+(1+Bildungsjahre_ges) + ((ed.abs + RiskMar)|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"

all_modc[[11]]      <- "glmer(accept.reject ~ RiskMar*group + (1+Age+Bildungsjahre_ges) + (RiskMar|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"
all_modc[[12]]      <- "glmer(accept.reject ~ (ed.abs + RiskMar)*group + (1+Age+Bildungsjahre_ges) + ((ed.abs + RiskMar)|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"

# ratio
all_modc[[13]]      <- "glmer(accept.reject ~ ratio +(1+Age) + (ratio|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"
all_modc[[14]]      <- "glmer(accept.reject ~ (ed.abs + ratio) +(1+Age) + ((ed.abs + ratio)|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"

all_modc[[15]]      <- "glmer(accept.reject ~ ratio +(1+Bildungsjahre_ges) + (ratio|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"
all_modc[[16]]      <- "glmer(accept.reject ~ (ed.abs + ratio) + (1+Bildungsjahre_ges) + ((ed.abs + ratio)|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"

all_modc[[17]]      <- "glmer(accept.reject ~ ratio + (1+Age+Bildungsjahre_ges)  + (ratio|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"
all_modc[[18]]      <- "glmer(accept.reject ~ (ed.abs + ratio) +(1+Age+Bildungsjahre_ges) + ((ed.abs + ratio)|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"


all_modc[[19]]      <- "glmer(accept.reject ~ ratio*group +(1+Age) + (ratio|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"
all_modc[[20]]      <- "glmer(accept.reject ~ (ed.abs + ratio)*group +(1+Age) + ((ed.abs + ratio)|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"

all_modc[[21]]      <- "glmer(accept.reject ~ ratio*group +(1+Bildungsjahre_ges) + (ratio|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"
all_modc[[22]]      <- "glmer(accept.reject ~ (ed.abs + ratio)*group +(1+Bildungsjahre_ges) + ((ed.abs + ratio)|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"

all_modc[[23]]      <- "glmer(accept.reject ~ ratio*group +(1+Age+Bildungsjahre_ges) + (ratio|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"
all_modc[[24]]      <- "glmer(accept.reject ~ (ed.abs + ratio)*group +(1+Age+Bildungsjahre_ges) + ((ed.abs + ratio)|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"

# diff
all_modc[[25]]      <- "glmer(accept.reject ~ diff +(1+Age) + (diff|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"
all_modc[[26]]      <- "glmer(accept.reject ~ (ed.abs + diff) +(1+Age) + ((ed.abs + diff)|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"

all_modc[[27]]      <- "glmer(accept.reject ~ diff +(1+Bildungsjahre_ges) + (diff|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"
all_modc[[28]]      <- "glmer(accept.reject ~ (ed.abs + diff) +(1+Bildungsjahre_ges) + ((ed.abs + diff)|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"

all_modc[[29]]      <- "glmer(accept.reject ~ diff +(1+Age+Bildungsjahre_ges) + (diff|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"
all_modc[[30]]      <- "glmer(accept.reject ~ (ed.abs + diff) +(1+Age+Bildungsjahre_ges) + ((ed.abs + diff)|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"


all_modc[[31]]      <- "glmer(accept.reject ~ diff*group +(1+Age) + (diff|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"
all_modc[[32]]      <- "glmer(accept.reject ~ (ed.abs + diff)*group +(1+Age) + ((ed.abs + diff)|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"

all_modc[[33]]      <- "glmer(accept.reject ~ diff*group+(1+Bildungsjahre_ges)  + (diff|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"
all_modc[[34]]      <- "glmer(accept.reject ~ (ed.abs + diff)*group+(1+Bildungsjahre_ges) + ((ed.abs + diff)|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"

all_modc[[35]]      <- "glmer(accept.reject ~ diff*group +(1+Age+Bildungsjahre_ges) + (diff|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"
all_modc[[36]]      <- "glmer(accept.reject ~ (ed.abs + diff)*group +(1+Age+Bildungsjahre_ges) + ((ed.abs + diff)|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"

# EV
all_modc[[37]]      <- "glmer(accept.reject ~ EV + (1+Age) + (EV|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"
all_modc[[38]]      <- "glmer(accept.reject ~ (ed.abs + EV) + (1+Age) + ((ed.abs + EV)|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"

all_modc[[39]]      <- "glmer(accept.reject ~ EV + (1+Bildungsjahre_ges)  + (EV|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"
all_modc[[40]]      <- "glmer(accept.reject ~ (ed.abs + EV) +(1+Bildungsjahre_ges) + ((ed.abs + EV)|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"

all_modc[[41]]      <- "glmer(accept.reject ~ EV +(1+Age+Bildungsjahre_ges) + (EV|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"
all_modc[[42]]      <- "glmer(accept.reject ~ (ed.abs + EV) +(1+Age+Bildungsjahre_ges) + ((ed.abs + EV)|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"


all_modc[[43]]      <- "glmer(accept.reject ~ EV*group + (1+Age) + (EV|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"
all_modc[[44]]      <- "glmer(accept.reject ~ (ed.abs + EV)*group + (1+Age) + ((ed.abs + EV)|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"

all_modc[[45]]      <- "glmer(accept.reject ~ EV*group + (1+Bildungsjahre_ges) + (EV|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"
all_modc[[46]]      <- "glmer(accept.reject ~ (ed.abs + EV)*group + (1+Bildungsjahre_ges) + ((ed.abs + EV)|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"

all_modc[[47]]      <- "glmer(accept.reject ~ EV*group +(1+Age+Bildungsjahre_ges) + (EV|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"
all_modc[[48]]      <- "glmer(accept.reject ~ (ed.abs + EV)*group +(1+Age+Bildungsjahre_ges) + ((ed.abs + EV)|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"

# complete LA model
# without group
# without ed
all_modc[[49]]      <- "glmer(accept.reject ~ (Gewinn + Verlust) +(1+Age) +((Gewinn + Verlust)|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"
all_modc[[50]]      <- "glmer(accept.reject ~ (Gewinn + Verlust) +(1+Bildungsjahre_ges) +((Gewinn + Verlust)|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"
all_modc[[51]]      <- "glmer(accept.reject ~ (Gewinn + Verlust) +(1+Age+Bildungsjahre_ges) +((Gewinn + Verlust)|subject), data=data.la, family='binomial',
                             control=glmerControl(optimizer='bobyqa'))"

# with ed
all_modc[[52]]     <- "glmer(accept.reject ~ (ed.abs + Gewinn + Verlust) +(1+Age) +((ed.abs + Gewinn + Verlust)|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"
all_modc[[53]]     <- "glmer(accept.reject ~ (ed.abs + Gewinn + Verlust) +(1+Bildungsjahre_ges) +((ed.abs + Gewinn + Verlust)|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"
all_modc[[54]]     <- "glmer(accept.reject ~ (ed.abs + Gewinn + Verlust) +(1+Age+Bildungsjahre_ges) +((ed.abs + Gewinn + Verlust)|subject), data=data.la, family='binomial',
                            control=glmerControl(optimizer='bobyqa'))"

# with group
# without ed
all_modc[[55]]     <- "glmer(accept.reject ~ (Gewinn + Verlust)*group +(1+Age)+((Gewinn + Verlust)|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"
all_modc[[56]]     <- "glmer(accept.reject ~ (Gewinn + Verlust)*group +(1+Bildungsjahre_ges)+((Gewinn + Verlust)|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"
all_modc[[57]]     <- "glmer(accept.reject ~ (Gewinn + Verlust)*group +(1+Age+Bildungsjahre_ges)+((Gewinn + Verlust)|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"
# with ed
all_modc[[58]]     <- "glmer(accept.reject ~ (ed.abs + Gewinn + Verlust)*group+(1+Age) +((ed.abs + Gewinn + Verlust)|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"
all_modc[[59]]     <- "glmer(accept.reject ~ (ed.abs + Gewinn + Verlust)*group+(1+Bildungsjahre_ges)+((ed.abs + Gewinn + Verlust)|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"
all_modc[[60]]     <- "glmer(accept.reject ~ (ed.abs + Gewinn + Verlust)*group+(1+Age+Bildungsjahre_ges) +((ed.abs + Gewinn + Verlust)|subject), data=data.la, family='binomial',
                           control=glmerControl(optimizer='bobyqa'))"



