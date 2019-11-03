## compute all candidate models
all_mod <- list()

all_mod[[1]]      <- 'glmer(accept.reject ~ RiskMar + (RiskMar|subject), data=data.la, family="binomial",
                           control=glmerControl(optimizer="bobyqa"))'
all_mod[[2]]      <- 'glmer(accept.reject ~ (ed.abs + RiskMar) + ((ed.abs + RiskMar)|subject), data=data.la, family="binomial",
                           control=glmerControl(optimizer="bobyqa"))'

all_mod[[3]]      <- 'glmer(accept.reject ~ RiskMar*group + (RiskMar|subject), data=data.la, family="binomial",
                           control=glmerControl(optimizer="bobyqa"))'
all_mod[[4]]      <- 'glmer(accept.reject ~ (ed.abs + RiskMar)*group + ((ed.abs + RiskMar)|subject), data=data.la, family="binomial",
                           control=glmerControl(optimizer="bobyqa"))'

all_mod[[5]]      <- 'glmer(accept.reject ~ ratio + (ratio|subject), data=data.la, family="binomial",
                           control=glmerControl(optimizer="bobyqa"))'
all_mod[[6]]      <- 'glmer(accept.reject ~ (ed.abs + ratio) + ((ed.abs + ratio)|subject), data=data.la, family="binomial",
                           control=glmerControl(optimizer="bobyqa"))'

all_mod[[7]]      <- 'glmer(accept.reject ~ ratio*group + (ratio|subject), data=data.la, family="binomial",
                           control=glmerControl(optimizer="bobyqa"))'
all_mod[[8]]      <- 'glmer(accept.reject ~ (ed.abs + ratio)*group + ((ed.abs + ratio)|subject), data=data.la, family="binomial",
                           control=glmerControl(optimizer="bobyqa"))'

all_mod[[9]]      <- 'glmer(accept.reject ~ diff + (diff|subject), data=data.la, family="binomial",
                           control=glmerControl(optimizer="bobyqa"))'
all_mod[[10]]      <- 'glmer(accept.reject ~ (ed.abs + diff) + ((ed.abs + diff)|subject), data=data.la, family="binomial",
                            control=glmerControl(optimizer="bobyqa"))'

all_mod[[11]]      <- 'glmer(accept.reject ~ diff*group + (diff|subject), data=data.la, family="binomial",
                            control=glmerControl(optimizer="bobyqa"))'
all_mod[[12]]      <- 'glmer(accept.reject ~ (ed.abs + diff)*group + ((ed.abs + diff)|subject), data=data.la, family="binomial",
                            control=glmerControl(optimizer="bobyqa"))'

all_mod[[13]]      <- 'glmer(accept.reject ~ EV + (EV|subject), data=data.la, family="binomial",
                            control=glmerControl(optimizer="bobyqa"))'
all_mod[[14]]      <- 'glmer(accept.reject ~ (ed.abs + EV) + ((ed.abs + EV)|subject), data=data.la, family="binomial",
                            control=glmerControl(optimizer="bobyqa"))'

all_mod[[15]]      <- 'glmer(accept.reject ~ EV*group + (EV|subject), data=data.la, family="binomial",
                            control=glmerControl(optimizer="bobyqa"))'
all_mod[[16]]      <- 'glmer(accept.reject ~ (ed.abs + EV)*group + ((ed.abs + EV)|subject), data=data.la, family="binomial",
                            control=glmerControl(optimizer="bobyqa"))'

# complete LA model
# without group
# without ed
all_mod[[17]]      <- 'glmer(accept.reject ~ (Gewinn + Verlust) +((Gewinn + Verlust)|subject), data=data.la, family="binomial",
                            control=glmerControl(optimizer="bobyqa"))'
# with ed
all_mod[[18]]     <- 'glmer(accept.reject ~ (ed.abs + Gewinn + Verlust) +((ed.abs + Gewinn + Verlust)|subject), data=data.la, family="binomial",
                           control=glmerControl(optimizer="bobyqa"))'

# with group
# without ed
all_mod[[19]]     <- 'glmer(accept.reject ~ (Gewinn + Verlust)*group +((Gewinn + Verlust)|subject), data=data.la, family="binomial",
                           control=glmerControl(optimizer="bobyqa"))'
# with ed
all_mod[[20]]     <- 'glmer(accept.reject ~ (ed.abs + Gewinn + Verlust)*group +((ed.abs + Gewinn + Verlust)|subject), data=data.la, family="binomial",
                           control=glmerControl(optimizer="bobyqa"))'