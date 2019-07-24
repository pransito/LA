## all candidate models for the 51 sample
prefix   = "glmer(accept.reject ~ "
suffix_1 = ", data=data.la, family='binomial',control=glmerControl(optimizer='bobyqa'))"
suffix_2 = ",data = data.la,family = 'binomial',nAGQ = 0,control=glmerControl(check.conv.grad='ignore',check.conv.singular='ignore',check.conv.hess='ignore',optCtrl=list(optimizer = 'nloptwrap',maxfun=500)))"
form = c()
form[1]  = "1 + (1|subject)"
form[2]  = "(Gewinn+Verlust)+ (Gewinn+Verlust|subject)"
form[3]  = "(Gewinn+Verlust)*group+ (Gewinn+Verlust|subject)"
form[4]  = "(Gewinn+Verlust)*(group+Age+Bildungsjahre_ges)+ (Gewinn+Verlust|subject)"
form[5]  = "(Gewinn+Verlust)*(Age+Bildungsjahre_ges)+ (Gewinn+Verlust|subject)"
form[6]  = "(Gewinn+Verlust)*(Age)+ (Gewinn+Verlust|subject)"
form[7]  = "(Gewinn+Verlust)*(Age+group)+ (Gewinn+Verlust|subject)"
form[8]  = "(Gewinn+Verlust+RiskMar)*group + (Gewinn+Verlust+RiskMar|subject)"


cform = list()
for (ii in 1:length(form)) {
  cform [[ii]] = paste0(prefix, form[ii],suffix_2)
}

all_modc = cform
