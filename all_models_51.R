## all candidate models for the 51 sample
prefix   = "glmer(accept.reject ~ "
suffix_1 = ", data=data.la, family='binomial',control=glmerControl(optimizer='bobyqa'))"
suffix_2 = ",data = data.la,family = 'binomial',nAGQ = 0,control=glmerControl(check.conv.grad='ignore',check.conv.singular='ignore',check.conv.hess='ignore',optCtrl=list(optimizer = 'nloptwrap',maxfun=500)))"
form = c()
form[1]  = "1 + (1|subject)"
form[2]  = "(Gewinn+Verlust+ed.abs)+ (Gewinn+Verlust+ed.abs|subject)"
form[3]  = "(Gewinn+Verlust+ed.abs)*group+ (Gewinn+Verlust+ed.abs|subject)"
form[4]  = "(Gewinn+Verlust+ed.abs)*(group+Age+Bildungsjahre_ges)+ (Gewinn+Verlust+ed.abs|subject)"
form[5]  = "(Gewinn+Verlust+ed.abs)*(Age+Bildungsjahre_ges)+ (Gewinn+Verlust+ed.abs|subject)"
form[6]  = "(Gewinn+Verlust+ed.abs)*(Age)+ (Gewinn+Verlust+ed.abs|subject)"
form[7]  = "(Gewinn+Verlust+ed.abs)*(Age+group)+ (Gewinn+Verlust+ed.abs|subject)"


cform = list()
for (ii in 1:length(form)) {
  if ((ii != 4) & (ii != 5) & (ii != 6) & (ii != 7)) {
    cform [[ii]] = paste0(prefix, form[ii],suffix_1)
  } else {
    cform [[ii]] = paste0(prefix, form[ii],suffix_2)
  }
  
}

all_modc = cform
