

out1 # CO2 / DOC model
out # DIC / CO2 / DOC model

Y1 = out1$Y
Y = out$Y

co2_doc_mod = apply(Y1[6,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12

dic_co2_doc_mod = apply(Y[7,1,,],MARGIN = 1,FUN=mean)/data2$epiVol*12


cor(co2_doc_mod, dic_co2_doc_mod)

plot(co2_doc_mod~ dic_co2_doc_mod)
abline(0,1, lty = 2, lwd = 2)
abline(lm(co2_doc_mod~dic_co2_doc_mod) )


AIC(logLik(lm(co2_doc_mod~dic_co2_doc_mod)))

