
#now Ct, make sure that the negatives don't get counted
tot.merge$GhV_ct[tot.merge$GhV_ct==0] <- NA
tot.merge$AngV_ct[tot.merge$AngV_ct==0] <- NA

#GhV
gamGhVCt <- gam(GhV_ct ~ 
                  s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                  s(year, bs="re"),
                method = "REML", 
                gamma = 1,
                select = TRUE,
                data= subset(tot.merge, bat_species =="Eidolon dupreanum"))

summary(gamGhVCt) #year and roost site Moramanga and Ankarana
plot(gamGhVCt, pages = 1)

gamAngVCt <- gam(AngV_ct ~ 
                   s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                   s(year, bs="re"),
                 method = "REML", 
                 gamma = 1,
                 select = TRUE,
                 data= subset(tot.merge, bat_species =="Eidolon dupreanum" & AngV_ct>0))
summary(gamAngVCt) #nothing

#rousettus

gamGhVCtRou <- gam(GhV_ct ~ 
                     s(as.numeric(doy), bs="cc", k=10, by=roost_site) +
                     s(year, bs="re"),
                   method = "REML", 
                   gamma = 1,
                   select = TRUE,
                   data= subset(tot.merge, bat_species =="Rousettus madagascariensis" & GhV_ct>0))
summary(gamGhVCtRou) #year and roost site Moramanga & Makira


gamAngVCtRou <- gam(AngV_ct ~ 
                      s(as.numeric(doy), bs="cc", k=9, by=roost_site) +
                      s(year, bs="re"),
                    method = "REML", 
                    gamma = 1,
                    select = TRUE,
                    data= subset(tot.merge, bat_species =="Rousettus madagascariensis" & AngV_ct>0))
summary(gamAngVCtRou) #year


#Ct predictions
pred.ct <- cbind.data.frame(doy=rep(1:365, 2), bat_species= rep(c("Eidolon dupreanum", "Rousettus madagascariensis"), each=365))
pred.ct$year <- 2019
pred.ct$roost_site <- "Moramanga"
pred.ct$antigen <- "GhV"
pred.ct2 <- pred.ct
pred.ct2$antigen <- "AngV"

pred.ct <- rbind(pred.ct, pred.ct2)

pred.ct$ct <- NA
pred.ct$ct_lci <- NA
pred.ct$ct_uci <- NA

sapply(gamGhVCt$smooth, function(x) x$label)

pred.ct$ct[pred.ct$bat_species=="Eidolon dupreanum" & pred.ct$antigen=="GhV"] <- predict(gamGhVCt, type = "response", newdata = pred.ct[pred.ct$bat_species=="Eidolon dupreanum"  & pred.ct$antigen=="GhV",], exclude = c("s(year)", "s(roost_site)"))
pred.ct$ct_lci[pred.ct$bat_species=="Eidolon dupreanum" & pred.ct$antigen=="GhV"] <- predict(gamGhVCt, type = "response", newdata = pred.ct[pred.ct$bat_species=="Eidolon dupreanum"  & pred.ct$antigen=="GhV",], exclude = c("s(year)", "s(roost_site)")) - 1.96*predict(gamGhVCt, type = "response", newdata = pred.ct[pred.ct$bat_species=="Eidolon dupreanum"  & pred.ct$antigen=="GhV",], exclude = c("s(year)"), se.fit = T)$se
pred.ct$ct_uci[pred.ct$bat_species=="Eidolon dupreanum" & pred.ct$antigen=="GhV"] <- predict(gamGhVCt, type = "response", newdata = pred.ct[pred.ct$bat_species=="Eidolon dupreanum"  & pred.ct$antigen=="GhV",], exclude = c("s(year)", "s(roost_site)")) + 1.96*predict(gamGhVCt, type = "response", newdata = pred.ct[pred.ct$bat_species=="Eidolon dupreanum"  & pred.ct$antigen=="GhV",], exclude = c("s(year)"), se.fit = T)$se


pred.ct$ct[pred.ct$bat_species=="Eidolon dupreanum" & pred.ct$antigen=="AngV"] <- predict(gamAngVCt, type = "response", newdata = pred.ct[pred.ct$bat_species=="Eidolon dupreanum"  & pred.ct$antigen=="AngV",], exclude = c("s(year)"))
pred.ct$ct_lci[pred.ct$bat_species=="Eidolon dupreanum" & pred.ct$antigen=="AngV"] <- predict(gamAngVCt, type = "response", newdata = pred.ct[pred.ct$bat_species=="Eidolon dupreanum"  & pred.ct$antigen=="AngV",], exclude = c("s(year)")) - 1.96*predict(gamAngVCt, type = "response", newdata = pred.ct[pred.ct$bat_species=="Eidolon dupreanum"  & pred.ct$antigen=="AngV",], exclude = c("s(year)"), se.fit = T)$se
pred.ct$ct_uci[pred.ct$bat_species=="Eidolon dupreanum" & pred.ct$antigen=="AngV"] <- predict(gamAngVCt, type = "response", newdata = pred.ct[pred.ct$bat_species=="Eidolon dupreanum"  & pred.ct$antigen=="AngV",], exclude = c("s(year)")) + 1.96*predict(gamAngVCt, type = "response", newdata = pred.ct[pred.ct$bat_species=="Eidolon dupreanum"  & pred.ct$antigen=="AngV",], exclude = c("s(year)"), se.fit = T)$se


pred.ct$ct[pred.ct$bat_species=="Rousettus madagascariensis" & pred.ct$antigen=="GhV"] <- predict(gamGhVCtRou, type = "response", newdata = pred.ct[pred.ct$bat_species=="Rousettus madagascariensis"  & pred.ct$antigen=="GhV",], exclude = c("s(year)"))
pred.ct$ct_lci[pred.ct$bat_species=="Rousettus madagascariensis" & pred.ct$antigen=="GhV"] <- predict(gamGhVCtRou, type = "response", newdata = pred.ct[pred.ct$bat_species=="Rousettus madagascariensis"  & pred.ct$antigen=="GhV",], exclude = c("s(year)")) - 1.96*predict(gamGhVCtRou, type = "response", newdata = pred.ct[pred.ct$bat_species=="Rousettus madagascariensis"  & pred.ct$antigen=="GhV",], exclude = c("s(year)"), se.fit = T)$se
pred.ct$ct_uci[pred.ct$bat_species=="Rousettus madagascariensis" & pred.ct$antigen=="GhV"] <- predict(gamGhVCtRou, type = "response", newdata = pred.ct[pred.ct$bat_species=="Rousettus madagascariensis"  & pred.ct$antigen=="GhV",], exclude = c("s(year)")) + 1.96*predict(gamGhVCtRou, type = "response", newdata = pred.ct[pred.ct$bat_species=="Rousettus madagascariensis"  & pred.ct$antigen=="GhV",], exclude = c("s(year)"), se.fit = T)$se


pred.ct$ct[pred.ct$bat_species=="Rousettus madagascariensis" & pred.ct$antigen=="AngV"] <- predict(gamAngVCtRou, type = "response", newdata = pred.ct[pred.ct$bat_species=="Rousettus madagascariensis"  & pred.ct$antigen=="AngV",], exclude = c("s(year)"))
pred.ct$ct_lci[pred.ct$bat_species=="Rousettus madagascariensis" & pred.ct$antigen=="AngV"] <- predict(gamAngVCtRou, type = "response", newdata = pred.ct[pred.ct$bat_species=="Rousettus madagascariensis"  & pred.ct$antigen=="AngV",], exclude = c("s(year)")) - 1.96*predict(gamAngVCtRou, type = "response", newdata = pred.ct[pred.ct$bat_species=="Rousettus madagascariensis"  & pred.ct$antigen=="AngV",], exclude = c("s(year)"), se.fit = T)$se
pred.ct$ct_uci[pred.ct$bat_species=="Rousettus madagascariensis" & pred.ct$antigen=="AngV"] <- predict(gamAngVCtRou, type = "response", newdata = pred.ct[pred.ct$bat_species=="Rousettus madagascariensis"  & pred.ct$antigen=="AngV",], exclude = c("s(year)")) + 1.96*predict(gamAngVCtRou, type = "response", newdata = pred.ct[pred.ct$bat_species=="Rousettus madagascariensis"  & pred.ct$antigen=="AngV",], exclude = c("s(year)"), se.fit = T)$se


tot.merge$month_doy <- tot.merge$month/12*365-30
tot.merge$month_doy <- yday(tot.merge$mid_date)
tot.merge$month <- month(tot.merge$mid_date) +.5

pred.ct$month <- (pred.ct$doy/365*12) + 1

#and plot
pCtSeas <- ggplot(data=subset(tot.merge, bat_species=="Eidolon dupreanum")) +
  geom_rect(data=subset(cal.dat,bat_species=="Eidolon dupreanum"), aes(xmin=doy_start, xmax=doy_stop, ymin=22, ymax=38, fill=season), alpha=.3) +
  geom_violin(aes(x=month, y=GhV_ct, group=month),  #color="coral",
              draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA)+
  geom_jitter(aes(x=month, y=GhV_ct,color=roost_site, group = month_doy), width=.1, size=.1) +
  #geom_point(aes(x=doy, y=prevGhV, size=Ntot, shape=year), color="#00C094") + 
  #geom_point(aes(x=doy, y=prevAngV, size=Ntot, shape=year), color="#F8766D") + 
  scale_fill_manual(values=filz) +
  ggnewscale::new_scale_fill() + 
  #scale_x_continuous(breaks=c(31, 93, 155, 218, 280, 344 ), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) +
  #geom_linerange(aes(x=doy, ymin=prevGhV_lci, ymax=prevGhV_uci),color="#00C094", linetype=2, alpha=.3, show.legend = F) +
  #geom_linerange(aes(x=doy, ymin=prevAngV_lci, ymax=prevAngV_uci),color="#F8766D", linetype=2, alpha=.3, show.legend = F) +
  geom_ribbon(data=subset(pred.ct, antigen=="GhV" & bat_species=="Eidolon dupreanum"), aes(x=month, ymin= ct_lci, ymax=ct_uci, fill=antigen), alpha=.3, show.legend = F) +
  geom_line(data=subset(pred.ct, antigen=="GhV" & bat_species=="Eidolon dupreanum"), aes(x=month, y= ct, color=antigen)) +
  facet_grid(~bat_species) + theme_bw() +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),, 
        strip.background = element_rect(fill="white"), 
        strip.text.y = element_text(size=16), strip.text.x = element_text(size=16, face="italic"),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "lines"),
        plot.margin = unit(c(.5,.5,.5,1), "lines")) +
  #scale_shape_manual(values=yearz) +
  coord_cartesian(xlim = c(0, 14), expand = FALSE) + ylab("Ct value") + xlab("month of year") +
  scale_y_reverse() +
  scale_fill_manual(values=colorz) + scale_color_manual(values=colorz) + 
  guides(size=guide_legend(override.aes = list(color="black")), shape=guide_legend(override.aes = list(color="black")))
pCtSeas

