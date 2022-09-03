rm(list = ls())
require("ggplot2");require("tidyverse");require("dplyr");require("ggpubr");require("Hmisc");require("png");require("grid");require("ggrepel")
theme_set(
  theme_bw()
)

mydf_raw<-read_csv("C:/Users/rmoge/OneDrive - Indiana University/Lab.Notebook/20180515_lifespan.mutants/data/cases_lifespan.mutants.csv")
mydf<-as_tibble(mydf_raw)
mydf$Selection<-factor(mydf$Selection, levels=c("No Long","Long","Anc"))
mydf$Strain_factor<-factor(mydf$Strain, levels=c("2","4","5","9","11","16","19","0"))
mydf$Strain<-recode(mydf$Strain_factor, "2" = "rpl19aΔ", "4" = "pmr1Δ", "5" = "sch9Δ", "9"="ypt6Δ","11"="tor1Δ","16"="tif2Δ","19"="LYS2", "0"="WT")
mydf$Strain<-factor(mydf$Strain, levels=c("rpl19aΔ","pmr1Δ","sch9Δ","ypt6Δ","tor1Δ","tif2Δ","LYS2","WT"))
#mydf$History<-factor(mydf$History, levels=c("Home","Away"))
mycols<-c("black","black")

avgdf_raw<-read_csv("C:/Users/rmoge/OneDrive - Indiana University/Lab.Notebook/20180515_lifespan.mutants/data/cases_lifespan.mutants_avgs.csv")
avgdf<-as_tibble(avgdf_raw)
avgdf$Selection<-factor(avgdf$Selection, levels=c("No Long","Long","Anc"))
avgdf$Strain_factor<-factor(avgdf$Strain, levels=c("2","4","5","9","11","16","19","0"))
avgdf$Strain<-recode(avgdf$Strain_factor, "2" = "rpl19aΔ", "4" = "pmr1Δ", "5" = "sch9Δ", "9"="ypt6Δ","11"="tor1Δ","16"="tif2Δ","19"="LYS2", "0"="WT")
avgdf$Strain<-factor(avgdf$Strain, levels=c("rpl19aΔ","pmr1Δ","sch9Δ","ypt6Δ","tor1Δ","tif2Δ","LYS2","WT"))


ance0df<-read_csv("C:/Users/rmoge/OneDrive - Indiana University/Lab.Notebook/20180515_lifespan.mutants/data/CLS/cases_anc.ANOVA_e0.csv")
ance0df$strain_factor<-factor(ance0df$strain, levels=c("2","4","5","9","11","16","19","0"))
ance0df$strain<-recode(ance0df$strain_factor, "2" = "rpl19aΔ", "4" = "pmr1Δ", "5" = "sch9Δ", "9"="ypt6Δ","11"="tor1Δ","16"="tif2Δ","19"="LYS2", "0"="WT")
ance0df$strain<-factor(ance0df$strain, levels=c("rpl19aΔ","pmr1Δ","sch9Δ","ypt6Δ","tor1Δ","tif2Δ","LYS2","WT"))

ancrrdf<-read_csv("C:/Users/rmoge/OneDrive - Indiana University/Lab.Notebook/20180515_lifespan.mutants/data/reproduction.assays/Ancestor/cases_anc.ANOVA_rr.csv")
ancrrdf$Strain_factor<-factor(ancrrdf$Strain, levels=c("2","4","5","9","11","16","19","0"))
ancrrdf$Strain<-recode(ancrrdf$Strain_factor, "2" = "rpl19aΔ", "4" = "pmr1Δ", "5" = "sch9Δ", "9"="ypt6Δ","11"="tor1Δ","16"="tif2Δ","19"="LYS2", "0"="WT")
ancrrdf$Strain<-factor(ancrrdf$Strain, levels=c("rpl19aΔ","pmr1Δ","sch9Δ","ypt6Δ","tor1Δ","tif2Δ","LYS2","WT"))


myydf_full<-as_tibble(mydf)
myydf<-filter(myydf_full, Selection != "Anc")

anc.colors<-c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey")

#anc.e0<-c(3.966526375,3.966526375,3.966526375,3.966526375,3.966526375,2.653540408,2.653540408,2.653540408,2.653540408,2.653540408,3.646673218,3.646673218,3.646673218,3.646673218,3.646673218,4.147585939,4.147585939,4.147585939,4.147585939,4.147585939,3.1103613,3.1103613,3.1103613,3.1103613,3.1103613,4.241732325,4.241732325,4.241732325,4.241732325,4.241732325,5.754412658,5.754412658,5.754412658,5.754412658,5.754412658,1.987391586,1.987391586,1.987391586,1.987391586,1.987391586,3.966526375,3.966526375,3.966526375,3.966526375,3.966526375,2.653540408,2.653540408,2.653540408,2.653540408,2.653540408,3.646673218,3.646673218,3.646673218,3.646673218,3.646673218,4.147585939,4.147585939,4.147585939,4.147585939,4.147585939,3.1103613,3.1103613,3.1103613,3.1103613,3.1103613,4.241732325,4.241732325,4.241732325,4.241732325,4.241732325,5.754412658,5.754412658,5.754412658,5.754412658,5.754412658,1.987391586,1.987391586,1.987391586,1.987391586,1.987391586)

xlleg<-readPNG("~/GitHub/XevoLong/figures/20211118_144448_legend_thicker.png")
xllegend<-rasterGrob(xlleg, interpolate = TRUE)

enaught <- ggplot(myydf, aes(x=Strain, y=e0))
enaught +
  #geom_errorbar(aes(ymax=anc.e0,ymin=anc.e0),color=anc.colors,lwd = 2.2,linetype=117) +
  geom_errorbar(data=myydf, mapping = aes(x=Strain,ymin=anc.e0,ymax=anc.e0),color=anc.colors,lwd=4.2,linetype=117) +
  geom_jitter(aes(shape = Selection, color = Selection),
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.55),
              size = 6.5, stroke = 2.25,
              show.legend = FALSE
  ) +
  stat_summary(
    aes(color = Selection),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 2.25, stroke=2, shape=95,
    position = position_dodge(0.35),
    show.legend=FALSE
  ) +
  scale_color_manual(values = c("black","black"), labels = c("No Long","Long")) +
  scale_shape_manual(values = c(15,0), labels = c("No Long", "Long")) +#12 is a square with a vertical cross inside it
  labs(x="\nStrain",y="Life expectancy (d)") +
  scale_y_continuous(limits = c(-0.04,16.8),breaks=c(0,3,6,9,12,15), expand = c(0,0)) +
  annotation_custom(grob = xllegend, xmin=7.4, xmax=8.4, ymin=15.06, ymax=16.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))



#anc.reprod<-c(0.293052583,0.293052583,0.293052583,0.293052583,0.293052583,0.26764054,0.26764054,0.26764054,0.26764054,0.26764054,0.21795,0.21795,0.21795,0.21795,0.21795,0.246061944,0.246061944,0.246061944,0.246061944,0.246061944,0.352697718,0.352697718,0.352697718,0.352697718,0.352697718,0.34632925,0.34632925,0.34632925,0.34632925,0.34632925,0.359498994,0.359498994,0.359498994,0.359498994,0.359498994,0.354181633,0.354181633,0.354181633,0.354181633,0.354181633,0.293052583,0.293052583,0.293052583,0.293052583,0.293052583,0.26764054,0.26764054,0.26764054,0.26764054,0.26764054,0.21795,0.21795,0.21795,0.21795,0.21795,0.246061944,0.246061944,0.246061944,0.246061944,0.246061944,0.352697718,0.352697718,0.352697718,0.352697718,0.352697718,0.34632925,0.34632925,0.34632925,0.34632925,0.34632925,0.359498994,0.359498994,0.359498994,0.359498994,0.359498994,0.354181633,0.354181633,0.354181633,0.354181633,0.354181633)


reprod_notmatched <- ggplot(myydf, aes(x=Strain, y=off.per.hr_50.50))
reprod_notmatched +
  #geom_errorbar(aes(ymin=anc.reprod,ymax=anc.reprod),color=anc.colors,lwd = 2.25,linetype=117) +
  geom_errorbar(data=myydf, mapping = aes(x=Strain,ymin=anc.off.per.hr,ymax=anc.off.per.hr),color=anc.colors,lwd=4.2,linetype=117) +
  geom_jitter(aes(shape = Selection, color = Selection),
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.65),
              size = 6.5, stroke = 2.25
  ) +
  stat_summary(
    aes(color = Selection),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 2.25, stroke=2, shape=95,
    position = position_dodge(0.35),
    show.legend=FALSE
  ) +
  scale_color_manual(values = c("black","black"), labels = c("No Long","Long")) +
  scale_shape_manual(values = c(15,0), labels = c("No Long", "Long")) +#12 is a square with a vertical cross inside it
  labs(x="\nStrain",y="No. offspring per hr\n") +
  scale_y_continuous(limits = c(.39,.87),breaks=c(.4,.6,.8), expand = c(0,0)) +
  annotation_custom(grob = xllegend, xmin=7.4, xmax=8.4, ymin=.82, ymax=.865)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

reprod <- ggplot(myydf, aes(x=Strain, y=off.per.hr.matched))
reprod +
  #geom_errorbar(aes(ymin=anc.reprod,ymax=anc.reprod),color=anc.colors,lwd = 2.25,linetype=117) +
  geom_errorbar(data=myydf, mapping = aes(x=Strain,ymin=anc.off.per.hr,ymax=anc.off.per.hr),color=anc.colors,lwd=4.2,linetype=117) +
  geom_jitter(aes(shape = Selection, color = Selection),
              position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.55),
              size = 6.5, stroke = 2.25,
              show.legend=FALSE
  ) +
  stat_summary(
    aes(color = Selection),
    fun.data = "mean_se", fun.args = list(mult = (1)), #mean_sdl add +/- standard deviation; mult=1 means that it is SD*1 that is drawn.  Mean_se draws the standard error of the mean
    geom = "pointrange", size = 2.25, stroke=2, shape=95,
    position = position_dodge(0.35),
    show.legend=FALSE
  ) +
  scale_color_manual(values = c("black","black"), labels = c("No Long","Long")) +
  scale_shape_manual(values = c(15,0), labels = c("No Long", "Long")) +#12 is a square with a vertical cross inside it
  labs(x="\nStrain",y="No. offspring per hr\n") +
  scale_y_continuous(limits = c(.39,.87),breaks=c(.4,.6,.8), expand = c(0,0)) +
  annotation_custom(grob = xllegend, xmin=7.4, xmax=8.4, ymin=.82, ymax=.865)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

#Now t-tests for T50 No Long vs. T50 Long life expectancy
var.test(c(filter(myydf, Selection == "Long",Genotype=="rpl19a")$e0),c(filter(myydf, Selection == "No Long",Genotype=="rpl19a")$e0))
#P = .2344. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="rpl19a")$e0),c(filter(myydf, Selection == "No Long",Genotype=="rpl19a")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)
#p = 0.01503 < 0.05

var.test(c(filter(myydf, Selection == "Long",Genotype=="pmr1")$e0),c(filter(myydf, Selection == "No Long",Genotype=="pmr1")$e0))
#P = .7278. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="pmr1")$e0),c(filter(myydf, Selection == "No Long",Genotype=="pmr1")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)
#p = 0.0002782 < 0.001

var.test(c(filter(myydf, Selection == "Long",Genotype=="sch9")$e0),c(filter(myydf, Selection == "No Long",Genotype=="sch9")$e0))
#P = .3613. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="sch9")$e0),c(filter(myydf, Selection == "No Long",Genotype=="sch9")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)
#p = 0.01071 < 0.05

var.test(c(filter(myydf, Selection == "Long",Genotype=="ypt6")$e0),c(filter(myydf, Selection == "No Long",Genotype=="ypt6")$e0))
#P = .2920. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="ypt6")$e0),c(filter(myydf, Selection == "No Long",Genotype=="ypt6")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)
#p = 0.01205 < 0.05

var.test(c(filter(myydf, Selection == "Long",Genotype=="tor1")$e0),c(filter(myydf, Selection == "No Long",Genotype=="tor1")$e0))
#P = .1120. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="tor1")$e0),c(filter(myydf, Selection == "No Long",Genotype=="tor1")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)
#p = 0.001116 < 0.01

var.test(c(filter(myydf, Selection == "Long",Genotype=="tif2")$e0),c(filter(myydf, Selection == "No Long",Genotype=="tif2")$e0))
#P = .2484. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="tif2")$e0),c(filter(myydf, Selection == "No Long",Genotype=="tif2")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)
#p = 0.006902 < 0.01

var.test(c(filter(myydf, Selection == "Long",Genotype=="JD174")$e0),c(filter(myydf, Selection == "No Long",Genotype=="JD174")$e0))
#P = .4977. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="JD174")$e0),c(filter(myydf, Selection == "No Long",Genotype=="JD174")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)
#p = 0.01094 < 0.05

var.test(c(filter(myydf, Selection == "Long",Genotype=="BY4742")$e0),c(filter(myydf, Selection == "No Long",Genotype=="BY4742")$e0))
#P = .2512. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="BY4742")$e0),c(filter(myydf, Selection == "No Long",Genotype=="BY4742")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)
#p = 0.000668 < 0.001

my.p.vec.e0<-c(t.test(c(filter(myydf, Selection == "Long",Genotype=="rpl19a")$e0),c(filter(myydf, Selection == "No Long",Genotype=="rpl19a")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="pmr1")$e0),c(filter(myydf, Selection == "No Long",Genotype=="pmr1")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="sch9")$e0),c(filter(myydf, Selection == "No Long",Genotype=="sch9")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="ypt6")$e0),c(filter(myydf, Selection == "No Long",Genotype=="ypt6")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="tor1")$e0),c(filter(myydf, Selection == "No Long",Genotype=="tor1")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="tif2")$e0),c(filter(myydf, Selection == "No Long",Genotype=="tif2")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="JD174")$e0),c(filter(myydf, Selection == "No Long",Genotype=="JD174")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="BY4742")$e0),c(filter(myydf, Selection == "No Long",Genotype=="BY4742")$e0),paired=F,alternative = "greater",mu=0,var.equal=T)$p.value)
print(my.p.vec.e0)

p.adjust(p=my.p.vec.e0, method = "BH")#cite BH&Y 2009
#We now have the following adjusted p-values
#rpl19a:  P = 0.0150 < 0.05
#pmr1:    P = 0.00223 < 0.01
#sch9:    P = 0.0138 < 0.05
#ypt6:    P = 0.0138 < 0.05
#tor1:    P = 0.00298 < 0.01
#tif2:    P = 0.0138 < 0.05
#JD174:   P = 0.0138 < 0.05
#BY4742:  P = 0.00267 < 0.01


#Now t-tests for T46 or 36 No Long vs. T50 Long reproduction rate
var.test(c(filter(myydf, Selection == "Long",Genotype=="rpl19a")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="rpl19a")$off.per.hr.matched))
#P = .2925. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="rpl19a")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="rpl19a")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=T)
#p = 3.726EE-6 < 0.05

var.test(c(filter(myydf, Selection == "Long",Genotype=="pmr1")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="pmr1")$off.per.hr.matched))
#P = .1377. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="pmr1")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="pmr1")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=T)
#p = 0.03467 < 0.05

var.test(c(filter(myydf, Selection == "Long",Genotype=="sch9")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="sch9")$off.per.hr.matched))
#P = .4903. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="sch9")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="sch9")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=T)
#p = 0.2101 --> n.s.

var.test(c(filter(myydf, Selection == "Long",Genotype=="ypt6")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="ypt6")$off.per.hr.matched))
#P = .0093. Heteroscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="ypt6")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="ypt6")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=F)
#p = 0.04991 < 0.05

var.test(c(filter(myydf, Selection == "Long",Genotype=="tor1")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="tor1")$off.per.hr.matched))
#P = .9625. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="tor1")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="tor1")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=T)
#p = 0.0006366 < 0.001

var.test(c(filter(myydf, Selection == "Long",Genotype=="tif2")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="tif2")$off.per.hr.matched))
#P = .05892. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="tif2")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="tif2")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=T)
#p = 0.007245 < 0.01

var.test(c(filter(myydf, Selection == "Long",Genotype=="JD174")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="JD174")$off.per.hr.matched))
#P = .02214. Heteroscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="JD174")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="JD174")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=F)
#p = 0.001932 < 0.01

var.test(c(filter(myydf, Selection == "Long",Genotype=="BY4742")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="BY4742")$off.per.hr.matched))
#P = .8486. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="BY4742")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="BY4742")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=T)
#p = 0.001091 < 0.01

my.p.vec.reprod<-c(t.test(c(filter(myydf, Selection == "Long",Genotype=="rpl19a")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="rpl19a")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="pmr1")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="pmr1")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="sch9")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="sch9")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="ypt6")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="ypt6")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=F)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="tor1")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="tor1")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="tif2")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="tif2")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="JD174")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="JD174")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=F)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="BY4742")$off.per.hr.matched),c(filter(myydf, Selection == "No Long",Genotype=="BY4742")$off.per.hr.matched),paired=F,alternative = "less",mu=0,var.equal=T)$p.value)

print(my.p.vec.reprod)
p.adjust(p=my.p.vec.reprod, method = "BH")#cite BH&Y 2009
#We now have the following adjusted p-values

#rpl19a:  P = 0.0000298 < 0.001
#pmr1:    P = 0.0462 < 0.05
#sch9:    P = 0.2101 --> n.s.
#ypt6:    P = 0.0570 < 0.1
#tor1:    P = 0.00255 < 0.01
#tif2:    P = 0.0116 < 0.05
#JD174:   P = 0.00386 < 0.01
#BY4742:  P = 0.00291 < 0.01

#Now one-sample t-tests to evaluate the fulfillment of criterion #1.
t.test(x=c(filter(myydf, Selection == "Long",Genotype=="rpl19a")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)
crit1p1<-t.test(x=c(filter(myydf, Selection == "Long",Genotype=="rpl19a")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)$p.value
t.test(x=c(filter(myydf, Selection == "Long",Genotype=="pmr1")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)
crit1p2<-t.test(x=c(filter(myydf, Selection == "Long",Genotype=="pmr1")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)$p.value
t.test(x=c(filter(myydf, Selection == "Long",Genotype=="sch9")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)
crit1p3<-t.test(x=c(filter(myydf, Selection == "Long",Genotype=="sch9")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)$p.value
t.test(x=c(filter(myydf, Selection == "Long",Genotype=="ypt6")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)
crit1p4<-t.test(x=c(filter(myydf, Selection == "Long",Genotype=="ypt6")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)$p.value
t.test(x=c(filter(myydf, Selection == "Long",Genotype=="tor1")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)
crit1p5<-t.test(x=c(filter(myydf, Selection == "Long",Genotype=="tor1")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)$p.value
t.test(x=c(filter(myydf, Selection == "Long",Genotype=="tif2")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)
crit1p6<-t.test(x=c(filter(myydf, Selection == "Long",Genotype=="tif2")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)$p.value
t.test(x=c(filter(myydf, Selection == "Long",Genotype=="JD174")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)
crit1p7<-t.test(x=c(filter(myydf, Selection == "Long",Genotype=="JD174")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)$p.value
t.test(x=c(filter(myydf, Selection == "Long",Genotype=="BY4742")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)
crit1p8<-t.test(x=c(filter(myydf, Selection == "Long",Genotype=="BY4742")$deltal),y=NULL,alternative = "g",mu=0,paired=F,var.equal=F)$p.value

crit1pvec<-c(crit1p1,crit1p2,crit1p3,crit1p4,crit1p5,crit1p6,crit1p7,crit1p8)
print(crit1pvec)
crit1pvec.adjust<-p.adjust(p=crit1pvec, method = "BH")#cite BH&Y 2009
crit1pvec.adjust
#[1] 0.159380982 0.011749436 0.387123569 0.075379888 0.009069475 0.011749436 0.233487463 0.009069475

#######################################################
#Now t-tests for T50 No Long vs. T50 Long reproduction rate
var.test(c(filter(myydf, Selection == "Long",Genotype=="rpl19a")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="rpl19a")$off.per.hr_50.50))
#P = .436. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="rpl19a")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="rpl19a")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=T)
#p = 8.03EE-06 < 0.0001

var.test(c(filter(myydf, Selection == "Long",Genotype=="pmr1")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="pmr1")$off.per.hr_50.50))
#P = .4583. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="pmr1")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="pmr1")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=T)
#p = 0.07809 < 0.1

var.test(c(filter(myydf, Selection == "Long",Genotype=="sch9")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="sch9")$off.per.hr_50.50))
#P = .1775. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="sch9")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="sch9")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=T)
#p = 0.3296 --> n.s.

var.test(c(filter(myydf, Selection == "Long",Genotype=="ypt6")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="ypt6")$off.per.hr_50.50))
#P = .003668. Heteroscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="ypt6")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="ypt6")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=F)
#p = 0.02674 < 0.05

var.test(c(filter(myydf, Selection == "Long",Genotype=="tor1")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="tor1")$off.per.hr_50.50))
#P = .1084. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="tor1")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="tor1")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=T)
#p = 0.008224 < 0.01

var.test(c(filter(myydf, Selection == "Long",Genotype=="tif2")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="tif2")$off.per.hr_50.50))
#P = .0201. Heteroscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="tif2")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="tif2")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=F)
#p = 0.005643 < 0.01

var.test(c(filter(myydf, Selection == "Long",Genotype=="JD174")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="JD174")$off.per.hr_50.50))
#P = .02185. Heteroscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="JD174")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="JD174")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=F)
#p = 0.005708 < 0.01

var.test(c(filter(myydf, Selection == "Long",Genotype=="BY4742")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="BY4742")$off.per.hr_50.50))
#P = .9. Homoscedastic.
t.test(c(filter(myydf, Selection == "Long",Genotype=="BY4742")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="BY4742")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=T)
#p = 0.002452 < 0.01

my.p.vec.reprod50.50<-c(t.test(c(filter(myydf, Selection == "Long",Genotype=="rpl19a")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="rpl19a")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="pmr1")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="pmr1")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="sch9")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="sch9")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="ypt6")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="ypt6")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=F)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="tor1")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="tor1")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="tif2")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="tif2")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=T)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="JD174")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="JD174")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=F)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="BY4742")$off.per.hr_50.50),c(filter(myydf, Selection == "No Long",Genotype=="BY4742")$off.per.hr_50.50),paired=F,alternative = "less",mu=0,var.equal=T)$p.value)

print(my.p.vec.reprod50.50)
p.adjust(p=my.p.vec.reprod50.50, method = "BH")#cite BH&Y 2009
#We now have the following adjusted p-values
#rpl19a:  P = 0.0000642 < 0.0001
#pmr1:    P = 0.0893 < 0.1
#sch9:    P = 0.330 --> n.s.
#ypt6:    P = 0.0356 < 0.05
#tor1:    P = 0.0132 < 0.05
#tif2:    P = 0.00654 < 0.01
#JD174:   P = 0.0114 < 0.05
#BY4742:  P = 0.00654 < 0.01



###############################################################################
#Now 2-way ANOVAs for all of the measured traits
e02way<-aov(myydf$e0 ~ myydf$Selection + myydf$Genotype + myydf$Selection*myydf$Genotype)
summary(e02way)
#TukeyHSD(e02way)
#Significant effect of selection regime (P = 4.42EE-14)
#Significant effect of strain           (P = 4.01EE-09)
#Significant interaction effect         (P = .00595)

reprod2way<-aov(myydf$off.per.hr.matched ~ myydf$Selection + myydf$Genotype + myydf$Selection*myydf$Genotype)
summary(reprod2way)
#TukeyHSD(reprod2way)

#Significant effect of selection regime (P = 7.48EE-13)
#Significant effect of strain           (P = 1.72EE-12)
#Significant interaction effect      (P = .00399)

reprod2way50.50<-aov(myydf$off.per.hr_50.50 ~ myydf$Selection + myydf$Genotype + myydf$Selection*myydf$Genotype)
summary(reprod2way50.50)
#TukeyHSD(reprod2way50.50)
#Significant effect of selection regime (P = 7.15EE-13)
#Significant effect of strain           (P = 4.95EE-12)
#Significant interaction effect         (P = .00179)


###############################################################################
#Now 2-way ANOVAs for all of the measured traits' DELTA VALUES. Changes vs. the ancor foeach strain
De02way<-aov(myydf$deltal ~ myydf$Selection + myydf$Genotype + myydf$Selection*myydf$Genotype)
summary(De02way)
#TukeyHSD(De02way)
#Significant effect of selection regime (P = 4.42EE-14)
#Significant effect of strain           (P = 5.99EE-12)
#Significant interaction effect         (P = .00595)

Dreprod2way<-aov(myydf$deltar ~ myydf$Selection + myydf$Genotype + myydf$Selection*myydf$Genotype)
summary(Dreprod2way)
#TukeyHSD(Dreprod2way)

#Significant effect of selection regime (P = 7.48EE-13)
#Significant effect of strain           (P = 2.07EE-14)
#No significant interaction effect      (P = .00399)

Dreprod2way50.50<-aov(myydf$deltar.50.50 ~ myydf$Selection + myydf$Genotype + myydf$Selection*myydf$Genotype)
summary(Dreprod2way50.50)
#TukeyHSD(Dreprod2way50.50)
#Significant effect of selection regime (P = 7.15EE-13)
#Significant effect of strain           (P = 4.18EE-16)
#Significant interaction effect         (P = .00179)


#########################################################################################################
#Statistics on strain ancestors

#First e0. ANOVA among all the strain ancestors
ance0aov<-aov(ance0df$e0 ~ ance0df$strain)
summary(ance0aov)
TukeyHSD(ance0aov)
#No significant differences, P = 0.127 (my favorite number is 127!)

#Second, reproduction rate. ANOVA among all the strain ancestors
ancrraov<-aov(ancrrdf$off.per.hr ~ ancrrdf$Strain)
summary(ancrraov)
TukeyHSD(ancrraov)
#Ssignificant differences, P < 0.05


#########################################################################################################
#Use one-sample t-tests to compare reproduction of 16, 19, and 0 L popns to their ancestors
t.test(c(filter(myydf, Selection == "Long",Genotype=="tif2")$off.per.hr.matched),alternative = "two.sided",mu=0.681405478)
#p = 0.09747 < 0.05
t.test(c(filter(myydf, Selection == "Long",Genotype=="JD174")$off.per.hr.matched),alternative = "two.sided",mu=0.706158003)
#p = 0.01377 < 0.05
t.test(c(filter(myydf, Selection == "Long",Genotype=="BY4742")$off.per.hr.matched),alternative = "two.sided",mu=0.694369123)
#p = 0.01829 < 0.05

my.p.vec.reprod.vs.anc.test<-c(t.test(c(filter(myydf, Selection == "Long",Genotype=="tif2")$off.per.hr.matched),alternative = "two.sided",mu=0.681405478)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="JD174")$off.per.hr.matched),alternative = "two.sided",mu=0.706158003)$p.value,t.test(c(filter(myydf, Selection == "Long",Genotype=="BY4742")$off.per.hr.matched),alternative = "two.sided",mu=0.694369123)$p.value)
print(my.p.vec.reprod.vs.anc.test)
p.adjust(p=my.p.vec.reprod.vs.anc.test, method = "BH")#cite BH&Y 2009



################################################################################
#####Bivariate plots: I need to make the bivariate plots that Farrah and Matt were interested in.
rpl19adf<-filter(myydf, Genotype=="rpl19a")
pmr1df<-filter(myydf, Genotype=="pmr1")
sch9df<-filter(myydf, Genotype=="sch9")
ypt6df<-filter(myydf, Genotype=="ypt6")
tor1df<-filter(myydf, Genotype=="tor1")
tif2df<-filter(myydf, Genotype=="tif2")
LYS2df<-filter(myydf, Genotype=="JD174")
WTdf<-filter(myydf, Genotype=="BY4742")

label2<- textGrob(label = "rpl19aΔ", x = 0.05, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))
label4<- textGrob(label = "pmr1Δ", x = 0.1, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))
label5<- textGrob(label = "sch9Δ", x = 0.05, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))
label9<- textGrob(label = "ypt6Δ", x = 0.05, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))
label11<- textGrob(label = "tor1Δ", x = 0.05, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))
label16<- textGrob(label = "tif2Δ", x = 0.05, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))
label19<- textGrob(label = "LYS2", x = 0.1, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))
label0<- textGrob(label = "WT", x = 0.05, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))



rpl19aDbiv <- ggplot(rpl19adf, aes(x=deltar, y=deltal))
rpl19aDbiv +
  geom_point(aes(shape = Selection, color = Selection),size = 6.5, stroke = 2.25) + # Show dots
  #geom_label_repel(label=rpl19a$Selection, nudge_x = 0.02, nudge_y = 1)+
  #geom_errorbar(data=myydf, mapping = aes(x=Strain,ymin=anc.off.per.hr,ymax=anc.off.per.hr),color=anc.colors,lwd=4.2,linetype=117) +
  #geom_jitter(aes(shape = Selection, color = Selection),position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.55),size = 6.5, stroke = 2.25,show.legend=TRUE) +
  #stat_summary(aes(color = Selection),fun.data = "mean_se", fun.args = list(mult = (1)),geom = "pointrange", size = 2.25, stroke=2, shape=95, position = position_dodge(0.35),show.legend=FALSE ) +
  scale_color_manual(values = c("black","black"), labels = c("No Long","Long")) +
  scale_shape_manual(values = c(15,0), labels = c("No Long", "Long")) +#12 is a square with a vertical cross inside it
  labs(x="\nΔR",y="ΔL\n") +
  #scale_y_continuous(limits = c(-3.75,3.75),breaks=c(-3,-2,-1,0,1,2,3), expand = c(0,0)) +
  #scale_x_continuous(limits = c(-.25,.25),breaks=c(-.2,-.1,0,.1,.2), expand = c(0,0)) +
  #annotation_custom(grob = xllegend, xmin=7.4, xmax=8.4, ymin=.82, ymax=.865)+
  annotation_custom(label2, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  #annotate("text",x=0.68,y=1.5,label="rpl19aΔ",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))


pmr1Dbiv <- ggplot(pmr1df, aes(x=deltar, y=deltal))
pmr1Dbiv +
  geom_point(aes(shape = Selection, color = Selection),size = 6.5, stroke = 2.25) + # Show dots
  scale_color_manual(values = c("black","black"), labels = c("No Long","Long")) +
  scale_shape_manual(values = c(15,0), labels = c("No Long", "Long")) +#12 is a square with a vertical cross inside it
  labs(x="\nΔR",y="ΔL\n") +
  scale_y_continuous(limits = c(-2.35,2.55),breaks=c(-2,-1,0,1,2), expand = c(0,0)) +
  scale_x_continuous(limits = c(-.08,.16),breaks=c(-.05,0,.05,.1,.15), expand = c(0,0)) +
  #annotation_custom(grob = xllegend, xmin=7.4, xmax=8.4, ymin=.82, ymax=.865)+
  annotation_custom(label4, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

sch9Dbiv <- ggplot(sch9df, aes(x=deltar, y=deltal))
sch9Dbiv +
  geom_point(aes(shape = Selection, color = Selection),size = 6.5, stroke = 2.25) + # Show dots
  scale_color_manual(values = c("black","black"), labels = c("No Long","Long")) +
  scale_shape_manual(values = c(15,0), labels = c("No Long", "Long")) +#12 is a square with a vertical cross inside it
  labs(x="\nΔR",y="ΔL\n") +
  #scale_y_continuous(limits = c(-3.75,3.75),breaks=c(-3,-2,-1,0,1,2,3), expand = c(0,0)) +
  scale_x_continuous(limits = c(-.01,.17),breaks=c(0,.05,.1,.15), expand = c(0,0)) +
  #annotation_custom(grob = xllegend, xmin=7.4, xmax=8.4, ymin=.82, ymax=.865)+
  annotation_custom(label5, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

ypt6Dbiv <- ggplot(ypt6df, aes(x=deltar, y=deltal))
ypt6Dbiv +
  geom_point(aes(shape = Selection, color = Selection),size = 6.5, stroke = 2.25) + # Show dots
  scale_color_manual(values = c("black","black"), labels = c("No Long","Long")) +
  scale_shape_manual(values = c(15,0), labels = c("No Long", "Long")) +#12 is a square with a vertical cross inside it
  labs(x="\nΔR",y="ΔL\n") +
  #scale_y_continuous(limits = c(-3.75,3.75),breaks=c(-3,-2,-1,0,1,2,3), expand = c(0,0)) +
  #scale_x_continuous(limits = c(-.01,.17),breaks=c(0,.05,.1,.15), expand = c(0,0)) +
  #annotation_custom(grob = xllegend, xmin=7.4, xmax=8.4, ymin=.82, ymax=.865)+
  annotation_custom(label9, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

tor1Dbiv <- ggplot(tor1df, aes(x=deltar, y=deltal))
tor1Dbiv +
  geom_point(aes(shape = Selection, color = Selection),size = 6.5, stroke = 2.25) + # Show dots
  scale_color_manual(values = c("black","black"), labels = c("No Long","Long")) +
  scale_shape_manual(values = c(15,0), labels = c("No Long", "Long")) +#12 is a square with a vertical cross inside it
  labs(x="\nΔR",y="ΔL\n") +
  #scale_y_continuous(limits = c(-3.75,3.75),breaks=c(-3,-2,-1,0,1,2,3), expand = c(0,0)) +
  #scale_x_continuous(limits = c(-.01,.17),breaks=c(0,.05,.1,.15), expand = c(0,0)) +
  #annotation_custom(grob = xllegend, xmin=7.4, xmax=8.4, ymin=.82, ymax=.865)+
  annotation_custom(label11, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

tif2Dbiv <- ggplot(tif2df, aes(x=deltar, y=deltal))
tif2Dbiv +
  geom_point(aes(shape = Selection, color = Selection),size = 6.5, stroke = 2.25) + # Show dots
  scale_color_manual(values = c("black","black"), labels = c("No Long","Long")) +
  scale_shape_manual(values = c(15,0), labels = c("No Long", "Long")) +#12 is a square with a vertical cross inside it
  labs(x="\nΔR",y="ΔL\n") +
  #scale_y_continuous(limits = c(-3.75,3.75),breaks=c(-3,-2,-1,0,1,2,3), expand = c(0,0)) +
  #scale_x_continuous(limits = c(-.01,.17),breaks=c(0,.05,.1,.15), expand = c(0,0)) +
  #annotation_custom(grob = xllegend, xmin=7.4, xmax=8.4, ymin=.82, ymax=.865)+
  annotation_custom(label16, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

LYS2Dbiv <- ggplot(LYS2df, aes(x=deltar, y=deltal))
LYS2Dbiv +
  geom_point(aes(shape = Selection, color = Selection),size = 6.5, stroke = 2.25) + # Show dots
  scale_color_manual(values = c("black","black"), labels = c("No Long","Long")) +
  scale_shape_manual(values = c(15,0), labels = c("No Long", "Long")) +#12 is a square with a vertical cross inside it
  labs(x="\nΔR",y="ΔL\n") +
  #scale_y_continuous(limits = c(-3.75,3.75),breaks=c(-3,-2,-1,0,1,2,3), expand = c(0,0)) +
  #scale_x_continuous(limits = c(-.01,.17),breaks=c(0,.05,.1,.15), expand = c(0,0)) +
  #annotation_custom(grob = xllegend, xmin=7.4, xmax=8.4, ymin=.82, ymax=.865)+
  annotation_custom(label19, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))


WTDbiv <- ggplot(WTdf, aes(x=deltar, y=deltal))
WTDbiv +
  geom_point(aes(shape = Selection, color = Selection),size = 6.5, stroke = 2.25) + # Show dots
  scale_color_manual(values = c("black","black"), labels = c("No Long","Long")) +
  scale_shape_manual(values = c(15,0), labels = c("No Long", "Long")) +#12 is a square with a vertical cross inside it
  labs(x="\nΔR",y="ΔL\n") +
  #scale_y_continuous(limits = c(-3.75,3.75),breaks=c(-3,-2,-1,0,1,2,3), expand = c(0,0)) +
  #scale_x_continuous(limits = c(-.01,.17),breaks=c(0,.05,.1,.15), expand = c(0,0)) +
  #annotation_custom(grob = xllegend, xmin=7.4, xmax=8.4, ymin=.82, ymax=.865)+
  annotation_custom(label0, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))


################################################################################
#####Now plot all 8 strain ancestors together.
ancdf<-filter(myydf_full, Selection=="Anc")
labelanc<- textGrob(label = "Strain ancestors", x = 0.5, y = 0.97,just = c("center", "top"),gp=gpar(col = "#800020",fontsize = 32))

ancbiv <- ggplot(ancdf, aes(x=anc.off.per.hr, y=anc.e0))
ancbiv +
  #geom_point(aes(shape = Selection, color = Selection),size = 6.5, stroke = 2.25) + # Show dots
  geom_label(label=ancdf$Strain, nudge_x = 0.0, nudge_y = 0.00,fontface="bold",size=8.5,fill="white",color="#800020")+
  #geom_errorbar(data=myydf, mapping = aes(x=Strain,ymin=anc.off.per.hr,ymax=anc.off.per.hr),color=anc.colors,lwd=4.2,linetype=117) +
  #geom_jitter(aes(shape = Selection, color = Selection),position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.55),size = 6.5, stroke = 2.25,show.legend=TRUE) +
  #stat_summary(aes(color = Selection),fun.data = "mean_se", fun.args = list(mult = (1)),geom = "pointrange", size = 2.25, stroke=2, shape=95, position = position_dodge(0.35),show.legend=FALSE ) +
  scale_color_manual(values = c("black","black"), labels = c("No Long","Long")) +
  scale_shape_manual(values = c(15,0), labels = c("No Long", "Long")) +#12 is a square with a vertical cross inside it
  labs(x="\nNo. offspring per hr",y="Life expectancy (d)\n") +
  #scale_y_continuous(limits = c(-3.75,3.75),breaks=c(-3,-2,-1,0,1,2,3), expand = c(0,0)) +
  #scale_x_continuous(limits = c(1.75,6.1),breaks=c(2,3,4,5,6), expand = c(0,0)) +
  #annotation_custom(grob = xllegend, xmin=7.4, xmax=8.4, ymin=.82, ymax=.865)+
  annotation_custom(labelanc, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))




################################################################################
#####Now make plots of the avg values of the evolved populations.
#####Bivariate plots: I need to make the bivariate plots that Farrah and Matt were interested in.
rpl19adf<-filter(myydf, Genotype=="rpl19a")
pmr1df<-filter(myydf, Genotype=="pmr1")
sch9df<-filter(myydf, Genotype=="sch9")
ypt6df<-filter(myydf, Genotype=="ypt6")
tor1df<-filter(myydf, Genotype=="tor1")
tif2df<-filter(myydf, Genotype=="tif2")
LYS2df<-filter(myydf, Genotype=="JD174")
WTdf<-filter(myydf, Genotype=="BY4742")

label2<- textGrob(label = "rpl19aΔ", x = 0.05, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))
label4<- textGrob(label = "pmr1Δ", x = 0.1, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))
label5<- textGrob(label = "sch9Δ", x = 0.05, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))
label9<- textGrob(label = "ypt6Δ", x = 0.05, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))
label11<- textGrob(label = "tor1Δ", x = 0.05, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))
label16<- textGrob(label = "tif2Δ", x = 0.05, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))
label19<- textGrob(label = "LYS2", x = 0.1, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))
label0<- textGrob(label = "WT", x = 0.05, y = 0.95,just = c("left", "top"),gp=gpar(col = "black",fontsize = 32))


evoavgdf<-filter(avgdf,Selection!="Anc")
evoavgdf$Selection<-factor(evoavgdf$Selection, levels=c("No Long","Long"))
NLavgdf<-filter(avgdf,Selection=="No Long")
Lavgdf<-filter(avgdf,Selection=="Long")



evoavgcol1<-c("black","black","black","black","black","black","black","black","white","white","white","white","white","white","white","white")
evoavgcol2<-c("white","white","white","white","white","white","white","white","black","black","black","black","black","black","black","black")
evoavgcol3<-c("black","black","black","black","black","black","black","black")
evoavgcol4<-c("white","white","white","white","white","white","white","white")


evoavgDbiv <- ggplot(evoavgdf, aes(x=avg_deltar, y=avg_deltal))
evoavgDbiv +
  #geom_point(aes(shape = Selection, color = Selection),size = 6.5, stroke = 2.25) + # Show dots
  geom_label(label=evoavgdf$Strain, nudge_x = 0.00, nudge_y = 0.0,fontface="bold",size=8.5,color=evoavgcol1,fill=evoavgcol2,segment.size=0)+
  #geom_errorbar(data=myydf, mapping = aes(x=Strain,ymin=anc.off.per.hr,ymax=anc.off.per.hr),color=anc.colors,lwd=4.2,linetype=117) +
  #geom_jitter(aes(shape = Selection, color = Selection),position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.55),size = 6.5, stroke = 2.25,show.legend=TRUE) +
  #stat_summary(aes(color = Selection),fun.data = "mean_se", fun.args = list(mult = (1)),geom = "pointrange", size = 2.25, stroke=2, shape=95, position = position_dodge(0.35),show.legend=FALSE ) +
  scale_color_manual(values = c("black","black"), labels = c("No Long","Long")) +
  scale_shape_manual(values = c(15,0), labels = c("No Long", "Long")) +#12 is a square with a vertical cross inside it
  labs(x="\nΔR",y="ΔL\n") +
  #scale_y_continuous(limits = c(-3.75,3.75),breaks=c(-3,-2,-1,0,1,2,3), expand = c(0,0)) +
  #scale_x_continuous(limits = c(-.25,.25),breaks=c(-.2,-.1,0,.1,.2), expand = c(0,0)) +
  #annotation_custom(grob = xllegend, xmin=7.4, xmax=8.4, ymin=.82, ymax=.865)+
  #annotation_custom(label2, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  #annotate("text",x=0.68,y=1.5,label="rpl19aΔ",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

NLavgDbiv <- ggplot(NLavgdf, aes(x=avg_deltar, y=avg_deltal))
NLavgDbiv +
  #geom_point(aes(shape = Selection, color = Selection),size = 6.5, stroke = 2.25) + # Show dots
  geom_label(label=NLavgdf$Strain, nudge_x = 0.00, nudge_y = 0.0,fontface="bold",size=8.5,color=evoavgcol4,fill=evoavgcol3,segment.size=0)+
  #geom_errorbar(data=myydf, mapping = aes(x=Strain,ymin=anc.off.per.hr,ymax=anc.off.per.hr),color=anc.colors,lwd=4.2,linetype=117) +
  #geom_jitter(aes(shape = Selection, color = Selection),position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.55),size = 6.5, stroke = 2.25,show.legend=TRUE) +
  #stat_summary(aes(color = Selection),fun.data = "mean_se", fun.args = list(mult = (1)),geom = "pointrange", size = 2.25, stroke=2, shape=95, position = position_dodge(0.35),show.legend=FALSE ) +
  scale_color_manual(values = c("black","black"), labels = c("No Long","Long")) +
  scale_shape_manual(values = c(15,0), labels = c("No Long", "Long")) +#12 is a square with a vertical cross inside it
  labs(x="\nΔR",y="ΔL\n") +
  #scale_y_continuous(limits = c(-3.75,3.75),breaks=c(-3,-2,-1,0,1,2,3), expand = c(0,0)) +
  #scale_x_continuous(limits = c(-.25,.25),breaks=c(-.2,-.1,0,.1,.2), expand = c(0,0)) +
  #annotation_custom(grob = xllegend, xmin=7.4, xmax=8.4, ymin=.82, ymax=.865)+
  #annotation_custom(label2, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  #annotate("text",x=0.68,y=1.5,label="rpl19aΔ",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))

LavgDbiv <- ggplot(Lavgdf, aes(x=avg_deltar, y=avg_deltal))
LavgDbiv +
  #geom_point(aes(shape = Selection, color = Selection),size = 6.5, stroke = 2.25) + # Show dots
  geom_label(label=Lavgdf$Strain, nudge_x = 0.00, nudge_y = 0.0,fontface="bold",size=8.5,color=evoavgcol3,fill=evoavgcol4,segment.size=0)+
  #geom_errorbar(data=myydf, mapping = aes(x=Strain,ymin=anc.off.per.hr,ymax=anc.off.per.hr),color=anc.colors,lwd=4.2,linetype=117) +
  #geom_jitter(aes(shape = Selection, color = Selection),position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.55),size = 6.5, stroke = 2.25,show.legend=TRUE) +
  #stat_summary(aes(color = Selection),fun.data = "mean_se", fun.args = list(mult = (1)),geom = "pointrange", size = 2.25, stroke=2, shape=95, position = position_dodge(0.35),show.legend=FALSE ) +
  scale_color_manual(values = c("black","black"), labels = c("No Long","Long")) +
  scale_shape_manual(values = c(15,0), labels = c("No Long", "Long")) +#12 is a square with a vertical cross inside it
  labs(x="\nΔR",y="ΔL\n") +
  #scale_y_continuous(limits = c(-3.75,3.75),breaks=c(-3,-2,-1,0,1,2,3), expand = c(0,0)) +
  #scale_x_continuous(limits = c(-.25,.25),breaks=c(-.2,-.1,0,.1,.2), expand = c(0,0)) +
  #annotation_custom(grob = xllegend, xmin=7.4, xmax=8.4, ymin=.82, ymax=.865)+
  #annotation_custom(label2, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  #annotate("text",x=0.68,y=1.5,label="rpl19aΔ",size=12,color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black")) +
  theme(legend.key=element_blank(),axis.text.x=element_text(size=34,margin = margin(t = 20, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=34,margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title=element_text(size=36),legend.text=element_text(size=22),legend.title = element_text(size=34), axis.ticks.y = element_line(color = "black", size = 3), axis.ticks.x = element_line(color = "black", size = 3), axis.ticks.length.x=unit(0.7, "cm"),axis.ticks.length.y=unit(0.7, "cm"), legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 3.5), axis.line = element_line(colour = "black"))