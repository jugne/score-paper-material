######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library(coda)
library(cowplot)

# clear workspace
rm(list = ls())

# out <- "root_filter_max"
out <- "segment_root_filter_max"
# out_masc <- "segment_filter_max"
filter <- T
# filter<- F
# out <- "no_filter"

dir.create(out)

this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)

rates = data.frame(method = character(), run=integer(), rea_A=c(), rea_G=c(), mig_GtA=c(), mig_AtG=c(), Ne_G=c(), Ne_A=c())

clock_rates <- data.frame(run=integer(), method=character(), clock=c())
rea_rates <- data.frame(run=integer(), state=character(), rea=c())
a_mig_rates <-data.frame(method=character(), run=character(), mig=c())
a_ne_rates <-data.frame(method=character(), run=character(), ne=c())

g_mig_rates <-data.frame(method=character(), run=character(), mig=c())
g_ne_rates <-data.frame(method=character(), run=character(), ne=c())

for (run in seq(1,10)){
if (filter){
  logs_score <- list.files(path=paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/output", "/",out), pattern="h5n1_combined.rootFilter.log$", full.names = TRUE)
  logs_mascot <- list.files(path=paste0("../FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_",run,"/output", "/",out_masc), pattern="h5n1_combined.rootFilter.log$", full.names = TRUE)
} else{
  logs_score <- list.files(path=paste0("../FullDatesInference/current/h5n1_reject_30DaysDiff_08_16/run_",run,"/combined"), pattern="h5n1_combined.log$", full.names = TRUE)
  logs_mascot <- list.files(path=paste0("../FullDatesInference/current/h5n1_mascot_30DaysDiff_08_16/run_",run,"/combined"), pattern="h5n1_combined.log$", full.names = TRUE)
}
  
# read in the log files
t_s <- read.table(logs_score[[1]], header=TRUE, sep="\t")
t_m <- read.table(logs_mascot[[1]], header=TRUE, sep="\t")

# get score rates
# s.rea_G = data.frame(estimated=mean(t_s$reassortmentRate.Galliformes), 
#                        upper=quantile(t_s$reassortmentRate.Galliformes,0.975), lower=quantile(t_s$reassortmentRate.Galliformes,0.025) )
# s.rea_A = data.frame(estimated=mean(t_s$reassortmentRate.Anseriformes), 
#                        upper=quantile(t_s$reassortmentRate.Anseriformes,0.975), lower=quantile(t_s$reassortmentRate.Anseriformes,0.025) )
# 
# s.mig_GtA = data.frame(estimated=mean(t_s$b_migrationRate.Galliformes_to_Anseriformes), 
#                        upper=quantile(t_s$b_migrationRate.Galliformes_to_Anseriformes,0.975), lower=quantile(t_s$b_migrationRate.Galliformes_to_Anseriformes,0.025) )
# s.mig_AtG = data.frame(estimated=mean(t_s$b_migrationRate.Anseriformes_to_Galliformes), 
#                        upper=quantile(t_s$b_migrationRate.Anseriformes_to_Galliformes,0.975), lower=quantile(t_s$b_migrationRate.Anseriformes_to_Galliformes,0.025) )
# 
# s.Ne_G = data.frame(estimated=mean(t_s$popSize.t.Galliformes), 
#                        upper=quantile(t_s$popSize.t.Galliformes,0.975), lower=quantile(t_s$popSize.t.Galliformes,0.025) )
# s.Ne_A = data.frame(estimated=mean(t_s$popSize.t.Anseriformes), 
#                        upper=quantile(t_s$popSize.t.Anseriformes,0.975), lower=quantile(t_s$popSize.t.Anseriformes,0.025) )

hpd.s.rea_G = HPDinterval(as.mcmc(t_s$reassortmentRate.Galliformes))
hpd.s.rea_A = HPDinterval(as.mcmc(t_s$reassortmentRate.Anseriformes))
hpd.s.mig_GtA = HPDinterval(as.mcmc(t_s$b_migrationRate.Galliformes_to_Anseriformes))
hpd.s.mig_AtG = HPDinterval(as.mcmc(t_s$b_migrationRate.Anseriformes_to_Galliformes))
hpd.s.Ne_G = HPDinterval(as.mcmc(t_s$popSize.t.Galliformes))
hpd.s.Ne_A = HPDinterval(as.mcmc(t_s$popSize.t.Anseriformes))
hpd.s.cl = HPDinterval(as.mcmc(t_s$clockRate.c))

# get mascot rates
# s.mig_GtA = data.frame(estimated=mean(t_s$b_migrationRate.Galliformes_to_Anseriformes), 
#                        upper=quantile(t_s$b_migrationRate.Galliformes_to_Anseriformes,0.975), lower=quantile(t_s$b_migrationRate.Galliformes_to_Anseriformes,0.025) )
# s.mig_AtG = data.frame(estimated=mean(t_s$b_migrationRate.Anseriformes_to_Galliformes), 
#                        upper=quantile(t_s$b_migrationRate.Anseriformes_to_Galliformes,0.975), lower=quantile(t_s$b_migrationRate.Anseriformes_to_Galliformes,0.025) )
# 
# s.Ne_G = data.frame(estimated=mean(t_s$popSize.t.Galliformes), 
#                     upper=quantile(t_s$popSize.t.Galliformes,0.975), lower=quantile(t_s$popSize.t.Galliformes,0.025) )
# s.Ne_A = data.frame(estimated=mean(t_s$popSize.t.Anseriformes), 
#                     upper=quantile(t_s$popSize.t.Anseriformes,0.975), lower=quantile(t_s$popSize.t.Anseriformes,0.025) )

hpd.m.mig_GtA = HPDinterval(as.mcmc(t_m$b_migrationRate.Galliformes_to_Anseriformes))
hpd.m.mig_AtG = HPDinterval(as.mcmc(t_m$b_migrationRate.Anseriformes_to_Galliformes))
hpd.m.Ne_G = HPDinterval(as.mcmc(t_m$popSize.t.Galliformes))
hpd.m.Ne_A = HPDinterval(as.mcmc(t_m$popSize.t.Anseriformes))
hpd.m.cl = HPDinterval(as.mcmc(t_m$clockRate.c))

# rates = rbind(rates, data.frame(method="score", run=run, rea_A=t_s$reassortmentRate.Anseriformes, rea_G=t_s$reassortmentRate.Galliformes ,
#                                 mig_GtA=t_s$b_migrationRate.Galliformes_to_Anseriformes, mig_AtG=t_s$b_migrationRate.Anseriformes_to_Galliformes,
#                                 Ne_G=t_s$popSize.t.Galliformes, Ne_A=t_s$popSize.t.Anseriformes))

rea_rates <-  rbind(rea_rates,data.frame(run=run, state="Galliformes", 
                                         rea=t_s$reassortmentRate.Galliformes[which(
                                           t_s$reassortmentRate.Galliformes < hpd.s.rea_G[2]
                                           & t_s$reassortmentRate.Galliformes > hpd.s.rea_G[1])]))
rea_rates <- rbind(rea_rates, data.frame(run=run, state="Anseriformes", 
                                         rea=t_s$reassortmentRate.Anseriformes[which(
                                           t_s$reassortmentRate.Anseriformes < hpd.s.rea_A[2]
                                           & t_s$reassortmentRate.Anseriformes > hpd.s.rea_A[1])]))

g_mig_rates <- rbind(g_mig_rates, data.frame(method="SCoRe", run=run,
                                             mig=t_s$b_migrationRate.Galliformes_to_Anseriformes[which(
                                               t_s$b_migrationRate.Galliformes_to_Anseriformes < hpd.s.mig_GtA[2]
                                               & t_s$b_migrationRate.Galliformes_to_Anseriformes > hpd.s.mig_GtA[1])]))
g_mig_rates <- rbind(g_mig_rates, data.frame(method="MASCOT", run=run, 
                                             mig=t_m$b_migrationRate.Galliformes_to_Anseriformes[which(
                                                 t_m$b_migrationRate.Galliformes_to_Anseriformes < hpd.m.mig_GtA[2]
                                                 & t_m$b_migrationRate.Galliformes_to_Anseriformes > hpd.m.mig_GtA[1])]))

a_mig_rates <- rbind(a_mig_rates, data.frame(method="SCoRe", run=run,
                                             mig=t_s$b_migrationRate.Anseriformes_to_Galliformes[which(
                                               t_s$b_migrationRate.Anseriformes_to_Galliformes < hpd.s.mig_AtG[2]
                                               & t_s$b_migrationRate.Anseriformes_to_Galliformes > hpd.s.mig_AtG[1])]))
a_mig_rates <- rbind(a_mig_rates, data.frame(method="MASCOT", run=run,
                                             mig=t_m$b_migrationRate.Anseriformes_to_Galliformes[which(
                                               t_m$b_migrationRate.Anseriformes_to_Galliformes < hpd.m.mig_AtG[2]
                                               & t_m$b_migrationRate.Anseriformes_to_Galliformes > hpd.m.mig_AtG[1])]))

a_ne_rates <- rbind(a_ne_rates, data.frame(method="SCoRe", run=run,
                                           ne=t_s$popSize.t.Anseriformes[which(
                                             t_s$popSize.t.Anseriformes < hpd.s.Ne_A[2]
                                             & t_s$popSize.t.Anseriformes > hpd.s.Ne_A[1])]))
a_ne_rates <- rbind(a_ne_rates, data.frame(method="MASCOT", run=run,
                                           ne=t_m$popSize.t.Anseriformes[which(
                                             t_m$popSize.t.Anseriformes < hpd.m.Ne_A[2]
                                             & t_m$popSize.t.Anseriformes > hpd.m.Ne_A[1])]))

g_ne_rates <- rbind(g_ne_rates, data.frame(method="SCoRe", run=run,
                                           ne=t_s$popSize.t.Galliformes[which(
                                             t_s$popSize.t.Galliformes < hpd.s.Ne_G[2]
                                             & t_s$popSize.t.Galliformes > hpd.s.Ne_G[1])]))
g_ne_rates <- rbind(g_ne_rates, data.frame(method="MASCOT", run=run,
                                           ne=t_m$popSize.t.Galliformes[which(
                                             t_m$popSize.t.Galliformes < hpd.m.Ne_G[2]
                                             & t_m$popSize.t.Galliformes > hpd.m.Ne_G[1])]))

clock_rates <- rbind(clock_rates, data.frame(method="SCoRe", run=run,
                                             clock=t_s$clockRate.c[which(
                                               t_s$clockRate.c < hpd.s.cl[2]
                                               & t_s$clockRate.c > hpd.s.cl[1])]))
clock_rates <- rbind(clock_rates, data.frame(method="MASCOT", run=run,
                                             clock=t_m$clockRate.c[which(
                                               t_m$clockRate.c < hpd.m.cl[2]
                                               & t_m$clockRate.c > hpd.m.cl[1])]))
}


#blue, green
myColors <- c("#0072B2", "#009E73")
names(myColors) <- c("Galliformes","Anseriformes")
colScale <- scale_colour_manual(name = "state",values = myColors)
fillScale <- scale_fill_manual(name = "state",values = myColors)

myFills <- c("#009E73", "#67CEB2", "#0072B2", "#6CB1D8")
# names(myFills) <- c("Galliformes","Anseriformes")
fillScale2 <- scale_fill_manual(values = myFills)
colScale2 <- scale_colour_manual(values = myFills)


myFills2 <- c("#939393", "#B5B5B5")
fillScale3 <- scale_fill_manual(values = myFills2)



GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             # Original function by Jan Gleixner (@jan-glx)
                             # Adjustments by Wouter van der Bijl (@Axeman)
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
                               quantiles <- create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           }
)

create_quantile_segment_frame <- function(data, draw_quantiles, split = FALSE, grp = NULL) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y)
  ys <- ecdf(draw_quantiles)
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  violin.xs <- (stats::approxfun(data$y, data$x))(ys)
  if (grp %% 2 == 0) {
    data.frame(
      x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  } else {
    data.frame(
      x = ggplot2:::interleave(violin.xminvs, violin.xs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  }
}

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, 
        show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


##### reassortment plot

rea_rates$run <- factor(rea_rates$run)
p_rea <- ggplot(rea_rates, aes(x=run, y=rea, fill=state)) + geom_split_violin(scale="width", trim=T, draw_quantiles = 0.5, lwd=0.1) + 
  labs(x="subset", y="reassortment rate") +
  facet_grid(~run, scales = 'free_x', switch = 'x') + ylim(0, 0.3) +  theme(plot.margin = unit(c(2,0,0,0),"line")) +
  fillScale + theme_minimal() + theme(legend.title = element_blank(), legend.position = c(0.8, 0.95), legend.box = "horizontal", text = element_text(size=18), axis.text.x = element_blank(), axis.ticks.x = element_blank())#+ theme(plot.margin = unit(c(1,1,1,1),"cm"))
ggsave(paste0(out,"/rea_rates.pdf"), p_rea, width=9, height=5)

##### migration plot

a_mig_rates$state = "Ans -> Gal"
g_mig_rates$state = "Gal -> Ans"

mig_rates <- rbind(a_mig_rates, g_mig_rates)
mig_rates$key <- factor(paste(mig_rates$state, mig_rates$method))

p_mig <- ggplot(mig_rates, aes(x=interaction(state, run), y=mig, fill=key)) + geom_split_violin(scale="width", trim=T, draw_quantiles = 0.5, lwd=0.1) + 
  labs(x="subset", y="backward in time migration rate")+
    facet_grid(~run, scales = 'free_x', switch = 'x') +
  theme_minimal() + fillScale2 +
  theme(legend.title = element_blank(), legend.position = "bottom", text = element_text(size=18), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  guides(fill = guide_legend(nrow = 2))  +  theme(plot.margin = unit(c(0,0,2,0),"line"))
ggsave(paste0(out,"/p_mig.pdf"), p_mig, width=9, height=5)


##### Ne plot

a_ne_rates$state = "Anseriformes"
g_ne_rates$state = "Galliformes"

ne_rates <- rbind(a_ne_rates, g_ne_rates)
ne_rates$key <- factor(paste(ne_rates$state, ne_rates$method))

p_ne <- ggplot(ne_rates, aes(x=interaction(state, run), y=ne, fill=key)) + geom_split_violin(scale="width", trim=T, draw_quantiles = 0.5, lwd=0.1) + 
  labs(x="subset", y="effective population size")+
  facet_grid(~run, scales = 'free_x', switch = 'x') +
  theme_minimal() + fillScale2 +
  theme(legend.title = element_blank(), legend.position = "bottom", text = element_text(size=18), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  guides(fill = guide_legend(nrow = 2))  +  theme(plot.margin = unit(c(0,2,2,0),"line"))
ggsave(paste0(out,"/p_ne.pdf"), p_ne, width=9, height=5)


##### clock rates


clock_rates$run <- factor(clock_rates$run)
p_clock<- ggplot(clock_rates, aes(x=run, y=clock, fill=method)) + geom_split_violin(scale="width", trim=T, draw_quantiles = 0.5, lwd=0.1) + #geom_split_violin(scale="width", trim=T, draw_quantiles = 0.5, lwd=0.1) + 
  labs(x="subset", y="evolutionary rate per site and year")+
  facet_grid(~run, scales = 'free_x', switch = 'x') +
  theme_minimal() + fillScale3 +
  guides(fill = guide_legend(nrow = 1)) + ylim(0.0044,0.0056) +  theme(plot.margin = unit(c(2,2,0,0),"line")) +
  theme(legend.title = element_blank(), legend.position = c(0.5, 0.95), legend.box = "horizontal", text = element_text(size=18), axis.text.x = element_blank(), axis.ticks.x = element_blank()) #+ theme(plot.margin = unit(c(1,1,1,1),"cm"))

ggsave(paste0(out,"/p_clock.pdf"), p_clock, width=9, height=5)

p_all<- plot_grid(p_ne, p_mig, p_clock, p_rea, labels = c("A", "B","C", "D"),  align = "h", scale = 0.95,
          label_size = 30, nrow=2, ncol=2,axis="l", hjust=0.5, vjust = -0.3) +  theme(plot.margin = margin(2.0,1.0,0,0.8, "cm")) 
ggsave(paste0(out,"/p_all.pdf"), p_all)

