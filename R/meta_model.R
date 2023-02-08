library(tornado)
library(modelr)
parms_disprog <- c('ontx.u5','ontx.o5','hivartOR:mn','hivartOR:sg','notx.u5','notx.o5','notxH.u5',
                   'notxH.o5','notxHA.u5','notxHA.o5','hivpi','artp','LTBI04','LTBI514','LTBI04hi',
                   'LTBI514hi','tptRRtstpos','tptRRhivpos','prog04','prog514')

parms_epicalc <- c('cfr.notx','cfr.tx', 'p.tbdx.1yr', 'tptRR', 'CDR', 'CDRi')
parms_inteff <- c('tpt.resultOR','tpt.initiationOR','tpt.completionOR','tb.resultOR','tb.diagnosisOR')

parms_cstssoc <- names(parmsout)[grepl('c.soc', names(parmsout))]
parms_cstsint <- names(parmsout)[grepl('c.int', names(parmsout))]

lambda <- 5e3
psaout[,c('lys.soc','lys.int'):=.(1-lyl.soc, 1-lyl.int)]
psaout[,c('Dlys','Dcost'):=.(lys.int-lys.soc, cost.int-cost.soc)]
psaout[,c('iNMB','iNHB'):=.(Dlys*lambda-Dcost,Dlys-(Dcost/lambda))]

psaout[,.(lys.soc, lys.int, Dlys, Dcost, iNHB, iNMB)]
psaout[,sum(iNMB<1)]
# using the 
library(dampack)
dim(psaout)
dim(parmsout)

parms <- names(parmsout)
parms <- parms[!parms %in% c('id','iso3', 'age')]

psa_obj <- list(make_psa_obj(
  cost=psaout[iso3=='CMR',.(cost.soc, cost.int)],
  effectiveness=psaout[iso3=='CMR',.(lys.soc, lys.int)],
  parameters = parmsout[iso3=='CMR',..parms],
  strategies = c('Control', 'Intervention'),
  currency = "$",
  other_outcome = NULL
),
make_psa_obj(
  cost=psaout[iso3=='UGA',.(cost.soc, cost.int)],
  effectiveness=psaout[iso3=='UGA',.(lys.soc, lys.int)],
  parameters = parmsout[iso3=='UGA',..parms],
  strategies = c('Control', 'Intervention'),
  currency = "$",
  other_outcome = NULL
))

# str(psa_obj)

head(psa_obj[[1]]$cost); head(psa_obj[[2]]$cost)
head(psa_obj[[1]]$effectiveness); head(psa_obj[[2]]$effectiveness)

plot1 <- plot(psa_obj[[1]]) + 
  xlab('DALYs') + ylab('Cost (US$)') +
  scale_y_continuous(limits = c(-250, 700)) +
  geom_vline(xintercept=0, size = 0.5) + geom_hline(yintercept=0, size = 0.5) + 
  theme(legend.position="top", legend.title = element_blank()); 
plot2 <- plot(psa_obj[[2]]) + 
  xlab('DALYs') + ylab('Cost (US$)') +
  scale_y_continuous(limits = c(-250, 700)) +
  geom_vline(xintercept=0) + geom_hline(yintercept=0) + 
  theme(legend.position="top", legend.title = element_blank()); 
combined_plot <- ggpubr::ggarrange(plot1, plot2 + xlab(''), ncol = 2)
combined_plot

ggsave(combined_plot,file=gh('plots/CEplot') + SA + '.png',w=7,h=5)

# sensitivity analysis
vars <- c(parms_epicalc, parms_inteff, parms_cstsint) 
cmr_onmb <- owsa(psa_obj[[1]], params = vars, outcome = c("eff"),strategies = 'Intervention', wtp = lambda)
uga_onmb <- owsa(psa_obj[[2]], params = vars, outcome = c("eff"),strategies = 'Intervention', wtp = lambda)

# plot(cmr_onmb,
#      n_x_ticks = 5)
# plot(uga_onmb,
#      n_x_ticks = 5)

cmr_ponmb <- owsa_tornado(cmr_onmb) + 
  ylab('Net monetary benefit') +
  theme(legend.position="top", legend.title = element_blank()); 

uga_ponmb <- owsa_tornado(uga_onmb) +
  ylab('Net monetary benefit') +
  theme(legend.position="top", legend.title = element_blank()); 

# Combine
ggpubr::ggarrange(cmr_ponmb, uga_ponmb, ncol = 2)

# custom de novo regression models 
psa_dat <- merge(parmsout,
                 psaout,
                 by = c('id', 'iso3', 'age'))

# library(DataExplorer)
# plot_missing(psa_dat)
# plot_histogram(psa_dat)
# 
# ggplot(psa_dat, aes(x=Dcost, y=iNMB, col=iso3)) +
#   geom_point()
# 
# library(psych)
# pairs.panels(psa_dat)
             
by_country <- psa_dat %>% 
  group_by(iso3) %>% 
  nest()

by_country
by_country$data[[1]]

# explanatory_vars <- c(parms_epicalc, parms_inteff, parms_cstssoc, parms_cstsint)
vars <- c(parms_epicalc, parms_inteff, parms_cstssoc, parms_cstsint)

# Identify the parameters that vary
paramvar <- parmsout[,..vars] %>% summarise_all(function(x) var(x))
explanatory_vars <- names(paramvar)[which(paramvar!=0)]

# inhb
## run a regression with all variables of interest 
fmla.unstd <- explanatory_vars %>%  ## begin with vector of explanatory column names
  str_c(collapse = "+") %>%     ## combine all names of the variables of interest separated by a plus
  str_c("iNHB ~ ", .) 

# fmla.unstd <- as.formula(paste0(YY,"~",paste0(explanatory_vars,collapse="+")))
fmla.std <- as.formula(paste0('iNHB',"~",paste0(paste0("scale(",explanatory_vars,")"),collapse="+")))

country_model <- function(df) {
  lm(fmla.unstd, data = df)
}

inhb <- by_country %>% 
  mutate(model = map(data, country_model))

inhb <- inhb %>% 
  mutate(
    resids = map2(data, model, add_residuals)
  )
inhb

resids <- unnest(inhb, resids)
resids

resids %>% 
  ggplot(aes(id, resid)) +
  geom_line(aes(group = age), alpha = 1 / 3) + 
  geom_smooth(se = FALSE)
#> `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

resids %>% 
  ggplot(aes(id, resid)) +
  geom_line(alpha = 1 / 3) + 
  facet_grid(iso3~age)

# Model quality
glance <- inhb %>% 
  mutate(glance = map(model, broom::glance)) %>% 
  unnest(glance, .drop = TRUE)
glance

glance %>% 
  arrange(r.squared)

glance %>% 
  ggplot(aes(iso3, r.squared)) + 
  geom_jitter(width = 0.5)

# using tornado package
# iNHB
torn1 <- inhb %>%
  pull(model) %>%
  pluck(1) %>% # Extract the coefficients and standard errors for each model
  tornado::tornado(type = "ranges", alpha = 0.05/2)
torn2 <- inhb %>%
  pull(model) %>%
  pluck(2) %>% 
  tornado::tornado(type = "ranges", alpha = 0.05/2)

# plotting
torn1_plot <- plot(torn1, geom_bar_control = list(width = 0.6),
                   sensitivity_colors = c("#66C2A5","#FC8D62"),
                   geom_point_control = list(size = 3, fill = "purple", col = "purple")) + 
  xlab('Parameter') + ylab('Incremental net health benefit') +
  theme(legend.position="top", legend.title = element_blank()); 

torn2_plot <- plot(torn2, geom_bar_control = list(width = 0.6),
                   sensitivity_colors = c("#66C2A5","#FC8D62"),
                   geom_point_control = list(size = 3, fill = "purple", col = "purple")) + 
  xlab('Parameter') + ylab('Incremental net health benefit') +
  theme(legend.position="top", legend.title = element_blank()); 

# Combine
ggpubr::ggarrange(torn1_plot, torn2_plot, ncol = 2)
ggsave(tornado_DLYSs,file=gh('plots/mwsa_inhb') + SA + '.png',w=7,h=5)

# inmb
## run a regression with all variables of interest 
fmla.unstd <- explanatory_vars %>%  ## begin with vector of explanatory column names
  str_c(collapse = "+") %>%     ## combine all names of the variables of interest separated by a plus
  str_c("iNMB ~ ", .) 

# fmla.unstd <- as.formula(paste0(YY,"~",paste0(explanatory_vars,collapse="+")))
fmla.std <- as.formula(paste0('iNMB',"~",paste0(paste0("scale(",explanatory_vars,")"),collapse="+")))

country_model <- function(df) {
  lm(fmla.unstd, data = df)
}

inmb <- by_country %>% 
  mutate(model = map(data, country_model))

torn1 <- inmb %>%
  pull(model) %>%
  pluck(1) %>% # Extract the coefficients and standard errors for each model
  tornado::tornado(type = "ranges", alpha = 0.05/2)
torn2 <- inmb %>%
  pull(model) %>%
  pluck(2) %>% 
  tornado::tornado(type = "ranges", alpha = 0.05/2)

# plotting
torn1_plot <- plot(torn1, geom_bar_control = list(width = 0.6),
                   sensitivity_colors = c("#66C2A5","#FC8D62"),
                   geom_point_control = list(size = 3, fill = "purple", col = "purple")) + 
  xlab('Parameter') + ylab('Incremental net monetary benefit') +
  theme(legend.position="top", legend.title = element_blank()); 

torn2_plot <- plot(torn2, geom_bar_control = list(width = 0.6),
                   sensitivity_colors = c("#66C2A5","#FC8D62"),
                   geom_point_control = list(size = 3, fill = "purple", col = "purple")) + 
  xlab('Parameter') + ylab('Incremental net monetary benefit') +
  theme(legend.position="top", legend.title = element_blank()); 

# Combine
ggpubr::ggarrange(torn1_plot, torn2_plot, ncol = 2)
ggsave(tornado_DLYSs,file=gh('plots/mwsa_inmb') + SA + '.png',w=7,h=5)

# cost
# explanatory_vars <- c(parms_epicalc[parms_epicalc!=c('CDRi')], parms_inteff, parms_cstssoc, parms_cstsint)
potential_vars <- c(parms_inteff, parms_cstssoc, parms_cstsint)
paramvar <- parmsout[,..potential_vars] %>% summarise_all(function(x) var(x))
explanatory_vars <- names(paramvar)[which(paramvar!=0)]

fmla.unstd <- as.formula(paste0('Dcost',"~",paste0(explanatory_vars,collapse="+")))
fmla.std <- as.formula(paste0('Dcost',"~",paste0(paste0("scale(",explanatory_vars,")"),collapse="+")))

cost_model <- function(df) {
  lm(fmla.unstd, data = df)
}

# cost model
dcost <- by_country %>% 
  mutate(model = map(data, cost_model))

# using tornado package
torn1 <- dcost %>%
  pull(model) %>%
  pluck(1) %>% tornado::tornado(type = "ranges", alpha = 0.05/2)
torn2 <- dcost %>%
  pull(model) %>%
  pluck(2) %>% tornado::tornado(type = "ranges", alpha = 0.05/2)

# plotting
torn1_plot <- plot(torn1, geom_bar_control = list(width = 0.6),
                   sensitivity_colors = c("#66C2A5","#FC8D62"),
                   geom_point_control = list(size = 3, fill = "purple", col = "purple")) + 
  xlab('Parameter') + ylab('Incremental costs') +
  theme(legend.position="top", legend.title = element_blank()); 

torn2_plot <- plot(torn2, geom_bar_control = list(width = 0.6),
                   sensitivity_colors = c("#66C2A5","#FC8D62"),
                   geom_point_control = list(size = 3, fill = "purple", col = "purple")) + 
  xlab('Parameter') + ylab('Incremental costs') +
  theme(legend.position="top", legend.title = element_blank()); 

# Combine
tornado_dcost <- ggpubr::ggarrange(torn1_plot, torn2_plot, ncol = 2)
tornado_dcost
ggsave(tornado_dcost,file=gh('plots/mwsa_costs') + SA + '.png',w=7,h=5)

# effects
# explanatory_vars <- c(parms_disprog,parms_epicalc[parms_epicalc!=c('CDRi')], parms_inteff, parms_cstssoc, parms_cstsint)
potential_vars <- c(parms_epicalc,parms_inteff)
paramvar <- parmsout[,..potential_vars] %>% summarise_all(function(x) var(x))
explanatory_vars <- names(paramvar)[which(paramvar!=0)]

fmla.unstd <- as.formula(paste0('Dlys',"~",paste0(explanatory_vars,collapse="+")))
fmla.std <- as.formula(paste0('Dlys',"~",paste0(paste0("scale(",explanatory_vars,")"),collapse="+")))

cost_model <- function(df) {
  lm(fmla.unstd, data = df)
}

# cost model
DLYS <- by_country %>% 
  mutate(model = map(data, cost_model))

# using tornado package
torn1 <- DLYS %>%
  pull(model) %>%
  pluck(1) %>% tornado::tornado(type = "ranges", alpha = 0.05/2)
torn2 <- DLYS %>%
  pull(model) %>%
  pluck(2) %>% tornado::tornado(type = "ranges", alpha = 0.05/2)

# plotting
torn1_plot <- plot(torn1, geom_bar_control = list(width = 0.4),
                   sensitivity_colors = c("#66C2A5","#FC8D62"),
                   geom_point_control = list(size = 3, fill = "purple", col = "purple")) + 
  xlab('Parameter') + ylab('Life years saved') +
  theme(legend.position="top", legend.title = element_blank()); 

torn2_plot <- plot(torn2, geom_bar_control = list(width = 0.4),
                   sensitivity_colors = c("#66C2A5","#FC8D62"),
                   geom_point_control = list(size = 3, fill = "purple", col = "purple")) + 
  xlab('Parameter') + ylab('Life years saved') +
  theme(legend.position="top", legend.title = element_blank()); 

# Combine
# ggpubr::ggarrange(cmr_p, uga_p + rremove("ylab") + rremove("y.text"), ncol = 2)
tornado_DLYSs <- ggpubr::ggarrange(torn1_plot, torn2_plot, ncol = 2)
tornado_DLYSs
ggsave(tornado_DLYSs,file=gh('plots/mwsa_effects') + SA + '.png',w=7,h=5)
