library(dplyr)
library(texreg)
library(fixest)
library(sp)
library(gstat)
library(rgdal)
library(ggplot2)
library(rgeos)
library(MatchIt)

setwd('C:/Users/Danie/OneDrive - Durham University/Daniel_phd_projects/CAA_II/')

data <- read.csv('data_final.csv')
carbon_data <- read.csv('./carbon/carbon_annual.csv') 
weather_data <- read.csv("weather_data/weather_data.csv")
disp_barry <- read.csv('disp_data.csv')
disp_neal <- read.csv('disp_data2.csv')
us <- readOGR('shapefiles/cb_2018_us_state_20m.shp')
us_cont <- us[!us$STUSPS %in% c('AK','PR','HI'), ]


# tests if cross-border sulfur is at least 1% of NAAQS
naaqs_test <- function(ExternalSulfur) {
  
  if (ExternalSulfur*1e+6 > 0.75 / 2.66) {
    
    return(1)
  }
  else {return(0)}
  
}


library("mice")

e_mids <- mice::mice(carbon_data, method = "cart", seed = 123)

carbon <- mice::complete(e_mids)

carbon <- carbon %>% group_by(id, year, state) %>%
  summarise(carbon = sum(carbon_short))

cair_states = c('TX', 'LA', 'MS', 'AL', 'GA', 'FL', 'SC', 'NC', 'VA', 'TN', 'KY', 'WV',
                'MO', 'IA','OH', 'IN', 'MI', 'PA', 'MD', 'DE', 'NJ', 'NY', 'IL', 'WI')

ozone_zone = c('AR', 'CT', 'MA')


dist_data <- data %>% 
  select(lat, lon, id, state) %>%
  distinct()

coordinates(dist_data) <- ~lon+lat

proj4string(dist_data) <- proj4string(us_cont)

dists <- as.data.frame(geosphere::dist2Line(p = dist_data, line = us_cont)) %>% 
  select(distance) %>%
  cbind(dist_data$id) %>%
  rename(id = 'dist_data$id')

data <- data %>% 
  left_join(carbon, by = c("id","year",'state')) %>%
  left_join(dists, by = 'id') %>%
  mutate(border_dist = distance/1000) %>%
  filter(op_time > 0)

data$drop <- 0
for (i in unique(data$id)) { 
  if (nrow(data[data$id == i, ]) != 24) {
    data[data$id == i, ]$drop <- 1
  }
}

data <- data[data$drop == 0, ] %>% select(-drop)

data_mids <- mice::mice(data, method = "pmm", seed = 123)

data_impute <- mice::complete(data_mids)


data_impute <- data_impute %>%
  mutate(ln_sulfur = log(sulfur),
         ln_carbon = log(carbon),
         did = treated*CAIR,
         cross_border = ifelse(ext_sulfur_sum > 0, 1, 0),
         naaqs_violate = 0,
         naaqs_domestic = 0,
         share = ext_sulfur_sum/total_sulfur*100)

for (i in 1:nrow(data_impute)) {
  # check if cross-border and within-state concentration exceed 1% NAAQS
  data_impute$naaqs_violate[i] = naaqs_test(data_impute$ext_sulfur_max[i])
  data_impute$naaqs_domestic[i] = naaqs_test(data_impute$total_sulfur[i]-data_impute$ext_sulfur_sum[i])
  
}

data <- data %>%
  mutate(ln_sulfur = log(sulfur),
         ln_carbon = log(carbon),
         did = treated*CAIR,
         cross_border = ifelse(ext_sulfur_mean > 0, 1, 0),
         did_cb = did*cross_border,
         share = ext_sulfur_sum/total_sulfur*100) %>%
  

# propensity score matching of 2004 groups
data_2004 <- data_impute %>% filter(year == 2004)

m_ps <- glm(formula = treated ~ heat_input + sulfur + generation + distance, family = binomial(), data = data_impute)

prs_df <- data.frame(pr_score = predict(m_ps, type = "response"), treated = m_ps$model$treated)

mod_match <- matchit(treated ~ heat_input + sulfur + generation + border_dist, method = "nearest", data = data_2004)
dta_m <- match.data(mod_match, distance = "nearest")

# filter matched IDs
data_matched <- data_impute[data_impute$id %in% dta_m$id, ]

# save data
write.csv(data, "data_noimpute.csv")
write.csv(data_impute, "data_impute.csv")

# Estimation

rm(list = ls())

library(plm)
library(lmtest)

data_impute <- read.csv('data_impute.csv')


# first stage regression

emit_lm <- feols(ln_sulfur ~ treated*CAIR +                 ## the key interaction: time × treatment status
                   heat_input + generation + op_time + sulfur_control + permits + allocated |  ## Other controls
                   id + year,                                                             ## FEs
                 cluster = ~id,                                                         ## Clustered SEs
                 data = data_impute)

cair_lm <- feols(ext_sulfur_mean ~ treated*CAIR +                      ## the key interaction: time × treatment status
                   sulfur + generation + heat_input + op_time + sulfur_control + permits + allocated |  ## Other controls
                   id + year,
                 cluster = ~id,                                             ## Clustered SEs
                 data = data_impute)

carbon_lm <- feols(ln_carbon ~ treated*CAIR +                                     ## the key interaction: time × treatment status
                   heat_input + generation + op_time + sulfur + sulfur_control + permits + allocated  |  ## Other controls
                   id + year,
                 cluster = ~id,                                             ## Clustered SEs
                 data = data_impute)

cair_lpm <- feols(naaqs_violate ~ treated*CAIR + 
                    sulfur + generation + heat_input + op_time + sulfur_control + permits + allocated |
                    id + year,
                    cluster = ~id, 
                    data = data_impute)

emit_ddd <- feols(ln_sulfur ~ treated*CAIR*cross_border 
                    + generation + heat_input + op_time + sulfur_control + permits + allocated |
                    id + year,
                  cluster = ~id, 
                  data = data_impute)

carbon_ddd <- feols(ln_carbon ~ treated*CAIR*cross_border +
                    generation + heat_input + op_time + sulfur_control + permits + allocated |
                    id + year,
                  cluster = ~id, 
                  data = data_impute)

ddd_sensitivity <- feols(ln_sulfur ~ treated*CAIR*cross_border + generation + heat_input + op_time + sulfur_control + permits + allocated |
                    id + year,
                  cluster = ~id, 
                  data = data_impute[data_impute$border_dist > 5/0.62, ])

ddd_sensitivity_carbon <- feols(ln_carbon ~ treated*CAIR*cross_border + 
                            sulfur + generation + heat_input + op_time + sulfur_control + permits + allocated |
                           id + year,
                         cluster = ~id, 
                         data = data_impute[data_impute$border_dist > 5/0.62, ])

naaqs_lpm <- feols(naaqs_domestic ~ treated*CAIR +                 ## the key interaction: time × treatment status
                   sulfur + heat_input + generation + op_time + sulfur_control + permits + allocated |  ## Other controls
                   id + year,                                                             ## FEs
                 cluster = ~id,                                                         ## Clustered SEs
                 data = data_impute)



etable(emit_lm, se = 'hetero')
etable(cair_lm, se = 'hetero')
etable(cair_lpm, se = 'hetero')
etable(naaqs_lpm, se = 'hetero')
etable(carbon_lm, se = 'hetero')
etable(emit_ddd, se = "hetero")
etable(carbon_ddd, se = "hetero")
etable(ddd_sensitivity, se = "hetero")
etable(ddd_sensitivity_carbon, se = "hetero")



ggplot(data_impute, aes(x = did, y = sulfur, color = as.factor(cross_border))) + geom_point()

# plots

# figure 1 (map of utilities)
al <- readOGR('shapefiles/cb_2018_01_bg_500k.shp')
ia <- readOGR('shapefiles/cb_2018_19_bg_500k.shp')
ne <- readOGR('shapefiles/cb_2018_31_bg_500k.shp')

us_cont_m <- spTransform(us_cont, CRSobj = "+proj=utm +zone=17 +datum=NAD83 +units=m")
al_m <- spTransform(al, CRSobj = "+proj=utm +zone=17 +datum=NAD83 +units=m")
ia_m <- spTransform(ia, CRSobj = "+proj=utm +zone=17 +datum=NAD83 +units=m")
ne_m <- spTransform(ne, CRSobj = "+proj=utm +zone=17 +datum=NAD83 +units=m")

al_df <- broom::tidy(al_m)
ia_df <- broom::tidy(ia_m)
ne_df <- broom::tidy(ne_m)

us_cont@data$program <- 0
us_cont@data[us_cont@data$STUSPS %in% cair_states, ]$program <- 1
us_cont@data[us_cont@data$STUSPS %in% ozone_zone, ]$program <- 2

us_cont@data <- us_cont@data %>% mutate(id = row.names(.))
us_cont_df <- broom::tidy(us_cont, region = "id") %>% 
  left_join(us_cont@data, by = c("id"="id"))

map_data <- data %>% 
  select(id, lat, lon, treated, ext_sulfur_mean, sulfur_rate, distance) %>%
  group_by(id, treated, lat, lon) %>%
  summarise(crossBorder = ifelse(mean(ext_sulfur_mean) > 0, 1, 0), distance = mean(distance), sulfur_rate = mean(sulfur_rate)) 

map_data$dummy <- 0
map_data[map_data$distance < 1000 & map_data$crossBorder == 0, ]$dummy <- 1
map_data[map_data$distance >= 1000 & map_data$crossBorder == 1, ]$dummy <- 2
map_data[map_data$distance < 1000 & map_data$crossBorder == 1, ]$dummy <- 3


map_data$sulfur_discrete <- 0
map_data[map_data$sulfur_rate >= 250 & map_data$sulfur_rate < 500,]$sulfur_discrete <- 1
map_data[map_data$sulfur_rate >= 500 & map_data$sulfur_rate < 750,]$sulfur_discrete <- 2
map_data[map_data$sulfur_rate >= 750,]$sulfur_discrete <- 3

stations <- weather_data %>% select(id, lat, lon) %>% distinct() %>% filter(lat < 50 & lat > 24)

figure0 <- ggplot() + 
  geom_polygon(data = us_cont_df, aes(x = long, y = lat, group = group, fill=factor(program)), colour = "black") +
  geom_point(data=map_data, aes(x=lon, y=lat),size=1.5, shape=22, fill="red", color='black') + 
  geom_point(data=stations, aes(x=lon, y=lat), size=1, shape=3, color='gold') +
  scale_fill_manual(values = c('steelblue4', 'steelblue2', "white"), name = 'Cap-and-trade program:', labels = c('No coverage', expression(paste(SO[2], sep="")), "Ozone")) +
  ggtitle("Coverage of US fossil-fired utilities under CAIR") +
  xlab('Longitude') +
  ylab('Latitude') +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightsteelblue1", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.position = c(0.16, 0.15),
        legend.box = "horizontal",
        legend.key= element_blank())

png(filename="figure0.png", 
    units="in", 
    width=7, 
    height=5.5, 
    pointsize=10, 
    res=600)
figure0
dev.off()


figure1 <- ggplot() + 
                geom_polygon(data = us_cont_df[us_cont_df$program==0 | us_cont_df$program==2,], aes(x = long, y = lat, group = group), colour = "black", fill="steelblue4") +
                geom_polygon(data = us_cont_df[us_cont_df$program==1,], aes(x = long, y = lat, group = group), colour = "black", fill="steelblue1") +
                geom_point(data=map_data, aes(x=lon, y=lat, shape=factor(dummy), fill=factor(sulfur_discrete)), size=1.5) + 
                scale_fill_manual(values = c('green', 'yellow', 'orange', 'red'), name = expression(paste(SO[2],' emit rate (g/s):')), labels = c('<250', '250-500', '500-750', '>750')) +
                scale_shape_manual(values = c(22,25,24,23), name = 'Geographic conditions\n(distance to state border / cross-border pollution):', 
                                   labels = c('Over 1km / No','Under 1km / No','Over 1km / Yes','Under 1km / Yes')) +
                xlab('Longitude') +
                ylab('Latitude') +
                xlim(-132, -70) +
                ggtitle(expression(paste("Average ", SO[2], " emission rates 1997-2020", sep=""))) +
                guides(fill = guide_legend(override.aes=list(shape=22, size=4), order=1),
                       shape = guide_legend(override.aes=list(shape=c(22,25,24,23), size=3, fill="steelblue4"), order=2)) +
                theme(axis.line = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill = "lightsteelblue1", colour = "black"),
                      panel.border = element_rect(colour = "black", fill=NA, size=2),
                      legend.position = c(0.18,0.25),
                      legend.background = element_rect(fill = "white", color = "black"),
                      legend.box = "vertical",
                      legend.key= element_blank())

png(filename="figure1.png", 
    units="in", 
    width=11, 
    height=7, 
    pointsize=10, 
    res=600)
figure1
dev.off()

library(RColorBrewer)

disp_barry <- disp_barry %>% 
  mutate(conc = conc*1e+6) %>%
  mutate(conc0 = cut(conc, breaks = c(0.0e+00,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1), include.lowest = TRUE))

pal <- c(NA, brewer.pal(n = 8, "YlOrRd"))

figure2 <- ggplot() + geom_polygon(data = al_df, aes(x = long, y = lat, group = group), colour = "gray60", fill='black') + 
  geom_tile(data = disp_barry, aes(x = lon, y = lat, fill=conc0), alpha=.6) + 
  coord_cartesian(ylim=c(3.445e+6,3.475e+6), xlim=c(-1.8e+5, -1.45e+5)) +
  scale_fill_manual(values = pal) +
  annotate(geom="text", y=3.47e+6, x=-1.5e+5, label="Alabama", col='white') +
  labs(title=expression(paste('Barry Power Plant: ',SO[2], " (",mu,"g/",m^{3},')', sep='')), x='', y='') +
  theme(title = element_text(size = 10),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightsteelblue1", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.key.height = unit(2, "cm"),
        legend.title = element_blank())

png(filename="figure2.png", 
    units="in", 
    width=5.5, 
    height=7, 
    pointsize=10, 
    res=600)
figure2
dev.off()

disp_neal <- disp_neal %>% 
  mutate(conc = conc*1e+6) %>%
  mutate(conc0 = cut(conc, breaks = c(0.0e+00,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e+1), include.lowest = TRUE)) %>%
  filter(!is.na(conc0))

figure2b <- ggplot() + geom_polygon(data = ia_df, aes(x = long, y = lat, group = group), colour = "gray60", fill='black') +
  geom_polygon(data = ne_df, aes(x = long, y = lat, group = group), colour = "gray60", fill='gray20') +
  geom_tile(data = disp_neal, aes(x = lon, y = lat, fill = conc0), alpha=.6) + 
  coord_cartesian(ylim=c(4.79e+6,4.823e+6), xlim=c(-7.9e+5, -7.5e+5)) +
  scale_fill_manual(values = c(pal,"violetred4")) +
  labs(title=expression(paste('Neal South Plant: ',SO[2], " (",mu,"g/",m^{3},')', sep='')), x='', y='') +
  annotate(geom='text', x=-7.55e+5, y=4.81e+6, label='Iowa', col='white') +
  annotate(geom='text', x=-7.85e+5, y=4.795e+6, label='Nebraska', col='white') +
  theme(title = element_text(size = 10),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightsteelblue1", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.key.height = unit(1.8, "cm"),
        legend.title = element_blank())

png(filename="figure2b.png", 
    units="in", 
    width=5.5, 
    height=7, 
    pointsize=10, 
    res=600)
figure2b
dev.off()



# figure 2 (timeline plots)

line_data_ <- data %>% filter(total_sulfur > 0) %>%
             group_by(id, year) %>%
             summarise(sulfur = sum(sulfur), ext = sum(ext_sulfur_sum), share = mean(share), treated = max(treated)) %>%
             group_by(year, treated) %>% 
             summarise(sulfur = mean(sulfur), ext = mean(ext), share=mean(share))

# sulfur dioxide emissions
figure3 <- ggplot(data = line_data_, aes(x = year, y = sulfur)) +
           geom_line(aes(color=as.factor(treated)), size = 1) + 
           scale_color_manual(values = c("steelblue4", "red"), name = "Program:", labels = c("Not covered by CAIR", "Covered by CAIR")) +
           ylab(expression(paste("Average annual  ",SO[2]," emissions (tonnes)", sep=""))) +
           theme(axis.line = element_blank(),
                 axis.title.x = element_blank(),
                 axis.title.y = element_text(size=10),
                 title = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_rect(fill = "lightsteelblue1", colour = "black"),
                 panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                 legend.background = element_rect(fill = "white", colour = "black"),
                 legend.position = c(0.61, 0.87),
                 legend.box = "vertical",
                 legend.key = element_blank()) 

# cross-state externalities
figure3b <- ggplot(data = line_data_, aes(x = year, y = ext*1e+6)) +
            geom_line(aes(color=as.factor(treated)), size = 1) + 
            scale_color_manual(values = c("steelblue4", "red"), name = "Program:", labels = c("Not covered by CAIR", "Covered by CAIR")) +
            ylab(expression(paste("Cross-border ", SO[2], " (", mu, "gs/", m^{3}, ")", sep=""))) +
            theme(axis.line = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size = 10),
                  title = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "lightsteelblue1", colour = "black"),
                  panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                  legend.position = 'none')

figure3c <- ggplot(data = line_data_, aes(x = year, y = share*100)) +
  geom_line(aes(color=as.factor(treated)), size = 1) + 
  scale_color_manual(values = c("steelblue4", "red"), name = "Program:", labels = c("Not covered by CAIR", "Covered by CAIR")) +
  ylab(expression(paste("Cross-border ", SO[2], " as share of total (%)", sep=""))) +
  ylim(0, 25) +
  theme(axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightsteelblue1", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        legend.position = 'none',)


png(filename="figure3a.png", 
    units="in", 
    width=3.5, 
    height=5, 
    pointsize=10, 
    res=600)
figure3
dev.off()

png(filename="figure3b.png", 
    units="in", 
    width=3.5, 
    height=5, 
    pointsize=10, 
    res=600)
figure3b
dev.off()

png(filename="figure3c.png", 
    units="in", 
    width=3.5, 
    height=5, 
    pointsize=10, 
    res=600)
figure3c
dev.off()

# figure 4 (permit demand 
market_data <- data %>% select(id, year, treated, sulfur, permits, allocated, clearing_price) %>%
  group_by(id, year) %>%
  summarize(treated = max(treated), sulfur = sum(sulfur), permits = sum(permits), allocated = sum(allocated), price = mean(clearing_price)) %>%
  group_by(year, treated) %>% 
  summarize(permits = mean(permits), 
            sulfur=mean(sulfur),
            allocated = mean(allocated),
            price = mean(price)) 

figure4 <-  ggplot(data = market_data, aes(x = year)) +
  #geom_bar(data = market_data[market_data$treated == 1, ], aes(y = permits), stat='identity', size = 1, alpha=0.3, fill='red') +
  #geom_bar(data = market_data[market_data$treated == 0, ], aes(y = permits), stat='identity', size = 1, alpha=.5, fill = 'steelblue4') +
  geom_line(aes(y = sulfur, color = factor(treated)), size=1) +
  geom_line(aes(y = price*100), linetype='dashed', color = 'black', size=1) +
  scale_color_manual(values = c("steelblue4", "red"), name = "Program:", labels = c("Not covered by CAIR", "Covered by CAIR")) +
  scale_y_continuous(sec.axis = sec_axis(trans=~./100, name="Permit price ($/tonne)")) +
  ylab(expression(paste(SO[2], " emissions (tonnes)", sep=""))) +
  ggtitle("Average emissions and allowances by group (1997 - 2020)") +
  theme(axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 11),
        title = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightsteelblue1", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        legend.background = element_rect(fill = "lightsteelblue1", colour = "grey50", size = 1),
        legend.position = 'bottom',
        legend.box = "vertical",
        legend.key= element_blank()) 

png(filename="figure4.png", 
    units="in", 
    width=6, 
    height=4, 
    pointsize=10, 
    res=600)
figure4
dev.off()


# figure 5 (event study plots)

data$time_to_cair <- data$year - 2005
data_impute$time_to_cair <- data_impute$year - 2005

data$id <- as.factor(data$id)
data$state <- as.factor(data$state)

ext_twfe <- feols(ext_sulfur_sum ~ i(time_to_cair, treated, ref = -1) +                           ## the key interaction: time × treatment status
                  sulfur_rate + sulfur_control + op_time + heat_input + allocated + permits   |   ## Other controls
                  id + year,                                                                              ## FEs
                  cluster = ~id,                                                                          ## Clustered SEs
                  data = data_impute)

ext_data <- data_impute %>% 
              dplyr::select(time_to_cair) %>% 
              distinct() %>%
              cbind(c(ext_twfe$coefficients[1:7], 0, ext_twfe$coefficients[8:23]),
                    c(ext_twfe$se[1:7]*1.96, 0, ext_twfe$se[8:23]*1.96))

naaqs_twfe <- feols(naaqs_violate ~ i(time_to_cair, treated, ref = -1) +                           ## the key interaction: time × treatment status
                    sulfur + sulfur_control + generation + op_time + heat_input + allocated + permits   |   ## Other controls
                    id + year,                                                                              ## FEs
                  cluster = ~id,                                                                          ## Clustered SEs
                  data = data_mathed)

naaqs_data <- data_mathed %>% 
  dplyr::select(time_to_cair) %>% 
  distinct() %>%
  cbind(c(naaqs_twfe$coefficients[1:7], 0, naaqs_twfe$coefficients[8:23]),
        c(naaqs_twfe$se[1:7]*1.96, 0, naaqs_twfe$se[8:23]*1.96))


emit_twfe <- feols(log(sulfur) ~ i(time_to_cair, treated, ref = -1) +                          ## the key interaction: time × treatment status
                      sulfur_control + op_time + heat_input + allocated + permits  |    ## Other controls
                      id + year,                                                                 ## FEs
                      cluster = ~id,                                                             ## Clustered SEs
                      data = data_impute)

emit_data <- data_impute %>% dplyr::select(time_to_cair) %>% 
  distinct() %>%
  cbind(c(emit_twfe$coefficients[1:7], 0, emit_twfe$coefficients[8:23]),
        c(emit_twfe$se[1:7]*1.96, 0, emit_twfe$se[8:23]*1.96))

carbon_twfe <- feols(log(carbon) ~ i(time_to_cair, treated, ref = -1) +                               ## the key interaction: time × treatment status
                     sulfur_control + op_time + heat_input + sulfur + generation + allocated + permits |  ## Other controls
                     id + year,                                                                              ## FEs
                   cluster = ~id,                                                                          ## Clustered SEs
                   data = data_impute)

carbon_data <- data_impute %>% dplyr::select(time_to_cair) %>% 
  distinct() %>%
  cbind(c(carbon_twfe$coefficients[1:7], 0, carbon_twfe$coefficients[8:23]),
        c(carbon_twfe$se[1:7]*1.96, 0, carbon_twfe$se[8:23]*1.96))

names(ext_data) <- c('time_to_treat', 'coef', 'CI')
names(naaqs_data) <- c('time_to_treat', 'coef', 'CI')
names(emit_data) <- c('time_to_treat', 'coef', 'CI')
names(carbon_data) <- c('time_to_treat', 'coef', 'CI')


figure5a <- ggplot(data = ext_data, aes(x = time_to_treat)) +
  geom_ribbon(aes(ymin=coef-CI, ymax=coef+CI), fill="steelblue4", alpha=0.3) +
  geom_line(aes(y = coef), color = 'red', linetype='dashed') +
  geom_vline(xintercept =  0, linetype = "dashed") +
  ylab("Event study coefficients and 95% confidence interval") +
  xlab("Years since CAIR announced in 2005") +
  ggtitle(expression(paste("Cross-border ",SO[2]," log mgs/", m^{3}, sep=""))) +
  theme(axis.line = element_blank(),
        plot.margin=unit(c(1,1,1,1), 'cm'),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        title = element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightsteelblue1", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5))

png(filename="figure5a.png", 
    units="in", 
    width=3.5, 
    height=5, 
    pointsize=10, 
    res=600)
figure5a
dev.off()


figure5b <- ggplot(data = emit_data, aes(x = time_to_treat)) +
  geom_ribbon(aes(ymin=coef-CI, ymax=coef+CI), fill="steelblue4", alpha=0.3) +
  geom_line(aes(y = coef), color = 'red', linetype='dashed') +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Event study coefficients and 95% confidence interval") +
  xlab("Years since CAIR announced in 2005") +
  ggtitle(expression(paste(SO[2]," dispersion: log mgs/", m^{3}, sep=""))) +
  theme(axis.line = element_blank(),
        plot.margin=unit(c(1,1,1,1), 'cm'),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        title = element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightsteelblue1", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5))

png(filename="figure5b.png", 
    units="in", 
    width=3.5, 
    height=5, 
    pointsize=10, 
    res=600)
figure5b
dev.off()


figure5c <- ggplot(data = carbon_data, aes(x = time_to_treat)) +
  geom_ribbon(aes(ymin=coef-CI, ymax=coef+CI), fill="steelblue4", alpha=0.3) +
  geom_line(aes(y = coef), color = 'red', linetype='dashed') +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Event study coefficients and 95% confidence interval") +
  xlab("Years since CAIR announced in 2005") +
  ggtitle(expression(paste("Annual ",CO[2]," emissions: log tonnes", sep=""))) +
  theme(axis.line = element_blank(),
        plot.margin=unit(c(1,1,1,1), 'cm'),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        title = element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightsteelblue1", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5))

png(filename="figure5c.png", 
    units="in", 
    width=3.5, 
    height=5, 
    pointsize=10, 
    res=600)
figure5c
dev.off()

figure5d <- ggplot(data = naaqs_data, aes(x = time_to_treat)) +
  geom_ribbon(aes(ymin=coef-CI, ymax=coef+CI), fill="steelblue4", alpha=0.3) +
  geom_line(aes(y = coef), color = 'red', linetype='dashed') +
  geom_vline(xintercept =  0, linetype = "dashed") +
  ylab("Event study coefficients and 95% confidence interval") +
  xlab("Years since CAIR announced in 2005") +
  ggtitle(expression(paste("Cross-border ",SO[2],"> 1% NAAQS", sep=""))) +
  theme(axis.line = element_blank(),
        plot.margin=unit(c(1,1,1,1), 'cm'),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        title = element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightsteelblue1", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5))

png(filename="figure5d.png", 
    units="in", 
    width=7, 
    height=4, 
    pointsize=10, 
    res=600)
figure5d
dev.off()


# figure 6

bar_cross = data_impute %>% 
  select(year, sulfur, treated, cross_border) %>%
  group_by(year, treated, cross_border) %>%
  summarise(sulfur = mean(sulfur)) %>%
  filter(cross_border == 1)%>%
  select(-cross_border) %>%
  rename(sulfur_c1 = sulfur)

bar_data = data_impute %>% 
  select(year, sulfur, treated, cross_border) %>%
  group_by(year, treated, cross_border) %>%
  summarise(sulfur = mean(sulfur)) %>%
  filter(cross_border == 0) %>%
  select(-cross_border) %>%
  rename(sulfur_c0 = sulfur) %>%
  left_join(bar_cross, by = c("year","treated")) %>%
  mutate(sulfur_c1_minus_c0 = sulfur_c1 - sulfur_c0)

figure6 <- ggplot(bar_data, aes(x = year, y = sulfur_c1_minus_c0)) +
  geom_bar(aes(fill = as.factor(treated)), position = "dodge", stat = "identity") +
  scale_fill_manual(values = c("steelblue4", "red"), name = "Program:", labels = c("Not covered by CAIR", "Covered by CAIR")) +
  geom_vline(xintercept = 2005, linetype = "dashed", size = 1) +
  ylab("Average annual sulfur (tonnes)") +
  ggtitle(expression(paste("Difference in ", SO[2]," emissions between cross-border polluting states and others", sep=""))) +
  theme(axis.line = element_blank(),
        plot.margin=unit(c(1,1,1,1), 'cm'),
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank(),
        title = element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightsteelblue1", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        legend.background = element_rect(fill = "white", colour = "black"),
        legend.position = c(0.85, 0.85),
        legend.box = "horizontal")

png(filename="figure6.png", 
    units="in", 
    width=7.8, 
    height=5, 
    pointsize=10, 
    res=600)
figure6
dev.off()




bar_10km = data_impute %>% 
  select(year, ext_sulfur_mean, treated, border_dist) %>%
  mutate(border_dist = ifelse(border_dist < 10, 1, 0)) %>%
  group_by(year, treated, border_dist) %>%
  summarise(sulfur = mean(ext_sulfur_mean)) %>%
  filter(border_dist == 1)%>%
  select(-border_dist) %>%
  rename(sulfur_under10km = sulfur)

bar_10km_data = data_impute %>% 
  select(year, ext_sulfur_mean, treated, border_dist) %>%
  mutate(border_dist = ifelse(border_dist < 10, 1, 0)) %>%
  group_by(year, treated, border_dist) %>%
  summarise(sulfur = mean(ext_sulfur_mean)) %>%
  filter(border_dist == 0) %>%
  select(-border_dist) %>%
  rename(sulfur_over10km = sulfur) %>%
  left_join(bar_10km, by = c("year","treated")) %>%
  mutate(sulfur_10km_diff = sulfur_under10km - sulfur_over10km)

figure7 <- ggplot(bar_10km_data, aes(x = year, y = sulfur_10km_diff)) +
  geom_bar(aes(fill = as.factor(treated)), position = "dodge", stat = "identity") +
  scale_fill_manual(values = c("steelblue4", "red"), name = "Program:", labels = c("Not covered by CAIR", "Covered by CAIR")) +
  geom_vline(xintercept = 2005, linetype = "dashed", size = 1) +
  ylab("Average annual sulfur (tonnes)") +
  ggtitle(expression(paste("Difference in cross-border ", SO[2]," between plants within 10km from border and others", sep=""))) +
  theme(axis.line = element_blank(),
        plot.margin=unit(c(1,1,1,1), 'cm'),
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank(),
        title = element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightsteelblue1", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        legend.background = element_rect(fill = "white", colour = "black"),
        legend.position = c(0.85, 0.85),
        legend.box = "horizontal")

png(filename="figure7.png", 
    units="in", 
    width=7.8, 
    height=5, 
    pointsize=10, 
    res=600)
figure7
dev.off()


