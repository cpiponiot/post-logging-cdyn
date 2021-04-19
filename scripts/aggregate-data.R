#### Open and aggregate Paracou data #### 

library(data.table)
library(BIOMASS)
sapply(list.files("R", full.names = TRUE), source)

paracou <- read.csv2("data/paracou.csv", stringsAsFactors = F)
setDT(paracou)

paracou$idtree <-
  paste("prc",
        paracou$n_parcelle,
        paracou$n_carre,
        paracou$i_arbre,
        sep = "_")
paracou$dbh <- paracou$circonf / pi
paracou <-
  subset(paracou, !campagne %in% c(1996, 1998, 2000, 2002, 2008, 2010, 2012))
paracou <- subset(paracou, n_parcelle %in% 1:12)

tab_tree = data.table(
  plot = paracou$n_parcelle,
  subplot = paracou$n_carre,
  idtree = paracou$idtree,
  name = paste(paracou$Genre, paracou$espece),
  logged = (paracou$code_mesure == 4 &
              paracou$code_vivant == 0),
  vern = paracou$nomPilote,
  t = paracou$campagne,
  xplot = paracou$x, yplot = paracou$y
)
tab_tree = tab_tree[, .(name=last(name), 
                        logged=logged[t == 1987], 
                        vern=last(vern), 
                        xplot=last(xplot),
                        yplot=last(yplot)), 
                    .(plot, subplot, idtree)]
tab_tree$logged = tab_tree$logged & !is.na(tab_tree$logged)

paracou = subset(paracou, code_vivant == 1)

data = data.table(
  idtree = paracou$idtree,
  dbh = paracou$dbh,
  year = paracou$campagne,
  ladder = paracou$code_vivant * (paracou$code_mesure %in% c(1:4, 10)),
  minDBH = 10
)
data = merge(data, tab_tree, by = 'idtree')

data = subset(data, !is.na(dbh))

plot_data_prc = data.table(
  plot = 1:12,
  plot.size = 6.25,
  treat = c(0:3, 2, 0, 1, 3, 1, 2, 0, 3)
)

# paracou$n_parcelle[grep("17", paracou$n_parcelle)] <- "17"
# paracou$n_parcelle[grep("18", paracou$n_parcelle)] <- "18"
# 
# subplot_prc = unique(data.frame(plot = as.numeric(as.character(paracou$n_parcelle)),
#                                 subplot = as.numeric(as.character(paracou$n_carre))))
# plot_data_prc = merge(plot_data_prc, subplot_prc, by="plot")
# plot_data_prc$treat[plot_data_prc$treat == 0] = "ctrl"
# plot_data_prc$treat[plot_data_prc$treat == 1] = "CL"
# plot_data_prc$treat[plot_data_prc$treat ==2] = "silv"
# plot_data_prc$treat[plot_data_prc$treat ==3] = "silv+"
# plot_data_prc$plot.size = 6.25
# plot_data = rbind(plot_data, plot_data_prc[,colnames(plot_data),with=FALSE])


### DBH correction ####

#### Check for missing data ####
missing = data[, .(add_missing(yr = year, id = idtree)), .(plot)]
missing$idtree = tstrsplit(missing$V1, split = "\\s+")[[1]]
missing$year = as.numeric(tstrsplit(missing$V1, split = "\\s+")[[2]])
missing$V1 = NULL
data = merge(data,
             missing,
             by = c("plot", "idtree", "year"),
             all = TRUE)
status = data[, .(year, status=tree_status(dbh)), .(idtree)]
data = merge(data, status, by = c("idtree", "year"))
# remove status = 0 (dead trees) + not recruited yet// spot overgrown recruits??
data = subset(data, !is.na(status) & status == 1)
data = data[order(idtree, year)]


#### (c)  Correct dbh ####
# first_census = data[,min(year),.(site,plot)]
# colnames(first_census)[3] = "first_census"
# data = merge(data, first_census,by=c("site","plot"))
## check for duplicates ##
data$key = paste(data$idtree, data$year)
setkey(data, key)
data = data[!duplicated(data)]
data$key = NULL
data = data[order(idtree, year)]

## correction ##
data[dbh > 9.5 & dbh < 10, dbh := 10]
data_corr = data[, .(year = year,
                     dbh_c = dbh_correction(
                       X = dbh,
                       tm = year,
                       limit = minDBH,
                     )), by = .(idtree)]
# colnames(data_corr) = c("idtree", "year", "dbh_c")
if (any(is.na(data_corr$dbh_c))) { 
  print(paste(sum(is.na(data_corr$dbh_c)), "corrected DBH have NA value."))
  data_corr = subset(data_corr,!is.na(data_corr$dbh_c))}
data = merge(data, data_corr, by = c("idtree", "year"))
data$dbh_c[data$dbh_c < data$minDBH & !is.na(data$dbh_c)] <- NA

## plot problematic trees and their correction
mod_id <-
  unique(subset(data, dbh != dbh_c &
                  !is.na(dbh_c) &
                  !is.na(dbh))$idtree)
length(mod_id) # nb trees modified
par(mfrow = c(2,4), mar=c(1,1,3,1))
sapply(mod_id[!is.na(mod_id)][1:24], plot_corr)

### get agb, awp and awm
data[, genus := tstrsplit(name, " ")[[1]]]
data[, species := tstrsplit(name, " ")[[2]]]
data[, wd := getWoodDensity(genus, species)$meanWD]
data[, agb := computeAGB(dbh_c, wd, coord = c(-52.883, 5.3))]
setorder(data, plot, idtree, year)
data[, dagb := c(diff(agb)/diff(year), NA), .(idtree)]

cdata = stand_dynamics(id = data$idtree, 
                           var = data$agb,
                           dvar = data$dagb,
                           year = data$year,
                           plot = data$plot)

## standardize by ha ##
cdata[, `:=`(plot = as.numeric(plot), 
             agb = stock/6.25, 
             awp = gain/6.25, 
             awm = loss/6.25)]
setorder(cdata, plot, year)
cdata[, dT := c(diff(year), NA)]

## kohyama correction
cdata[, BS0 := agb - awm*dT]
cdata[, BT := agb + (awp-awm)*dT]
cdata[, awp := (log(BT/BS0)*(BT-agb))/(dT*log(BT/agb))]
cdata[, awm := (log(agb/BS0)*(BT-agb))/(dT*log(BT/agb))]


cdata = cdata[!is.na(awp), c("plot", "year", "agb", "awp", "awm")]

## add treatment info
cdata = merge(cdata, plot_data_prc, by = "plot")
ggplot(cdata, aes(x=year, y=awp, group = plot, color=as.factor(treat)))+geom_point()+geom_line()

###### SAVE DATA #######
save(cdata, file = "data/cdata.rda")

