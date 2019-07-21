library(geiger)
library(phytools)

#load in trees and data (tab delimited file)
setwd("/Users/joycechery/Desktop/Science/PNAS/New_Analysis/")
tree <- read.nexus("infile.nex.run1_run6Mfrontandremovedrunswithsplitfreq02andmoreMburninRemoved_mcc.tre")
tree2<- root(tree, "Matayba_guianensis", resolve.root=TRUE)
tree2$edge.length <- tree2$edge.length * 6596
chronogram <- chronos(tree2, lambda=0.5, model="relaxed", control= chronos.control())

#pruned tree
species <- c("Paullinia_acutangula_a"
             , "Paullinia_acutangula_b"
             , "Paullinia_alata_a"
             , "Paullinia_alata_b"
             , "Paullinia_alata_c"
             , "Paullinia_alata_d"
             , "Paullinia_allenii_a"
             , "Paullinia_allenii_b"
             , "Paullinia_alsmithii_a"
             , "Paullinia_alsmithii_b"
             , "Paullinia_baileyi_a"
             , "Paullinia_baileyi_b"
             , "Paullinia_bilobulata"
             , "Paullinia_boliviana_a"
             , "Paullinia_bracteosa_b"
             , "Paullinia_bracteosa_c"
             , "Paullinia_bracteosa_d"
             , "Paullinia_caloptera_a"
             , "Paullinia_carpopodea_a"
             , "Paullinia_carpopodea_b"
             , "Paullinia_carpopodea_c"
             , "Paullinia_cf.boliviana_b"
             , "Paullinia_cf.nobilis"
             , "Paullinia_cf.stellata_a"
             , "Paullinia_clathrata_a"
             , "Paullinia_costaricensis"
             , "Paullinia_dasystachya"
             , "Paullinia_elegans"
             , "Paullinia_elongata_a"
             , "Paullinia_elongata_b"
             , "Paullinia_eriocarpa"
             , "Paullinia_exalata"
             , "Paullinia_faginea"
             , "Paullinia_fibrigera_a"
             , "Paullinia_fibrigera_b"
             , "Paullinia_fimbriata"
             , "Paullinia_fuscescens_a"
             , "Paullinia_fuscescens_b"
             , "Paullinia_glomerulosa_a"
             , "Paullinia_glomerulosa_b"
             , "Paullinia_hystrix_a"
             , "Paullinia_hystrix_b"
             , "Paullinia_hystrix_c"
             , "Paullinia_hystrix_d"
             , "Paullinia_imberbis"
             , "Paullinia_ingifolia_a"
             , "Paullinia_ingifolia_b"
             , "Paullinia_ingifolia_c"
             , "Paullinia_ingifolia_d"
             , "Paullinia_jamaicensis"
             , "Paullinia_josecuatrii"
             , "Paullinia_killipii_a"
             , "Paullinia_killipii_b"
             , "Paullinia_largifolia"
             , "Paullinia_latifolia"
             , "Paullinia_leiocarpa"
             , "Paullinia_microneura"
             , "Paullinia_neglecta_a"
             , "Paullinia_neglecta_b"
             , "Paullinia_obovata_a"
             , "Paullinia_obovata_b"
             , "Paullinia_olivacea"
             , "Paullinia_paullinoides_a"
             , "Paullinia_paullinoides_b"
             , "Paullinia_pinnata_a"
             , "Paullinia_pinnata_b"
             , "Paullinia_pinnata_c"
             , "Paullinia_pinnata_d"
             , "Paullinia_pseudota_a"
             , "Paullinia_pseudota_b"
             , "Paullinia_pseudota_c"
             , "Paullinia_rubiginosa_a"
             , "Paullinia_rubiginosa_b"
             , "Paullinia_rubiginosa.subsp.setosa"
             , "Paullinia_rufescens"
             , "Paullinia_rugosa_a"
             , "Paullinia_rugosa_b"
             , "Paullinia_rugosa_c"
             , "Paullinia_selenoptera"
             , "Paullinia_serjaniifolia"
             , "Paullinia_simulans"
             , "Paullinia_sp._a"
             , "Paullinia_sp._b"
             , "Paullinia_sphaerocarpa"
             , "Paullinia_spicata"
             , "Paullinia_sprucei"
             , "Paullinia_stellata_b"
             , "Paullinia_stipitata_a"
             , "Paullinia_stipitata_b"
             , "Paullinia_tomentosa"
             , "Paullinia_trigonia"
             , "Paullinia_turbacensis_a"
             , "Paullinia_turbacensis_b"
             , "Paullinia_turbacensis_c"
             , "Paullinia_turbacensis_d"
             , "Paullinia_turbacensis_d")

pruned_chronogram<-drop.tip(chronogram,chronogram$tip.label[!(chronogram$tip.label %in% species)])
pruned_chronogram<-ladderize(pruned_chronogram, F )
class(pruned_chronogram)<-"phylo"

#Primary growth character evolution
primary<-read.delim("primary.txt" , sep = "\t", row.names = 1)
primary<-setNames(primary[,1],rownames(primary))

#fit the primary growth data model of evolution to the tree  -ASR
er<-make.simmap(pruned_chronogram, primary, model ="ER", pi="estimated")  #$1df
sym<-make.simmap(pruned_chronogram, primary, model ="SYM", pi="estimated")#15
ard<-make.simmap(pruned_chronogram, primary, model ="ARD", pi="estimated") #30
1-pchisq(2*abs(ard$logL- er$logL), 1)
#ER is prefer yall.

#simulate character history
colrs<-c("black", "olivedrab4")
cols<-setNames((colrs)[1:length(unique(primary))],sort(unique(primary)))
ER_100<-make.simmap(pruned_chronogram, primary, model ="ER", nsim=100, pi="estimated")
densityTree(ER_100,method="plotSimmap",lwd=7,nodes="intermediate", 
            colors=cols,ylim=c(3,92),compute.consensus=FALSE,
            fsize=.65, show.axis = F, direction="rightwards")
add.simmap.legend(colors=cols,prompt=FALSE,x=90, y=90 )
describe.simmap(ER_100)
countSimmap(ER_100)
describe.simmap(er)
