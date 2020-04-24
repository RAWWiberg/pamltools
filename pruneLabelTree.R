#!/usr/bin/Rscript
# Prune Macrostomum species tree to include only tips for orthogroups.
library("ape")
library("phytools")
library("optparse")
#---------------------------------------------------#
# Testing Variables
#set<-"ovaries2"
#sp<-c("Maclig_20180705","Maclig_37v3","Mac112_3462_20180705","Macmin_3429_20180705","Macmed_2294_20180705","Mac003_3326_20180705","Mac088_2936_20180705","Mac101_3113_20180705","Mac108_3220_20180705","Mac094_2994_20180705","Mac103_3118_20180705")
#---------------------------------------------------#
treefile<-"~/Projects/macrostomum/data/transcriptome_assemblies/2018_assemblies/phylogenomics/SM_Dez18_astral_subset.tre"
# Options List
option_list = list(
  make_option(c("-s", "--set"), type="character", default=NULL, 
              help="which gene set", metavar="character"),
  make_option(c("-l", "--labels"), type="character", default=NULL, 
              help="file of table of branches to label", metavar="character")
  );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# Program
cat("Pruning trees for: ",opt$set,"...\n")
# List the files for each orthogroup
#files<-list.files(paste("~/Projects/macrostomum/data/arbore_2014/orthofinder_ogs/",set,"/cds",sep=""),pattern = ".list")
files<-list.files(paste("~/Projects/macrostomum/data/reproduction_related_genes/",
                        opt$set,"/alns_noRedBris",sep=""),
                  pattern = "_sp.list")
dat<-data.frame(og=gsub("_sp.list","",files),
                nsp=vector(length=length(files)),
                labs=vector(length=length(files)),
                n_labsp=vector(length=length(files)),
                medlen=vector(length=length(files)),
                maxlen=vector(length=length(files)),
                labsp=vector(length=length(files)))

# Read in the full tree
mac_tree<-read.tree(treefile)

if(is.null(opt$labels)){
  for(f in files){
    # Get base name for OrthoGroup
    f<-gsub(pattern = ".list",replacement = "",f)
    # Scan the species list
    sp<-scan(paste("~/Projects/macrostomum/data/reproduction_related_genes/",opt$set,"/alns_noRedBris/",f,".list",sep=""),what = "character")
    dat$nsp[dat$og==gsub("_sp","",f)]<-length(sp)
    # If nr species < 4; don't do anything
    if(length(sp) >= 4){
      # Find the species to drop
      toDrop<-mac_tree$tip.label[which(!(mac_tree$tip.label %in% sp))]
      # Drop tips
      mac_tree_sub<-drop.tip(mac_tree,tip = toDrop)
      # Make a figure of the tree
      #png(filename = paste("~/Projects/macrostomum/data/arbore_2014/orthofinder_ogs/",opt$set,"/cds/",f,"_tree.png",sep=""))
      #ggtree(mac_tree) + geom_tiplab(size = 4,linesize = 0.5,align = TRUE)+ggplot2::xlim(0,0.2)
      #dev.off()
      # Calculate the median and maximum pairwise distance for each species to Mlig_37v3
      not_mlig<-mac_tree_sub$tip.label[grep("Maclig_37v3",mac_tree_sub$tip.label,invert=TRUE)]
      #medlen<-median(unlist(lapply(not_mlig,function(x){return(fastDist(tree=mac_tree,sp1="Maclig_37v3",sp2=x))})))
      #dat$medlen[dat$og==gsub("_sp","",f)]<-medlen
      #maxlen<-max(unlist(lapply(not_mlig,function(x){return(fastDist(tree=mac_tree,sp1="Maclig_37v3",sp2=x))})))
      #dat$maxlen[dat$og==gsub("_sp","",f)]<-maxlen
      
      # Write an unrooted nexus tree
      mac_tree_sub$node.label<-NULL
      mac_tree_sub$edge.length<-NULL

      write.tree(unroot(mac_tree_sub),
                 file = paste("~/Projects/macrostomum/data/reproduction_related_genes/",
                              opt$set,"/alns/",f,".tree",sep=""))
    }
  }
  write.table(dat,
              paste("~/Projects/macrostomum/data/reproduction_related_genes/",opt$set,"/alns_noRedBris/og_species_data_null.tab",sep=""),
              row.names = FALSE,col.names=TRUE,quote=FALSE,sep="\t")
}else{
  # Read in the labelling table
  lab_spp<-read.table(paste("~/Projects/macrostomum/data/morphology_behaviour/",opt$labels,".csv",sep=""),
                      header=TRUE,sep=",")
  # If tip labels include species in lab_spp, add " #X" to the name.
  # number of labels depends on number of columns in the lab_spp table
  nspp<-vector(length=ncol(lab_spp))
  labs<-paste("#",0:(ncol(lab_spp)-1),sep="")
  ids<-paste(labs,"=",colnames(lab_spp),sep="")
  write.table(cbind(labs,ids),
              paste("~/Projects/macrostomum/data/reproduction_related_genes/",
                    opt$set,
                    "/alns_noRedBris/",
                    opt$labels,"_treelabels.tab",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep=",")
  for(f in files){
    # Get base name for OrthoGroup
    f<-gsub(pattern = ".list",replacement = "",f)
    # Scan the species list
    sp<-scan(paste("~/Projects/macrostomum/data/reproduction_related_genes/",opt$set,"/alns_noRedBris/",f,".list",sep=""),
             what = "character")
    dat$nsp[dat$og==gsub("_sp","",f)]<-length(sp)
    if(length(sp)>=4){
      # Find the species to drop
      toDrop<-mac_tree$tip.label[which(!(mac_tree$tip.label %in% sp))]
      # Drop tips
      mac_tree_sub<-drop.tip(mac_tree,tip = toDrop)
      # Make a figure of the tree
      #png(filename = paste("~/Projects/macrostomum/data/arbore_2014/orthofinder_ogs/",opt$set,"/cds/",f,"_tree.png",sep=""))
      #ggtree(mac_tree) + geom_tiplab(size = 4,linesize = 0.5,align = TRUE)+ggplot2::xlim(0,0.2)
      #dev.off()
      # Calculate the median and maximum pairwise distance for each species to Mlig_37v3
      not_mlig<-mac_tree_sub$tip.label[grep("Maclig_37v3",mac_tree_sub$tip.label,invert=TRUE)]
      #medlen<-median(unlist(lapply(not_mlig,function(x){return(fastDist(tree=mac_tree,sp1="Maclig_37v3",sp2=x))})))
      #dat$medlen[dat$og==gsub("_sp","",f)]<-medlen
      #maxlen<-max(unlist(lapply(not_mlig,function(x){return(fastDist(tree=mac_tree,sp1="Maclig_37v3",sp2=x))}))) 
      #dat$maxlen[dat$og==gsub("_sp","",f)]<-maxlen
      
      # Write an unrooted nexus tree (null tree)
      mac_tree_sub$node.label<-NULL
      mac_tree_sub$edge.length<-NULL
      
      write.tree(unroot(mac_tree_sub),
                 file = paste("~/Projects/macrostomum/data/reproduction_related_genes/",
                              opt$set,"/alns_noRedBris/",f,".tree",sep=""))
      n_labsp<-rep(0,ncol(lab_spp))
      labsp<-vector(length = ncol(lab_spp))
      for(col in 1:ncol(lab_spp)){
        lab<-labs[col]
        n_labsp[col]<-length(which(mac_tree_sub$tip.label %in% lab_spp[,col]))
        #dat$n_injsp[dat$og==gsub("_sp","",f)]<-n_inj
        labsp[col]<-paste(mac_tree_sub$tip.label[
          which(mac_tree_sub$tip.label %in% lab_spp[,col])],collapse = "|")
        spp<-mac_tree_sub$tip.label[which(mac_tree_sub$tip.label %in% lab_spp[,col])]
        spp<-paste(spp," ",lab,sep="")
        mac_tree_sub$tip.label[which(mac_tree_sub$tip.label %in% lab_spp[,col])]<-spp
      }
      #cat(n_labsp,"\n")
      #cat(all(n_labsp >= 1),"\n")
      dat$n_labsp[dat$og==gsub("_sp","",f)]<-paste(n_labsp,collapse=":")
      dat$labsp[dat$og==gsub("_sp","",f)]<-paste(labsp,collapse=":")
      if(all(n_labsp >= 1)){
        # Write an unrooted nexus tree (labelled tree)
        mac_tree_sub$node.label<-NULL
        mac_tree_sub$edge.length<-NULL
        write.tree(unroot(mac_tree_sub),
                   file = paste("~/Projects/macrostomum/data/reproduction_related_genes/",
                                opt$set,"/alns_noRedBris/",f,"_",opt$labels,".tree",sep=""))
      }
    }
  }
  write.table(dat,
              paste("~/Projects/macrostomum/data/reproduction_related_genes/",opt$set,
                    "/alns_noRedBris/og_species_data_",opt$labels,".tab",sep=""),
              row.names = FALSE,col.names=TRUE,quote=FALSE,sep="\t")
}


