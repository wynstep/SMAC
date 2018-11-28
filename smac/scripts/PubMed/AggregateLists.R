#title           :AggregateLists.R
#description     :Aggregate ranked lists of MeSH terms
#author          :s.pirro
#date            :20180123
#version         :0.1
#notes           :
#==============================================================================

# loading packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(RankAggreg))
suppressMessages(library(data.table))
for (f in list.files(path="RankAggreg/", pattern="*.R")) {
  source(paste0("RankAggreg/",f))
}

# setup arguments
option_list <- list(
  make_option(c("-l", "--lists"),
              help="lists of MeSH terms to combine"),
  make_option(c("-o", "--output"),
              help="output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

# defining output folders
output = opt$output
data_dir = paste0(output,"/data/")
graph_dir = paste0(output,"/graphs/entropy_minimisation")

# loading MeSH lists
meshs = unlist(strsplit(opt$lists, ","))
mesh_level = fread(meshs[1],sep="\t",header=FALSE, stringsAsFactors = FALSE, data.table=FALSE, fill=TRUE)
mesh_laboundance = fread(meshs[2],sep="\t",header=FALSE, stringsAsFactors = FALSE, data.table=FALSE, fill=TRUE)
mesh_gaboundance = fread(meshs[3],sep="\t",header=FALSE, stringsAsFactors = FALSE, data.table=FALSE, fill=TRUE)

# find common meSH terms
common_mesh = Reduce(intersect, list(mesh_level[,1],mesh_laboundance[,1],mesh_gaboundance[,1]))

# subsetting MeSH lists according to common meSH terms
mesh_level.sub = mesh_level[which(mesh_level[,1] %in% common_mesh),]
mesh_laboundance.sub = mesh_laboundance[which(mesh_laboundance[,1] %in% common_mesh),]
mesh_gaboundance.sub = mesh_gaboundance[which(mesh_gaboundance[,1] %in% common_mesh),]
# Naming columns subsets
colnames(mesh_level.sub) = colnames(mesh_laboundance.sub) = colnames(mesh_gaboundance.sub) = c("mh","aboundance")

# compute ratio between local and global terms
tmp_bind = merge(mesh_laboundance.sub,mesh_gaboundance.sub, by = "mh")
mesh_ratio = data.frame(mh = tmp_bind[,1], ratio = tmp_bind[,2]/tmp_bind[,3], stringsAsFactors = F)

# sort MeSH lists
mesh_level.sorted = mesh_level.sub[order(-mesh_level.sub[,2]),] # descending
mesh_laboundance.sorted = mesh_laboundance.sub[order(-mesh_laboundance.sub[,2]),] # descending
mesh_gaboundance.sorted = mesh_gaboundance.sub[order(mesh_gaboundance.sub[,2]),] # ascending
mesh_ratio = mesh_ratio[order(-mesh_ratio$ratio),] # descending

# aggregate lists
all_lists = t(cbind(mesh_level.sorted[,1], mesh_laboundance.sorted[,1], mesh_gaboundance.sorted[,1], mesh_ratio[,1]))
aggregated_list = RankAggreg(all_lists, length(common_mesh), method="CE", distance="Spearman",
                  N=100, verbose=TRUE, output_dir=graph_dir)

# save aggregated list
write(aggregated_list$top.list, file = paste0(data_dir,"/aggregated_mesh.txt"), append = FALSE, sep = "\n")

# save informations
sink(paste0(data_dir,"/merging_log.txt"))
cat("### Merging log -- SMAC 2.0 ### \n
      Method for merging: ",aggregated_list$method,"\n
      Distance: ",aggregated_list$distance,"\n
      Number of iterations until convergence: ",aggregated_list$num.iter,"\n
      Minimum reached value: ",aggregated_list$optimal.value,"")
sink()