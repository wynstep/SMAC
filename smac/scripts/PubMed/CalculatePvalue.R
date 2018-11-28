#title           :CalculatePvalue.R
#description     :Calculate p-value fluctuations
#author          :s.pirro
#date            :20180123
#version         :0.1
#notes           :
#==============================================================================

# loading packages
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))

# setup arguments
option_list <- list(
  make_option(c("-l", "--list"),
              help="lists of combined meSH term"),
  make_option(c("-o", "--output"),
              help="output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

list_file = opt$list
outdir = opt$output
mesh_dict_file = "src/mesh_dict.txt"

# loading list of merged terms
merged_list <- readLines(list_file)

# loading list of all mesh terms
mesh_dict <- readLines(mesh_dict_file)

# creating ordered list of positions (first vector)
query <- c(1:length(merged_list))
# iterating 1000 times -- parameter can be changed later
# create tau and p-value dataframes
tau_container = data.frame(iteration=numeric(), tau=numeric())
pvalue_container = data.frame(iteration=numeric(), pvalue=numeric())

for (i in 1:1000) {
  Sys.sleep(0.1)
  cat("\r", paste0("Calculating tau and p-value for iteration ", i))
  # create vector of random positions (from mesh dictionary, as background)
  background <- sample(x = c(1:length(mesh_dict)), size = length(query), replace = FALSE)
  # calculate the kendall distance and the p-value
  distance <- cor.test(query,background,method = "kendall")
  # report and store tau distance and p-values (for each iteration)
  tau_container[i,] <- c(i, distance$estimate)
  pvalue_container[i,] <- c(i, distance$p.value)
}

# create plots for p-value and tau fluctuations
png(filename = paste0(outdir,"/graphs/tau_fluctuations.png"))
  mean_tau = round(mean(tau_container$tau), digits = 3)
  tau_plot <- ggplot(data=tau_container, aes(x=iteration, y=tau, group=1, colour = tau)) + geom_line(alpha = 0.1) + geom_point(alpha = 0.2)
  tau_plot + scale_color_gradient(low="blue", high="red") +
  geom_hline(yintercept = mean_tau) +
  annotate("text", min(tau_container$iteration+200), mean_tau, vjust = -1, label = mean_tau)
dev.off()

png(filename = paste0(outdir,"/graphs/pvalue_fluctuations.png"))
  mean_pvalue = round(mean(pvalue_container$pvalue), digits = 3)
  pvalue_plot <- ggplot(data=pvalue_container, aes(x=iteration, y=pvalue, group=1, colour = pvalue)) + geom_line(alpha = 0.1) + geom_point(alpha = 0.2)
  pvalue_plot + scale_color_gradient(low="blue", high="red") +
  geom_hline(yintercept = mean_pvalue) +
  annotate("text", min(pvalue_container$iteration+200), mean_pvalue, vjust = -1, label = mean_pvalue)
dev.off()

# save p_values and fluctuations details
write.table(pvalue_container, file = paste0(outdir,"/data/p_value_fluctuations.csv"), sep = ",", quote = FALSE, col.names = T, row.names = F)
write.table(tau_container, file = paste0(outdir,"/data/tau_fluctuations.csv"), sep = ",", quote = FALSE, col.names = T, row.names = F)
