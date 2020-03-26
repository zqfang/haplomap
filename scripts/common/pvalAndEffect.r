#Susan's modified version of the the following r script, found in CURIS/DILL

# Try to pick out the best arrays using Susan's pvalues and
# effect size of genes (computed on dippy with my maxEffects perl
# script)
EFFECT_SIZES="effectSizes_reprodMoreSubs.txt"
SORTED_BY_PVAL="sorted_by_pvalue_reprodMoreSubs.txt"
SORT_BY_EFFECT_OUTPUT="sort_by_effect_reprodMoreSubs.txt"
SORT_BY_GENE_OUTPUT="sort_by_gene_reprodMoreSubs.txt"
GENE_AND_FREQ="gene_and_freq_reprodMoreSubs.txt"

options(stringsAsFactors = FALSE)

effectSize <- read.delim(EFFECT_SIZES, header=FALSE)
names(effectSize) <- c("varnameAndMeasnum", "sex", "gene", "aasubst", "pvalue", "effect")

sbpv <- read.csv(SORTED_BY_PVAL)

sbpv1 <- sbpv[sbpv$p_vals < 1.0e-5,]

sbpv1$varnameAndMeasnum <-paste(sbpv1$varname,"_",sbpv1$measnum,sep="")
sbpv1$varname <- NULL
sbpv1$measnum <- NULL


xx <- merge(sbpv1, effectSize, by=c("varnameAndMeasnum", "sex"))
xxx <- xx[xx$effect >= 0.85,]

o <- order(xxx$effect, decreasing=TRUE)
write.table(xxx[o, c('effect', 'varnameAndMeasnum', 'sex', 'gene', 'aasubst', 'p_vals') ], row.names=FALSE, quote=FALSE, sep="\t",file=SORT_BY_EFFECT_OUTPUT)

# sort by genes.
o1 <- order(xxx$gene)
write.table(xxx[o1, c('effect', 'varnameAndMeasnum', 'sex', 'gene', 'aasubst', 'p_vals') ], row.names=FALSE, quote=FALSE, sep="\t",file=SORT_BY_GENE_OUTPUT)

#make a list of all the genes and how many times they appear
genes <-data.frame(table(xxx$gene))
names(genes) <-c('gene','Freq')
o2 <-order(genes$Freq,decreasing=TRUE)
write.table(merge(genes[o2,],xxx,by=c('gene'),sort=FALSE),row.names=FALSE,quote=FALSE,sep="\t",file=GENE_AND_FREQ)
