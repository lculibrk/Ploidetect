#!/usr/bin/Rscript

' plot_regions.R

Usage:
plot_regions.R -d ploidetect_dir -o outfile -t type -b cyto_path (--pos position | --chrom chromosome --start start --end end)

Options:
-d --dir ploidetect_dir     directory of ploidetect output
-o --out outfile            output file
-t --type type              color scheme (loh vs cnv)
-b --band cyto_path         path to cytobands file            
-p --pos position           positions formatted as chr:start-end. Either this or all of chrom, start, and end must be provided
    --chrom chromosome      chromosome
    --start start           start
    --end end               end

' -> doc


#
# Load libraries
library(docopt)
library(devtools)
#library(dplyr)
#library(GenomicRanges)
args = docopt(doc)
library(Ploidetect)

data = load_ploidetect_data(args$dir)

print(str(data))

cnv_data = data$cna

cn_positions = data$cn_positions

if(!is.null(args$pos)){
    posstring = args$pos
    chr = gsub("([^\\:]*)\\:([^-]*)-([^-]*)", "\\1", posstring)
    start = as.numeric(gsub("([^\\:]*)\\:([^-]*)-([^-]*)", "\\2", posstring))
    end = as.numeric(gsub("([^\\:]*)\\:([^-]*)-([^-]*)", "\\3", posstring))
}else{
    chr = args$chr
    start = args$start
    end = args$end
}

type = args$type

cytobands = fread(args$band)
names(cytobands) =  c("chr", "pos", "end", "band", "type")

p = focus_view(cnv_data, chr, start, end, type, cytobands = cytobands)

png(args$out, type = "cairo")
p
dev.off()
