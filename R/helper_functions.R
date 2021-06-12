

#' Load data from Ploidetect's output
#' 
#' \code{load_ploidetect_data} takes the path to the ploidetect output and returns a named \code{list} containing \code{cna},
#' a data.table of the bin-level copy number data, \code{segments}, a data.table of the 
#' segment-level copy number data, and \code{raw}, a data.table of the raw data.
#' @param path the path to the ploidetect outputs
#' @return a named list oif data.tables
#' @export
load_ploidetect_data = function(path){
	if(!grepl("/$", path)){
		path = paste0(path, "/")
	}
	cna_raw = fread(paste0(path, "cna.txt"))
	cna_condensed = fread(paste0(path, "cna_condensed.txt"))
	original_data = readRDS(paste0(path, "segmented.RDS"))$all_data
	out = list("cna" = cna_raw, "segments" = cna_condensed, "raw" = original_data)
	return(out)
	}



#' Produce ploidetect-style plots for a copy number profile
#' 
#' \code{plot_ploidetect} takes a data.table containing the bin-level copy number data for a
#' locus of interest and an optional argument of whether to indicate segment breakpoints. It returns a 
#' \code{ggplot2} plot of the data
#' @param cna a data.table containing the regions of interest. Single chromosome only
#' @param seg_lines a boolean of whether to include vertical dashed lines at each segment breakpoint
#' @return a ggplot2 plot
#' @export
plot_ploidetect = function(cna, seg_lines = F){
	CN_palette <- c("0" = "#000000", 
									"1" = "#000066", 
									"2" = "#26d953", 
									"3" = "#609f70", 
									"4" = "#cccc00", 
									"5" = "#ADAD47",
									"6" = "#cc6600", 
									"7" = "#856647", 
									"8" = "#cc0000"
	)
	plot_labs = c("0" = "HOMD", 
								"1" = "CN = 1", 
								"2" = "CN = 2 HET", 
								"3" = "CN = 2 HOM", 
								"4" = "CN = 3 HET", 
								"5" = "CN = 3 HOM", 
								"6" = "CN = 4 HET", 
								"7" = "CN = 4 HOM", 
								"8" = "CN = 5+")
	plot_shapes <- c("0" = 4, 
									 "1" = 19, 
									 "2" = 19, 
									 "3" = 19, 
									 "4" = 19, 
									 "5" = 19,
									 "6" = 19, 
									 "7" = 19, 
									 "8" = 19)
	
	CNA_plot = ggplot(cna, aes(x = pos, y = log(corrected_depth, base = 2), color = as.character(state))) + 
			geom_point_rast(size = 0.5, aes(shape = factor(state))) +
			scale_shape_manual(values = plot_shapes, 
												 labels = plot_labs,
												 name = "State") +
			scale_color_manual(name = "State",
												 values = CN_palette, 
												 labels = plot_labs) + 
			ylab("log(Read Depth)") + 
			xlab("Position") + 
			ggtitle(paste0("Copy number profile")) + 
			theme_bw()
	
	if(seg_lines){
		breakpoints = cna[segment != data.table::shift(segment)]$pos
		CNA_plot = CNA_plot + geom_vline(xintercept = breakpoints, alpha = 0.3, linetype = 2)
	}
	return(CNA_plot)
}

#' Produce ploidetect-style plots for a copy number profile at the given chromosome and positions
#' 
#' \code{view_segments} is a wrapper for plot_ploidetect to make filtering for specific regions simple
#' @param cna a data.table containing the bin-level copy number calls
#' @param chrom the chromosome of interest
#' @param plot_pos the start position of the plot
#' @param plot_end the end position of the plot
#' @return a ggplot2 plot
#' @export
view_segments = function(cna, chrom, plot_pos, plot_end){
	plot_ploidetect(cna[chr == chrom][pos >= plot_pos & end <= plot_end])
}

#' Map all segments within an interval to the dominant CN state 
#' 
#' \code{map_noisy_segments} takes in a genomic region on chromosome \code{chr} between \code{start_pos}
#' and \code{end_pos} and produces a named list of vectors to be passed to \code{curate_segments} 
#' indicating a mapping between each segment and the most common segment in the region
#' 
#' @param cna a data.table containing the bin-level copy number calls
#' @param chrom the chromosome of interest
#' @param start_pos the start position of the mapping
#' @param end_pos the end position of the mapping - the segments of interest must have end <= this number
#' @return a named list
#' @export
map_noisy_segments = function(cna, chrom, start_pos, end_pos){
	region = cna[chr == chrom][pos >= start_pos & end <= end_pos]
	top_cn = region[,.(s = sum(end - pos)), by = CN][order(s, decreasing = T)]$CN[1]
	top_segment = region[CN == top_cn][,.(s = sum(end - pos)), by = segment][order(s, decreasing = T)]$segment[1]
	mapping = unique(region$segment)
	regions = list(mapping)
	names(regions) = top_segment
	return(regions)
	}

#' @export
most_common = function(v){
	res = names(sort(table(v),decreasing=TRUE))[1]
	return(res)

}


#' Merge noisy segments with their locally dominant segment 
#' 
#' \code{curate_segments} takes in \code{cna} of the bin-level copy number calls and
#' \code{regions} from \code{map_noisy_segments}. Multiple outputs from \code{map_noisy_segments}
#' can be combined with \code{c()}
#' 
#' @param cna a data.table containing the bin-level copy number calls
#' @param regions the output from map_noisy_segments
#' @return a \code{cna} data.table with the modified segments.
#' @export
curate_segments = function(cna, regions){
	regions_list =  unlist(regions)
	if(any(duplicated(regions_list))){
		repeated = regions_list[duplicated(regions_list)]
		stop(paste0("Segments in regions must not be re-used! \nOffending segments: ", paste0(repeated, collapse = ", ")))
	}
	for(seg in names(regions)){
		from = regions[[seg]]
		to = as.integer(seg)
		extract = cna[segment %in% from]
		others = cna[!segment %in% from]
		extract[,segment:=to]
		extract[,segment_depth := median(corrected_depth), by = segment]
		extract[,CN := as.numeric(most_common(CN)), by = segment]
		extract[,state := as.numeric(most_common(state)), by = segment]
		extract[,zygosity:=most_common(zygosity), by = segment]
		cna = rbind(others, extract)[order(chr, pos)]
	}
	return(cna)
}

#' @export
collapse_segments = function(cna){
	cna[,.(pos = first(pos), end = last(end), CN = first(CN), state = first(state), zygosity = first(zygosity), segment_depth = first(segment_depth)), by = list(chr, segment)]
}

#' @export
plot_all_segments = function(cna){
	data(centromeres)
	CN_palette <- c("0" = "#000000", 
									"1" = "#000066", 
									"2" = "#26d953", 
									"3" = "#609f70", 
									"4" = "#cccc00", 
									"5" = "#ADAD47",
									"6" = "#cc6600", 
									"7" = "#856647", 
									"8" = "#cc0000"
	)
	plot_labs = c("0" = "HOMD", 
								"1" = "CN = 1", 
								"2" = "CN = 2 HET", 
								"3" = "CN = 2 HOM", 
								"4" = "CN = 3 HET", 
								"5" = "CN = 3 HOM", 
								"6" = "CN = 4 HET", 
								"7" = "CN = 4 HOM", 
								"8" = "CN = 5+")
	plot_shapes <- c("0" = 4, 
									 "1" = 19, 
									 "2" = 19, 
									 "3" = 19, 
									 "4" = 19, 
									 "5" = 19,
									 "6" = 19, 
									 "7" = 19, 
									 "8" = 19)
	
	cna = split(cna, f = cna$chr)
	CNA_plot <- lapply(cna, function(x){

		chr = x$chr[1]
		x %>% filter(end < centromeres$pos[which(centromeres$chr %in% chr)[1]] | pos > centromeres$end[which(centromeres$chr %in% chr)[2]]) %>% ggplot(aes(x = pos, y = log(corrected_depth, base = 2), color = as.character(state))) + 
			geom_point_rast(size = 0.5, aes(shape = factor(state))) +
			scale_shape_manual(values = plot_shapes, 
												 labels = plot_labs,
												 name = "State") +
			scale_color_manual(name = "State",
												 values = CN_palette, 
												 labels = plot_labs) + 
			ylab("log(Read Depth)") + 
			xlab("position") + 
			ggtitle(paste0("Chromosome ", chr, " copy number profile")) + 
			theme_bw()
	})
	vaf_plot <- lapply(cna, function(x){
		chr = x$chr[1]
		x %>% filter(end < centromeres$pos[which(centromeres$chr %in% chr)][1] | pos > centromeres$end[which(centromeres$chr %in% chr)][2]) %>% filter(!is.na(maf)) %>% ggplot(aes(x = pos, y = unlist(unmerge_mafs_grouped(maf, flip = T)), color = as.character(state))) + 
			geom_point_rast(size = 0.5, aes(shape = factor(state))) + 
			scale_shape_manual(values = plot_shapes, 
												 labels = plot_labs,
												 name = "State") +
			scale_color_manual(name = "State",
												 values = CN_palette, 
												 labels = plot_labs) + 
			ylab("Major allele frequency") + 
			xlab("position") + 
			ggtitle(paste0("Chromosome ", chr, " allele frequency profile")) + 
			scale_y_continuous(limits = c(0.5, 1)) +
			theme_bw()
	})
	
	
	
	cna_plots <- list()
	
	for(i in 1:length(CNA_plot)){
		cna_plots[i] <- list(plot_grid(CNA_plot[[i]], vaf_plot[[i]], align = "v", axis = "l", ncol = 1))
	}
	chrs = suppressWarnings(as.numeric(names(cna)))
	sortedchrs = sort(chrs)
	chrs = c(sortedchrs, names(CN_calls)[is.na(chrs)])
	
	cna_plots = cna_plots[order(order(chrs))]

	return(cna_plots)
	
}
#' @export
adjust_model_ploidy = function(models, model_ind = 1, new_ploidy){
	model = models[model_ind,]
	tp = model$tp
	rtp = 1 - model$tp
	diff = model$modeling_tp
	position_2 = 2*diff/(1-rtp)
	old_ploidy = model$ploidy
	## Calculate maxpeak position
	mp = position_2 + diff * (old_ploidy - new_ploidy)
	new_hd = mp - new_ploidy*diff
	new_pos2 = mp - diff * (new_ploidy - 2)
	new_tp = 1 - new_hd/new_pos2
	model2 = model
	model2$tp = new_tp
	model2$ploidy = new_ploidy
	return(model2)
}
