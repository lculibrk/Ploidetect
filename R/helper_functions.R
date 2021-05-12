#' @export
load_ploidetect_data = function(path){
	if(!grepl("/$", path)){
		path = paste0(path, "/")
	}
	cna_raw = fread(paste0(path, "cna.txt"))
	cna_condensed = fread(paste0(path, "cna_condensed.txt"))
	original_data = readRDS(paste0(path, "segmented.RDS"))$all_data
}

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

#' @export
view_segments = function(cna, chrom, plot_pos, plot_end){
	plot_ploidetect(cna[chr == chrom][pos >= plot_pos & end <= plot_end])
}

#' @export
map_noisy_segments = function(cna, chrom, start_pos, end_pos, exclude){
	region = cna[chr == chrom][pos >= start_pos & end <= end_pos]
	top_cn = region[,.(s = sum(end - pos)), by = CN][order(s, decreasing = T)]$CN[1]
	top_segment = region[CN == top_cn][,.(s = sum(end - pos)), by = segment][order(s, decreasing = T)]$segment[1]
	mapping = unique(region$segment)
	regions = list(mapping)
	names(regions) = top_segment
	}

#' @export
most_common = function(v){
	res = names(sort(table(v),decreasing=TRUE))[1]
	return(res)

}


#' @export
curate_segments = function(cna, cna_condensed, regions){
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

