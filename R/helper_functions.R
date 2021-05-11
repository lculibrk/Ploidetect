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

map_noisy_segments = function(cna, chrom, start_pos, end_pos, exclude){
	region = cna[chr == chrom][pos >= start_pos & end <= end_pos]
	top_cn = region[,.(s = sum(end - pos)), by = CN][order(s, decreasing = T)]$CN[1]
	top_segment = region[CN == top_cn][,.(s = sum(end - pos)), by = segment][order(s, decreasing = T)]$segment[1]
	mapping = unique(region$segment)
	regions = list(mapping)
	names(regions) = top_segment
	}

#' @export
curate_segments = function(cna, cna_condensed, regions){
	
}



