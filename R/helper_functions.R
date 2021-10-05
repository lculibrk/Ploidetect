

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
	
	if(file.exists(paste0(path, "segmented.RDS"))){
		original_data = readRDS(paste0(path, "segmented.RDS"))$all_data
	}else if(file.exists(gsub("output", "data", paste0(path, "segmented.RDS")))){
		original_data = readRDS(gsub("output", "data", paste0(path, "segmented.RDS")))
	}
	
	cn_positions = readRDS(paste0(path, "metadat.rds"))[[1]]
	
	out = list("cna" = cna_raw, "segments" = cna_condensed, "raw" = original_data, cn_positions = cn_positions)
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
plot_ploidetect = function(cnv_data, cn_positions, cytobands, mode = "all", seg_lines = F, colors = "cnv"){
	if(!mode %in% c("all", "zoomed")){
		stop("Invalid mode selected")
	}
	
	if(colors == "cnv"){
		CN_palette <- c("#000000FF", # Black
										"#2aa4caFF", # Blue
										"#979797FF", # Light grey
										"#757575FF", # grey
										"#f6b75eFF", # Light Orange
										"#f39420FF", # Orange
										"#89b99aFF", # Light Green
										"#529977FF", # Green
										"#ec0077FF" # pink
		)
		names(CN_palette) = c(0:8)
		plot_labs = c("0" = "HOMD", 
									"1" = "CN = 1", 
									"2" = "CN = 2 HET", 
									"3" = "CN = 2 HOM", 
									"4" = "CN = 3 HET", 
									"5" = "CN = 3 HOM", 
									"6" = "CN = 4 HET", 
									"7" = "CN = 4 HOM", 
									"8" = "CN = 5+")

	}else if(colors == "loh"){
		CN_palette <- c(
			"#f39420FF", # Orange/LOH
			"#000000FF"  # Black/HET
		)
		names(CN_palette) = c("HOM", "HET")
		plot_labs = c("HOM" = "HOM",
									"HET" = "HET")
	}else{
		stop("must set \"cnv\" or \"loh\" for colors")
	}
	
	if(!"chr" %in% names(cytobands)){
		names(cytobands) = c("chr", "pos", "end", "lab", "stain")
	}

	plot_shapes <- c("0" = 4, 
									 "1" = 19, 
									 "2" = 19, 
									 "3" = 19, 
									 "4" = 19, 
									 "5" = 19,
									 "6" = 19, 
									 "7" = 19, 
									 "8" = 19)

	plot_calls = cnv_data
	
	plt_positions = cn_positions[which(as.numeric(names(cn_positions)) == round(as.numeric(names(cn_positions))))]
	plt_positions = log2(plt_positions[names(plt_positions) %in% c(0:5)])
	
	if(min(2^plt_positions) - diff(2^plt_positions)[1] <= 0){
		minbound = min(log2(plot_calls[corrected_depth > 0]$corrected_depth))
	}else{
		minbound = log2(min(2^plt_positions) - diff(2^plt_positions)[1])
	}
	maxbound = log2(max(2^plt_positions) +  15 * diff(2^plt_positions)[1])
	
	plot_calls[,c("overflow", "underflow"):=F]
	plot_calls[suppressWarnings(log2(corrected_depth)) > maxbound, overflow:=T]
	plot_calls[suppressWarnings(log2(corrected_depth)) < minbound, underflow:=T]
	plot_calls[,corrected_depth:=pmin(2^maxbound, corrected_depth)]
	plot_calls[,corrected_depth:=pmax(2^minbound, corrected_depth)]
	
	## Lower bound on CN
	states = c(0:8, 8)
	states = data.table("state" = states, state_cn = c(0, 1, 2, 2, 3, 3, 4, 4, 5, 5), zygosity = c("HOM", "HOM", "HET", "HOM", "HET", "HOM", "HET", "HOM", "HET", "HOM"))
	legend_plot_fn = function(states, colors = colors){
		if(colors == "loh"){
			plt_df = data.table("state" = c("HOM", "HOM", "HET"), "cols" = CN_palette[c(1, 1, 2)], "plot_labs" = plot_labs[c(1,1,2)], "plot_shapes" = rep(19, 3))
			plt_df$y = c(NA, 1, 1)
			plt_df$x = c(NA, 1, 2)
			plt_df[,x:=pmin(as.numeric(x), 1.15)]
			plt_df[,y:=(y*0.2) + 0.2]
		}else{
			plt_df = data.table("state" = 0:8, "cols" = CN_palette, plot_labs, plot_shapes)
			pairings = data.table("col1" = c(0, 1, 3, 5, 7, NA), "col2" = c(NA, NA, 2, 4, 6, 8))
			pairings = pairings[col1 %in% plt_df$state | col2 %in% plt_df$state]
			matches = pairings[,lapply(.SD, function(x)x %in% states)]
			keep_states = unique(unlist(pairings[rowSums(matches) >= 1]))
			plt_df[,y:=which(pairings$col1 %in% state | pairings$col2 %in% state), by = 1:length(state)]
			plt_df[,x:=as.numeric(which(unlist(pairings[y,]) %in% state)), by = 1:length(state)]
			plt_df[,x:=pmin(as.numeric(x), 1.15)]
			plt_df[,y:=(y*0.2) + 0.2]
		}


		#    plt_df = plt_df[state %in% states]

		#plt_df = plt_df[state %in% keep_states]


		col_vec = plt_df$cols
		names(col_vec) = plt_df$state
		label_df = data.table(y = unique(plt_df$y), x = 0.83)
		cn_tbl = list("0" = 0, "1" = 1, "2" = 2, "3" = 2, "4" = 3, "5" = 3, "6" = 4, "7" = 4, "8" = 5)
		label_df$label = as.character(unique(unlist(cn_tbl[as.character(plt_df$state)])))
		if(5 %in% label_df$label){
			label_df[label == 5]$label = "5+"
		}
		#label_df[,label:=paste0(" ", label)]
		label_df[,just:=0]
		top_labs = data.table("y" = max(label_df$y, na.rm = T) + max(as.numeric(na.omit(c(diff(plt_df$y), 0.2)))), "x" = c(0.85, 1, 1.15), "label" = c("CN", "HOM", "HET"), "just" = 0.5)
		label_df = rbind(label_df, top_labs)
		legend_obj = ggplot(plt_df[state != 0], aes(x = x, y = y, color = cols, shape = plot_shapes)) + 
			geom_point(size = 5) + 
			theme_void() + 
			scale_x_continuous(limits = c(0.75, max(plt_df$x, na.rm = T) + 0.1)) + 
			scale_color_identity() +
			theme(legend.position = "none") + 
			geom_text(data = label_df, mapping = aes(x = x, y = y, label = label, hjust = just), 
								inherit.aes = F, 
								size = 3, 
								fontface = "bold") + 
			scale_shape_identity()
		if(colors != "loh"){
			legend_obj = legend_obj + 
				geom_point(plt_df[state == 0], mapping = aes(x = x, y = y, color = cols, shape = plot_shapes, stroke = 2)) + 
				scale_y_continuous(limits = c(0, max(plt_df$y, na.rm = T) + 2.5))
		}else{
			legend_obj = legend_obj + scale_y_continuous(limits = c(-2, max(plt_df$y, na.rm = T) + 2.4))
		}
		return(legend_obj)
	}
	
	legend = legend_plot_fn(plot_calls$state, colors = colors)
	
	cna_plot_fn = function(cnv_calls, colors = colors){
		chr = cnv_calls$chr[1]
		map_pos = lapply(names(plt_positions), function(k){
			cn = k
			k = plt_positions[k]
			poss_states = unique(states[state_cn == cn]$state)
			#if(!any(poss_states %in% x$state)){
			#  return(NA)
			#}
			if(length(poss_states) == 1){
				return(poss_states)
			}
			poss_states = cnv_calls[state %in% poss_states][,.(.N), by = state][order(N, decreasing = T)]$state[1]
			return(poss_states)
		})
		map_pos = unlist(map_pos)
		lines = data.table(state = as.character(map_pos), y = plt_positions)
		lines = lines[!is.na(state)]
		lines = rbind(lines, data.table("state" = 8, y = maxbound))
		ideal_ord = cnv_calls[,.N, by = state]
		cnv_calls$state = factor(cnv_calls$state, levels = ideal_ord$state)
		cnv_calls = cnv_calls[order(state)]
		cnv_calls$state = as.numeric(as.character(cnv_calls$state))
		cnv_calls$size = 1
		cnv_calls[state == 0, size:=2]
		clipped = cnv_calls[underflow == T | overflow == T]
		clipped[,size:=2]
		if(colors == "cnv"){
			plot_obj = ggplot(data = cnv_calls[end < cytobands$pos[cytobands$chr %in% chr][1] | pos > cytobands$end[which(cytobands$chr %in% chr)[2]]][state != 0], 
												aes(x = pos/1e+06, y = log(corrected_depth, base = 2), color = as.character(state), size = size)) + 
				geom_point_rast(aes(shape = factor(state)), alpha = 0.5) +
				geom_point_rast(cnv_calls[state == 0], mapping = aes(x = pos/1e+06, y = log(corrected_depth, base = 2), shape = factor(state)), stroke = 1) + 
				geom_point_rast(clipped[overflow == T], mapping = aes(x = pos/1e+06, y = log(corrected_depth, base = 2), color = as.character(state)), shape = 8) +
				scale_y_continuous(sec.axis = sec_axis(trans=~.*1, breaks = lines$y, labels = c(states[1:9,][state %in% lines$state]$state_cn, expression("">=20)), name = "Copy Number"), limits = c(minbound, maxbound)) + 
				scale_x_continuous(limits = c(min(cnv_calls$pos)/1e+06, max(cnv_calls$end)/1e+06)) +
				scale_size_identity() +
				#geom_rug(data = x[segment != shift(segment)], mapping = aes(x = pos/1e+06), inherit.aes = F) + 
				geom_hline(data = lines, mapping = aes(yintercept = y, color = state), size = 0.8, linetype = 1, alpha = 0.5) +
				scale_shape_manual(values = plot_shapes, 
													 labels = plot_labs,
													 name = "State") +
				scale_color_manual(name = "State",
													 values = CN_palette, 
													 labels = plot_labs) + 
				ylab("log(Read Depth)") + 
				xlab("Position (Mb)") + 
				ggtitle(paste0("Chromosome ", chr, " copy number profile")) + 
				theme_bw() + 
				theme(legend.position = "none")
		}else{
			plot_obj = ggplot(data = cnv_calls[end < cytobands$pos[cytobands$chr %in% chr][1] | pos > cytobands$end[which(cytobands$chr %in% chr)[2]]][state != 0], 
						 aes(x = pos/1e+06, y = log(corrected_depth, base = 2), color = zygosity, size = size)) + 
				geom_point_rast(aes(shape = factor(state)), alpha = 0.5) +
				geom_point_rast(cnv_calls[state == 0], mapping = aes(x = pos/1e+06, y = log(corrected_depth, base = 2), shape = factor(state)), stroke = 1) + 
				geom_point_rast(clipped[overflow == T], mapping = aes(x = pos/1e+06, y = log(corrected_depth, base = 2), color = zygosity), shape = 8) +
				scale_y_continuous(sec.axis = sec_axis(trans=~.*1, breaks = lines$y, labels = c(states[1:9,][state %in% lines$state]$state_cn, expression("">=20)), name = "Copy Number"), limits = c(minbound, maxbound)) + 
				scale_x_continuous(limits = c(min(cnv_calls$pos)/1e+06, max(cnv_calls$end)/1e+06)) +
				scale_size_identity() +
				#geom_rug(data = x[segment != shift(segment)], mapping = aes(x = pos/1e+06), inherit.aes = F) + 
				#geom_hline(data = lines, mapping = aes(yintercept = y, color = state), size = 0.8, linetype = 1, alpha = 0.5) +
				scale_shape_manual(values = plot_shapes, 
													 labels = plot_labs,
													 name = "State") +
				scale_color_manual(name = "Zygosity",
													 values = CN_palette, 
													 labels = plot_labs) + 
				ylab("log(Read Depth)") + 
				xlab("Position (Mb)") + 
				ggtitle(paste0("Chromosome ", chr, " copy number profile")) + 
				theme_bw() + 
				theme(legend.position = "none")
		}
		return(plot_obj)
	}
	
	vaf_plot_fn = function(cnv_calls, colors = colors){
		chr = cnv_calls$chr[1]
		cnv_calls[,size:=1]
		cnv_calls[state == 0, size:=2]
		ggplot(data = cnv_calls[end < cytobands$pos[cytobands$chr %in% chr][1] | pos > cytobands$end[which(cytobands$chr %in% chr)[2]]][state != 0][!is.na(maf)], 
					 aes(x = pos/1e+06, y = unlist(unmerge_mafs_grouped(maf, flip = T)), color = factor(state), size = size)) + 
			geom_point_rast(size = 1, alpha = 0.5, aes(shape = factor(state))) + 
			geom_point_rast(cnv_calls[state == 0], mapping = aes(shape = factor(state)), stroke = 1) +
			scale_size_identity() +
			scale_shape_manual(values = plot_shapes, 
												 labels = plot_labs,
												 name = "State") +
			scale_color_manual(name = "State",
												 values = CN_palette, 
												 labels = plot_labs) + 
			ylab("Major allele frequency") + 
			xlab("Position (Mb)") + 
			ggtitle(paste0("Chromosome ", chr, " allele frequency profile")) + 
			scale_y_continuous(limits = c(0.5, 1)) +
			theme_bw() +
			theme(legend.position = "none",
						plot.margin = unit(c(5.5,43,5.5,5.5), "pt"))
	}
	
	### Ideograms
	## TODO: take in genome ver and use specific one
	cyto_plot_fn = function(cnv_calls, text = F){
		cytoband_dat = data.table::copy(cytobands)
		ideo_colors = c("gneg" = "white", "gpos25" = "black", "gpos50" = "black", "gpos75" = "black", "gpos100" = "black", "acen" = "red", "gvar" = "black", "stalk" = "black")
		ideo_alphas = c("gneg" = 0, "gpos25" = .25, "gpos50" = .5, "gpos75" = .75, "gpos100" = 1, "acen" = 1, "gvar" = 1, "stalk" = 0.5)
		chr = cnv_calls$chr[1]
		xlim = max(cnv_calls$end)
		xmlim = min(cnv_calls$pos)
		ylim = 1
		if(grepl("chr", cytoband_dat$chr[1])){
			cytoband_dat[,chr:=gsub("chr", "", chr)]
		}
		## used because chr variable name is same as the variable 
		ideo_1 = cytoband_dat[eval(cytoband_dat[, chr %in% ..chr])]
		out_plt = ggplot(ideo_1, aes(xmin = pos, xmax = end, ymin = 0, ymax = ylim, fill = stain, alpha = stain)) + 
			geom_rect(color = "black") + 
			geom_rect(aes(xmin = min(pos), xmax = max(end), ymin = 0, ymax = ylim), inherit.aes = F, fill = NA) +
			scale_fill_manual(values = ideo_colors) +
			scale_alpha_manual(values = ideo_alphas) +
			scale_x_continuous(limits = c(0, xlim), position = "top", breaks = seq(from = 0, to = 2.5e+08, by = 2.5e+07), labels = c(0, paste0(seq(25, 250, by = 25), "Mb"))) +
			scale_y_continuous(limits = c(-1, 1)) +
			theme_void() + 
			theme(legend.position = "none", plot.title = element_text(size = 10), plot.margin = unit(c(5.5,30,5.5,5.5), "pt")) #+ 
		if(text == T){
			out_plt = ggplot(ideo_1, aes(xmin = pos, xmax = end, ymin = 0, ymax = ylim, fill = stain, alpha = stain)) + 
				geom_text(aes(x = (pos + end)/2, y = -0.5, label = lab), inherit.aes = F, angle = 90, vjust = 0.5, check_overlap = T) +
				geom_rect(color = "black") + 
				geom_rect(aes(xmin = pmax(min(pos), xmlim), xmax = pmin(max(end), xlim), ymin = 0, ymax = ylim), inherit.aes = F, fill = NA, alpha = 1, color = "black") +
				
				scale_fill_manual(values = ideo_colors) +
				scale_alpha_manual(values = ideo_alphas) +
				scale_x_continuous(limits = c(0, xlim), position = "top", breaks = seq(from = 0, to = 2.5e+08, by = 2.5e+07), labels = c(0, paste0(seq(25, 250, by = 25), "Mb"))) +
				scale_y_continuous(sec.axis = sec_axis(trans=~.*1, breaks = NULL, labels = NULL, name = "Copy Number"),
													 limits = c(-1, 1)) +
				theme_void() + 
				theme(legend.position = "none", plot.title = element_text(size = 10), plot.margin = unit(c(5.5,30,5.5,5.5), "pt"))
				
		}
		return(out_plt)
	}
	if(length(cytobands) == 1 & all(cytobands == F)){
		if(mode == "all"){
			CNA_plot = plot_grid(plot_grid(cna_plot_fn(cnv_data, colors = colors), vaf_plot_fn(cnv_data, colors = colors), align = "v", axis = "l", ncol = 1, rel_heights = c(1, 0.5)), legend_plot_fn(plot_calls$state, colors = colors), rel_widths = c(1, 0.2))
		}else if(mode == "zoomed"){
			CNA_plot = plot_grid(plot_grid(cna_plot_fn(cnv_data, colors = colors), align = "v", axis = "l", ncol = 1, rel_heights = c(1, 0.5)), legend_plot_fn(plot_calls$state, colors = colors), rel_widths = c(1, 0.2))
		}
	}else{
		if(mode == "all"){
			CNA_plot = plot_grid(plot_grid(cna_plot_fn(cnv_data, colors = colors), vaf_plot_fn(cnv_data, colors = colors), cyto_plot_fn(cnv_data), align = "v", axis = "l", ncol = 1, rel_heights = c(1, 0.5, 0.05)), legend_plot_fn(plot_calls$state, colors = colors), rel_widths = c(1, 0.2))
		}else if(mode == "zoomed"){
			CNA_plot = plot_grid(plot_grid(cna_plot_fn(cnv_data, colors = colors), cyto_plot_fn(cnv_data, text = T), align = "v", axis = "l", ncol = 1, rel_heights = c(1, 0.5, 0.05)), legend_plot_fn(plot_calls$state, colors = colors), rel_widths = c(1, 0.2))
		}
	}

	return(CNA_plot)
}

#' @export
focus_view = function(cnv_calls, chrom, start_position, end_position, colors, cytobands){
	filt_data = cnv_calls[chrom == chr]
	filt_data = filt_data[chrom == chr & pos >= start_position & end <= end_position]
	print(chr)
	print(filt_data)
	plot_ploidetect(cnv_data = filt_data, cn_positions = cn_positions, cytobands = cytobands, mode = "zoomed", colors = colors)
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
