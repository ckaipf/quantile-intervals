#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.tag = "example_plot"

// Files
params.genome = "example/genome_file"
params.variable_bedgraph = ["example/example_rep_2.bg"]
//params.variable_bedgraph = ["example/example_rep_1.bg", "example/example_rep_2.bg"] // If multiple bedgraphs are provided, they can be merged by a function, e.g. addition
params.target = "example/example.gff"

// Bedtools slop parameters
params.slop_left = 400
params.slop_right = 400

// R parameters
params.normalization_mode = "function" // function, position, z_score
params.normalize_function = "(function(x)  max(x, na.rm = T))"// Use constant function to ignore normalization (function(x) 1)
params.normalize_position = -50
params.quantiles = "seq(0,1,0.1)"
params.groups = 'c("group")'
params.facet_formula = "group ~ ."

workflow  {
  main:
  quantilize(params.genome, params.variable_bedgraph, params.target)
}


workflow quantilize {
  take: genome
        variable_bedgraph
        target
  main:
  variable_bedgraph_ch = Channel.fromPath(variable_bedgraph) | \
  sortBed | \
  toList //| \
  //mergeReplicates

  genome_file = Channel.fromPath(genome) 
  target_ch = Channel.fromPath(target) | \
      sortGff
  slopped = slop(target_ch, genome_file) 
  intersect(slopped, variable_bedgraph_ch) |
  maultaschify 

  maultaschify.out.intervals | removeLastLine | plotQuantiles
  maultaschify.out.warnings.toList().forEach { if(it) print it.trim() }
}

process intersect {
  //  publishDir ".", mode: "copy"
    input:
    path slop
    path bedGraph
    output:
    path "*.flanked"
   
    """
    bedtools intersect -a "${slop}" -b ${bedGraph} -wb -s -wa | \
    awk -v OFS="," '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$11,\$12,\$14}' > "${slop}.flanked"
    
    # Add dummy line at the end for subsequent process
    echo "-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1" >> "${slop}.flanked"
    """
}

process maultaschify {
    input:
    path csv
    output:
    path "intervals.tsv", emit: intervals
      stdout emit: warnings
    
    """
    #!/usr/bin/env python3
    import os, csv
    # Takes a list of intervals with assoc. values and returns a 
    # flattend and 1-binned list of values 

    def flat_n_chop_n_fill(li, acc=[]):
      if len(li) > 1:
        (i,j,x), (p,q,y),*tail = li
        # Problem when variable interval exceeds target interval 
        if p < i:
           p = i
        if j > q:
          j = q
        return flat_n_chop_n_fill([(p,q,y)] + tail, acc+([x]*(j-i))+([0]*(p-j)))
      else:
        return acc
    def strList(li):
      return [str(el) for el in li]

    # Assumes input is sorted
    # Collect intersected intervals of target and variable, then flat chop and fill
    with open('${csv}', mode='r') as in_, open('intervals.tsv', mode="w") as out:
      l, pre, wr, rows = [], [-1]*12, csv.writer(out, delimiter="\t"), csv.reader(in_)
      for row in rows:
        target, variable = row[:9], (int(row[9]),int(row[10]),float(row[11]))
        if target != pre:
          wr.writerow(pre+[",".join(strList(flat_n_chop_n_fill([(int(pre[3]),int(pre[3]),0)]+l+[(int(pre[4]),int(pre[4]),0)])))]) 
          pre, l = target, []
        l += [variable]
    """
  }

  process removeLastLine {
    input:
    path file
    output:
    path "${file}.tail"
    """
    tail -n +2 ${file} > ${file}.tail 
    """
  }


process mergeReplicates {
 // publishDir ".", mode: "copy"
  input:
  val files
  output:
  path("acc")

  
  script:
  if (files.size() > 1) {
    head = files[0]
    tail = files[1..-1]
    input_list = tail.toString().minus("[").replace(",", "").minus("]") 
  }
   
  if (files.size() > 1) 
    """
    cat ${head} > acc
    
    for i in ${input_list}; do
      paste acc <(cut -f 5 \${i}) | awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,\$4,\$5+\$7,\$6}' > tmp 
      mv tmp acc
    done
    """
  
  else 
    """
    cat ${files[0]} > acc
    """
}


process sortBed {
  input:
  path(a)
  output:
  path("${a}.sorted")

  """
  sort -k1,1 -k2,2n -k3,3n ${a} > ${a}.sorted
  """
}



process sortGff {
  input:
  path a
  output:
  path "${a}.sorted"

  """
  sort -k1,1 -k4,4n -k5,5n ${a} > ${a}.sorted
  """
}


process slop {
//  publishDir ".", mode: "copy"

  input:
  path target
  path genome

  output:
  path "*.slop"

  """
  bedtools slop -l ${params.slop_left} -r ${params.slop_right} -g ${genome} -i ${target} > ${target}.slop
  """
}

process plotQuantiles {
    publishDir ".", mode: "copy"
    input:
    path csv
    output:
    path "*.png"
   
  """
  #!/usr/bin/env Rscript 
  library(tidyverse)

  groups <- ${params.groups}
  id <- c("start", "end", "strand")
  normalization_mode <- "${params.normalization_mode}"   

  #---
  # Set quantiles
  quantiles <- ${params.quantiles}
  quantiles_functions <- lapply(quantiles, function(p) {
    f <- function(x) quantile(x, probs = p, na.rm = T)
    return(f)
  }) 
  quantiles_functions <- quantiles_functions %>% 
    set_names((quantiles * 100))

  #---
  x <- read_delim("${csv}", delim = "\t",
             col_types = "cccddccccc", 
             col_names = c("seqname",
                           "source", 
                           "feature", 
                           "start", 
                           "end", 
                           "score", 
                           "strand", 
                           "frame", 
                           "attribute",
                           "target_variable"
                           )
             )

  cs <- str_count(x\$attribute, ";")
  ps_max <- max(str_count(x\$target_variable, ",")) +1

str_split(x\$attribute, ";") %>%
     map(function(s) map(groups, function(g) (str_subset(s, pattern = paste0("^", g))))) %>%
     map(set_names, nm = groups) %>%
     bind_rows() %>% 
     mutate(across(everything(), .fns =  ~ str_remove(., ".*(?=[=])") %>% str_remove("=") %>% str_remove(";"))) %>%
     bind_cols(x) %>% 
     mutate(l = str_count(target_variable, ",") + 1) %>% 
     mutate(normalization_mode = normalization_mode) %>%
  #---
  # Calculate sample sizes for groups
    group_by(across(all_of(groups))) %>% 
    mutate(n = n_distinct(across(everything()))) %>%
    ungroup() %>% 

  #---
  # Get positions variables
    separate(target_variable, 
             into = paste0("target_variable_", as.character(c(1:max(unlist(chuck(., "l")))))), 
             sep = ",", 
             convert = T
             ) %>%  
      pivot_longer(starts_with("target_variable_"), names_to = "position") %>% 
      mutate(position = as.numeric(str_remove(position, "target_variable_"))) %>% 
  # ---
  # Flip negative strand
       group_by(across(all_of(id))) %>%
       mutate(position = if_else(strand == "+", position, l-position)) %>% 

  # ---
  # Align windows by their centers (ambiguity +-1bp, floor, ceiling)  
       filter(!is.na(value)) %>%
       mutate(max_l = max(l)) %>%
       group_by(across(c(-position, value))) %>%
       mutate(position = as.numeric(position))  %>% 
       mutate(position = position + floor((max_l-l)/2))  %>% 
       mutate(position = position - (max_l/2)) %>%
  #---
  # Pivot_wider does not support margins, reshape2::dcast also adds cols if values don't exists for variable   
       mutate(position = paste0("target_variable_", position)) %>%
       reshape2::dcast(... ~ position, value.var = "value")  %>%
       as_tibble() %>% 
  #---

  #---
  # Calculate parameters for normalization within window
    group_by(., across(all_of(id))) %>%
    
     { if (normalization_mode == "position") {
         mutate(., parameters =  list(list(
         p_1 =.data[[ paste0("target_variable_", "${params.normalize_position}" )]], 
         p_2 = NA
         )))
     }
       else if (normalization_mode == "function") {
         rowwise(.) %>% 
           mutate(parameters = list(list(
             p_1 = ${params.normalize_function}(c_across(starts_with("target_variable_"))), 
             p_2 = NA
             )))
       }
         else if (normalization_mode == "z_score") {
           rowwise(.) %>% 
             mutate(
               parameters = list(
                 list(
                p_1 = mean(c_across(starts_with("target_variable_")), na.rm = T),
                p_2 = sd(c_across(starts_with("target_variable_")), na.rm = T)
                )
               )
             ) %>% ungroup
         } else {.}
     } %>% 
    pivot_longer(starts_with("target_variable_"), names_to = "position") %>% 
    mutate(position = as.numeric(str_remove(position, "target_variable_"))) %>% 


  # ---
  # Apply normalization within window
    group_by(across(all_of(id))) %>% 
    
       mutate(
      x_norm = 
        case_when(normalization_mode == "function" ~ value / (chuck(parameters, 1, "p_1")),
                  normalization_mode == "position" ~ value / (chuck(parameters, 1, "p_1")),
                  normalization_mode == "z_score" ~ ((as.numeric(value) - chuck(parameters, 1, "p_1")) / (chuck(parameters, 1, "p_2"))) 
            )
    ) %>% 
       
  # ---
  # Reduce and tidy data
    filter(!is.infinite(x_norm), !is.na(x_norm)) %>% 
    select(all_of(c(groups, id, "position", "n", "x_norm"))) %>% 
    distinct(across(everything())) %>% 
    
  # ---
  # Summarize variable for each positions by quantiles
    group_by(position, n, across(all_of(groups))) %>% 
    summarise(across(.cols = x_norm, .fns = quantiles_functions, .names = "m_{.fn}")) %>% 
    pivot_longer(starts_with("m_")) %>% 
    mutate(quantile = as.factor(as.numeric(str_remove(name, "m_")))) %>%

  # ---
  # Cosmetics
    mutate(n = paste0("italic(n)==", n)) %>% 
    mutate(quantile = fct_reorder(quantile, as.numeric(quantile), .desc = T) ) %>%

  # ---
  # Plot
    ggplot(aes(x = position, y = value, fill = quantile)) +
    #geom_area(position = "identity", alpha = .9) +
    stat_smooth(
    geom = 'area',
    method = 'loess', 
    span = .01,
    alpha = .5
    ) + 

    scale_fill_viridis_d(option = "B", labels = ~ paste0(., "%")) +
    geom_vline(xintercept = 0) +
       
    ylab(bquote(paste(italic(v)[normalized]))) +
    
    facet_grid(formula("${params.facet_formula}"), switch = "y") + 
    geom_text(parse = T, mapping = aes(
      x = -Inf, 
      y = Inf, 
      label = n,
      hjust = -0.1, 
      vjust = 1.5
    ), check_overlap = T, show.legend = F) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(
      strip.background.x = element_rect(color = "white"),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      strip.placement = "outside",
      legend.position = "bottom"
  ) -> p

  ggsave("${params.tag}.png", plot = p, device = "png")
  """
}
