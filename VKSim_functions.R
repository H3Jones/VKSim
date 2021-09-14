#VKSim_functions

#helper functions
check_CHO_file<-function(file_path){
  if(all(colnames(readr::read_csv(file_path, n_max = 5)) %in% c("C", "H", "O"))) {
    return(NULL)
  } else {
    return("Ensure file is the correct type")
  }
}

check_Composer64<-function(file_path){
  if(!is.na(KairosMSfunctions::determine_file_format(file_path))) {
    return(NULL)
  } else {
    return("Ensure file is the correct type")
  }
}

check_file_input<-function(file_path, file_type){
  switch(file_type,
         "Composer64" = check_Composer64(file_path),
         "CHO" =  check_CHO_file(file_path),
         "Ensure file is the correct type"
  )
  
}

string_to_list<-function(string_vec){
  if(any(is.na(as.numeric(string_vec[1:3])))) return(NULL)
  
  raw_rep_value <-  as.numeric(string_vec[4])
  rep_value <- if (is.na(raw_rep_value) | raw_rep_value <= 0) {
    1
  } else {
    raw_rep_value
  }
  
  vec <- as.numeric(string_vec[1:3])
  list(
    raw_vec = vec,
    rep_value = rep_value
  )
}

extract_inputs<-function(string_list){
  vec <- set_names(string_list$raw_vec, c("C","H","O"))
  out <- tibble(!!!vec)
  out$rep <- string_list$rep_value
  out$reaction <- createReaction(c("C","H","O"), vec)
  return(out)
}

tidy_reaction <- function(starting_data, reactions, min_dbe = 0){
  
  input_data <- starting_data %>%
    transpose() %>%
    simplify_all()
  
  output_data <- generateReactions(input_data, reactions, min_dbe)
  
  tidy_output <- map_df(flatten(output_data),~bind_rows(.x, .id = "ID"), .id = "Iter") %>%
    bind_rows(.,rowid_to_column(starting_data, var = "ID") %>% mutate(ID = as.character(ID), "Iter" = "0")) %>%
    mutate(H_C = H/C,
           O_C = O/C,
           Iter = as.numeric(Iter),
           ID = as.numeric(ID))
  
  return(tidy_output)
  
  
}

filter_iter <-function(sim_data){
  max_Iter <- max(sim_data$Iter)
  if(max_Iter > 10){
    data <- sim_data %>%
      filter(Iter %in% seq(0,max_Iter,round(max_Iter/10)))
  } else {
    data <- sim_data    
  }
  return(data)
}

#Custom ggplot2 theme
theme_custom <- function (
  base_size = 11, 
  base_family = "",
  base_line_size = base_size/22, 
  base_rect_size = base_size/22
) {
  half_line <- base_size/2
  theme_bw(base_size = base_size,
           base_family = base_family, 
           base_line_size = base_line_size,
           base_rect_size = base_rect_size) %+replace% 
    theme(
      # legend.position = "none",
      # legend.title=element_blank(),
      strip.background =element_rect(fill= "#4169E1"),
      strip.text = element_text(colour = "white", 
                                size = rel(0.8),
                                face = "bold",
                                margin = margin(0.8 * half_line,
                                                0.8 * half_line,
                                                0.8 * half_line,
                                                0.8 * half_line)),
      # axis.text = element_text(size = 10),
      axis.title.y = element_text(face = "bold", angle = 90, vjust = 2),
      axis.title.x = element_text(face = "bold"))}




multi_save_plot<-function(file,plot,type,width,height,units,dpi){
  switch(type,
         "png" = ggsave(filename = file,plot = plot,width = width,height = height,units = units, device = ragg::agg_png, res = dpi),
         "pdf" = ggsave(filename = file,plot = plot,width = width,height = height,units = units, device = grDevices::cairo_pdf, dpi = dpi),
         "eps" = ggsave(filename = file,plot = plot,width = width,height = height,units = units, device = grDevices::cairo_ps,dpi = dpi,fallback_resolution = 600),
         "tiff" = ggsave(filename = file,plot = plot,width = width,height = height,units = units, device = ragg::agg_tiff,res = dpi),
         "svg" = ggsave(filename = file,plot = plot,width = width,height = height,units = units, dpi = dpi)
  )
  
}

#wrapper for plot type selector
select_plot_file_type<-function(inputId, label = "Select the file type"){
  radioButtons(inputId = inputId,label = label, choices = list("png", "pdf", "eps", "tiff","svg"), selected = "png")
}

