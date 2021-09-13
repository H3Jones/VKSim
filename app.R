

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
    
    generateReactions(input_data, reactions, min_dbe)
    
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

# Packages to Load --------------------------------------------------------
suppressPackageStartupMessages(c(
    library(shiny),
    library(shinydashboard),
    library(shinyWidgets),
    library(magrittr),
    library(tidyverse),
    library(glue),
    library(plotly),
    library(gganimate)
))

#Rcpp::sourceCpp("VKSim_cpp_funs.cpp")
Rcpp::sourceCpp("VKSim_improved.cpp")

theme_set(theme_custom())

# Sidebar -----------------------------------------------------------------

sidebar <- dashboardSidebar(
    width = 250,
    hr(),
    sidebarMenu(id="tabs",
                menuItem("Input",tabName = "Input",selected = TRUE),
                menuItem("Simulate",tabName = "Sim"),
                menuItem("Other plots",tabName = "Plots"),
                menuItem("Gif",tabName = "Gif"),
                menuItem("Comparison", tabName = "Comparison"),
                hr(),
                    materialSwitch("set_limits","Set plot limits"),
                    conditionalPanel(condition = 'input.set_limits',
                                     numericRangeInput("set_limits_x", "Set O/C axis limits", value = c(0,3)),
                                     numericRangeInput("set_limits_y", "Set H/C axis limits", value = c(0,3))
                                     ),
                hr(),
                fluidRow(
                    column(width = 3, offset = 3,
                h5("Graphics options")
                    )
                ),
                radioButtons(inputId = "unitsGraph", label = "Select the units", choices = list("mm", "cm","in"), selected = "in"),
                numericInput(inputId = "dpiGraph",label = "DPI of the graphics:", value = 600, max = 2000),
                numericInput(inputId = "widthGraph",label = "Width of the graphics:", value = 12),
                numericInput(inputId = "heightGraph",label = "Height of the graphics:", value = 8),
                numericInput(inputId = "textSize",label = "Text size", value = 16),
                select_plot_file_type(inputId = "plot_format",label = "Select the file type")
    ))


# Body --------------------------------------------------------------------

body <- dashboardBody(
    tabItems(
        
        # Main Tab ----------------------------------------------------------------
        tabItem(tabName = "Input",
                fluidPage(
                    box(width = NULL,
                        selectInput("file_type","Select input type",choices = c("Composer64","C,H,O values" = "CHO")),
                        fileInput(inputId = "file_input",
                                  label = "Upload data",
                                  multiple= TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        selectizeInput("masslist_selector","Choose an uploaded masslist to use", choices = "")
                        ),
                    box(width = NULL,
                        materialSwitch("filter_mass","Filter by m/z"),
                        conditionalPanel(condition = "input.filter_mass",
                                         numericRangeInput("filter_mz","Filter m/z range", value = c(0,1500))
                        ),
                        materialSwitch("filter_data","Filter Data"),
                        conditionalPanel(condition = "input.filter_data",
                                         numericRangeInput("filter_carbon","Filter Carbon Number range", value = c(0,30)),
                                         numericRangeInput("filter_dbe","Filter DBE range", value = c(0,10)),
                                         numericRangeInput("filter_HC","Filter H/C range", value = c(0,4)),
                                         numericRangeInput("filter_OC","Filter O/C range", value = c(0,2))
                        ),
                        div(style = 'overflow-x: scroll', DT::dataTableOutput("neutral_data")),
                        collapsible = TRUE, title = "Neutral Data", status = "primary", solidHeader = TRUE
                        )
                )
                ),
        tabItem(tabName = "Sim",
                fluidPage(
                    tags$head(
                        tags$style(HTML("hr {border-top: 1px solid #000000;}"))
                    ),
                    box(width = NULL,
                        fluidRow(
                        column(
                            width = 6,
                            textAreaInput("list_input", "Reaction Vector Input",
                                          placeholder = "-1,0,-2 (6) #6*CO2 losses \n")
                        ),
                        column(
                            width = 6,
                            div(style = 'overflow-x: scroll', DT::dataTableOutput("reaction_vectors"))
                        )
                        ),
                        fluidRow(
                            column(
                                width = 12,
                                actionButton("simulate","Simulate Reaction"),
                                hr(),
                                textInput("sim_names","ID to use", width = '30%'),
                                actionButton("save_sim","Save simulation data for comparison")
                            )
                        )
                        ),
                    box(width = NULL,
                        downloadButton("download_vk_plot","Download Plot"),
                        plotlyOutput("vk_plot", height = "600px"),
                        div(style = 'overflow-x: scroll', DT::dataTableOutput("vk_datatable")),
                        collapsible = TRUE, title = "Van Krevelen Plot", status = "primary", solidHeader = TRUE
                        ),
                    box(width = NULL,
                        downloadButton("download_HC_plot","Download H/C Plot"),
                        downloadButton("download_OC_plot","Download O/C Plot"),
                        plotOutput("HC_plot", height = "300px"),
                        plotOutput("OC_plot", height = "300px"),
                        collapsible = TRUE, title = "Density Plot", status = "primary", solidHeader = TRUE
                        )
                    
                )
        ),
        tabItem(tabName = "Plots",
                fluidPage(
                    box(width = NULL,
                        column(width = 12,
                               div(style = 'overflow-x: scroll', DT::dataTableOutput("iteration_details"))
                               ),
                        title = "Iteration details",
                        status = "primary", solidHeader = TRUE
                        ),
                    box(width = NULL,
                        column(width = 12,
                               downloadButton("download_iter_plot","Download Plot"),
                               plotlyOutput("iter_plot", height = "600px")
                               ),
                        title = "Iteration plot",
                        status = "primary", solidHeader = TRUE
                        ),
                    box(width = NULL,
                        column(width = 12,
                               downloadButton("download_class_plot","Download Plot"),
                               plotlyOutput("class_plot", height = "600px")
                               ),
                        title = "Class distribution",
                        status = "primary", solidHeader = TRUE
                        ),
                    box(width = NULL,
                        column(width = 12,
                               downloadButton("download_violin_plot","Download Plot"),
                               materialSwitch("violin_filter_classes","Filter classes"),
                               conditionalPanel(condition = 'input.violin_filter_classes',
                                                selectizeInput("violin_class_filter","Select classes", choices = c("HC", "O1"), multiple = TRUE)
                                                ),
                               plotOutput("violin_plot", height = "600px")
                        ),
                        title = "Violin plot",
                        status = "primary", solidHeader = TRUE
                    ),
                    box(width = NULL,
                        column(width = 12,
                               downloadButton("download_HC_val_data","Download Data"),
                               div(style = 'overflow-x: scroll', DT::dataTableOutput("HC_datatable")),
                               downloadButton("download_HC_val_plot","Download Plot"),
                               splitLayout(
                                 numericInput("HC_val_line_size","Point size", value = 2, width = '150px'),
                                 numericInput("HC_val_point_size","Point size", value = 3, width = '150px')
                               ),
                               plotOutput("HC_val_plot", height = "600px")
                        ),
                        title = "Hydrocarbon H/C values",
                        status = "primary", solidHeader = TRUE
                    ),
                )),
        tabItem(tabName = "Gif",
                fluidPage(
                    imageOutput("vk_gif")
                )),
        tabItem(tabName = "Comparison",
                fluidPage(
                    box(width = NULL,
                        tableOutput("stored_data"),
                        pickerInput("selected_entries","Select data to use",choices = "", multiple = TRUE),
                        collapsible = TRUE, title = "Comparison Options", status = "primary", solidHeader = TRUE
                        ),
                    box(width = NULL,
                        downloadButton("download_combined_vk_plot","Download Plot"),
                        plotOutput("combined_vk_plot", height = "600px"),
                        downloadButton("download_combined_class_plot","Download Plot"),
                        plotlyOutput("combined_class_plot", height = "600px"),
                        materialSwitch("combined_violin_filter_classes","Filter classes"),
                        conditionalPanel(condition = 'input.combined_violin_filter_classes',
                                         selectizeInput("combined_violin_class_filter","Select classes", choices = c("HC", "O1"), multiple = TRUE)
                        ),
                        downloadButton("download_combined_violin_plot","Download Plot"),
                        plotOutput("combined_violin_plot", height = "600px"),
                        downloadButton("download_combined_HC_val_plot","Download Plot"),
                        plotOutput("combined_HC_val_plot", height = "600px"),
                        collapsible = TRUE, title = "Comparison Plots", status = "primary", solidHeader = TRUE
                    )
                    
                ))
        
        
        
        
        
        
        # Closing brackets --------------------------------------------------------
        
    ))

# Setup UI ----------------------------------------------------------------

ui <- dashboardPage(
    dashboardHeader(title = "VKSim",
                    titleWidth = 250),
    sidebar,
    body
)

# Setup Server ------------------------------------------------------------

server <- function(input, output, session) {
    
    reactive_values <- reactiveValues()  
   
    # raw_masslist<-reactive({
    #     req(!is.null(input$file_input))
    #     validate(
    #         check_file_input(file_path = input$file_input$datapath,
    #                          file_type = input$file_type)
    #     )
    #     switch(input$file_type,
    #                     "Composer64" = KairosMSfunctions::Master_ReadeR(file_to_read = input$file_input$datapath) %>%
    #                       filter(.,isotope == FALSE) %>%
    #                       filter(str_detect(class,"N|S|B|Na", negate = TRUE)) %>%
    #                       mutate(electron_status = ifelse(dbe != as.integer(dbe),'even','odd')),
    #                     "CHO" = readr::read_csv(input$file_input$datapath, col_types = cols(
    #                       C = col_double(),
    #                       H = col_double(),
    #                       O = col_number()
    #                     )) %>%
    #                       mutate(.,
    #                              assigned.mz = C*12 + H + O*16,
    #                              dbe = C -H/2 +1
    #                       )
    #              )
    #     
    # })
    
    # output$debug <- renderText({
    #     input$list_input
    #     })
    
    observeEvent(input$file_input,{
      req(!is.null(input$file_input))
      validate(
        check_file_input(file_path = input$file_input$datapath,
                         file_type = input$file_type)
      )
      file<-switch(input$file_type,
             "Composer64" = KairosMSfunctions::Master_ReadeR(file_to_read = input$file_input$datapath) %>%
               filter(.,isotope == FALSE) %>%
               filter(str_detect(class,"N|S|B|Na", negate = TRUE)) %>%
               mutate(electron_status = ifelse(dbe != as.integer(dbe),'even','odd')),
             "CHO" = readr::read_csv(input$file_input$datapath, col_types = cols(
               C = col_double(),
               H = col_double(),
               O = col_number()
             )) %>%
               mutate(.,
                      assigned.mz = C*12 + H + O*16,
                      dbe = C -H/2 +1
               )
      )
      
      filename <- str_extract(input$file_input$name,'.+(?=.csv)')
      reactive_values$raw_masslists[[length(reactive_values$raw_masslists)+1]] <- file
      reactive_values$masslist_names[[length(reactive_values$masslist_names)+1]] <- filename

    }
                 )
    

# Reaction details input --------------------------------------------------

    
    raw_vector_input <-reactive({
        req(!is.null(input$list_input))
        
        str_split(input$list_input,"\n") %>%
            unlist() %>% 
            str_extract_all("(#?)(-?)[:digit:]+") %>%
            keep(~!any(str_detect(.x, "#"))) %>%
            keep(~length(.x) >= 3) %>%
            map(~string_to_list(.x)) %>%
            keep(~!is_null(.x))
        
    })
    
    reaction_details <- reactive({
        req(length(raw_vector_input()) > 0)
        map_df(raw_vector_input(), ~extract_inputs(.x)) 
    })
    
    output$reaction_vectors<-DT::renderDataTable({
        req(reaction_details())
        reaction_details()
        })
    
    

# Prepare data ------------------------------------------------------------

    observe({
      req(length(reactive_values$raw_masslists) > 0)
      updateSelectInput(session = shiny::getDefaultReactiveDomain(), inputId = "masslist_selector", label = "Choose an uploaded masslist to use", choices = unlist(reactive_values$masslist_names))
    })
    
    
    filtered_data <- reactive({
        req(length(reactive_values$raw_masslists) > 0)
      reactive_values$raw_masslists[[match(input$masslist_selector,reactive_values$masslist_names)]] %>%
            {if(input$filter_mass) filter(., between(assigned.mz,input$filter_mz[[1]],input$filter_mz[[2]])) else .} %>%
            {if(input$filter_data) filter(., 
                                          between(C,input$filter_carbon[[1]],input$filter_carbon[[2]]) &
                                              between(dbe,input$filter_dbe[[1]],input$filter_dbe[[2]]) &
                                              between(H/C,input$filter_HC[[1]],input$filter_HC[[2]]) &
                                              between(O/C,input$filter_OC[[1]],input$filter_OC[[2]])
                                              
                                          ) else .}
    })
    
    starting_data <- reactive({
        req(filtered_data())
        
        filtered_data() %>%
            select(C, H, O) %>%
            distinct(.keep_all = TRUE)
    })
    
    output$neutral_data<-DT::renderDataTable({
        req(filtered_data())
        filtered_data()
    })
    
    
    simulated_data <- eventReactive(input$simulate,{
        req(starting_data(), raw_vector_input())
        test_reaction <-  raw_vector_input()
        
        input_data <- starting_data() %>%
            transpose() %>%
            simplify_all()
        
        generateReactions(input_data , test_reaction, min_dbe = 0)
        
        
    })
    
    tidy_data <-reactive({
        req(simulated_data(),starting_data())
        
        map_df(flatten(simulated_data()),~bind_rows(.x, .id = "ID"), .id = "Iter") %>%
            bind_rows(.,rowid_to_column(starting_data(), var = "ID") %>% mutate(ID = as.character(ID), "Iter" = "0")) %>%
            mutate(H_C = round(H/C, digits = 6),
                   O_C = round(O/C, digits = 6),
                   Iter = as.numeric(Iter),
                   ID = as.numeric(ID)) 
    })
    

# Save simulation data ----------------------------------------------------

    observeEvent(input$save_sim,{
        req(tidy_data())
        
        sim_names <- if (input$sim_names != ""){
          input$sim_names
        }else{
          paste0("simulation", length(reactive_values$sim_data))
        }
        
        validate(
          need(!sim_names %in% reactive_values$sim_names, "Name has already been used")
        )
        
        reactive_values$sim_data[[length(reactive_values$sim_data)+1]] <- tidy_data()
        #reactive_values$reaction_data[[length(reactive_values$reaction_data)+1]] <- reaction_details
        reactive_values$sim_names[[length(reactive_values$sim_names)+1]] <- sim_names
        
        save_message <- glue("Data saved with name {sim_names},
                             currently {length(reactive_values$sim_data)} simulations stored")
        
        showNotification(ui = save_message)
    })
    
    observeEvent(input$clear_list,{
        reactive_values$sim_data <- NULL
        reactive_values$reaction_data <- NULL
        reactive_values$sim_names <- NULL
    })
    
      

# Iteration details -------------------------------------------------------

    
    iteration_details <- reactive({
        req(tidy_data(), reaction_details(), simulated_data())

        rxns <- reaction_details()$reaction
        rxn_length <- simulated_data() %>%
            map_int(length)

        reaction_details <- tibble(
            reaction = rep(rxns,rxn_length),
            Iter = 1:sum(rxn_length)
        ) %>%
            add_row(reaction = "Starting", Iter = 0)

        n_molecules <- nrow(filtered_data())

        tidy_data() %>%
            mutate(formula = paste0("C",C,"H",H,"O",O)) %>%
            arrange(Iter) %>%
            group_by(ID) %>%
            mutate(unchanged = lag(formula) == formula) %>%
            ungroup() %>%
            group_by(Iter) %>%
            summarise(unreacted = sum(unchanged), .groups = "drop_last") %>%
            left_join(., reaction_details, by = "Iter") %>%
            mutate(reacted = n_molecules - unreacted) 
    })
    
    output$iteration_details<-DT::renderDataTable({
        req(iteration_details())
        iteration_details() 
    })
    
    #Iteration plot
    iter_plot <- reactive({
        req(iteration_details())
        
        iteration_details() %>%
            pivot_longer(cols = c(unreacted, reacted)) %>%
            ggplot(aes(Iter, value, fill = paste(name, reaction))) +
            geom_col() +
            scale_fill_viridis_d(end = .8) +
            labs(
                y = "Number of molecules",
                x = "Iteration",
                fill = ""
            ) +
            theme(text = element_text(size = input$textSize))
            
    })
    
    output$iter_plot <- renderPlotly({
        ggplotly(iter_plot())
    })
    
    output$download_iter_plot <- downloadHandler(
        filename =  function() {
            paste("Iteration_plot", input$plot_format, sep=".")
        },
        # content is a function with argument file. content writes the plot to the device
        content = function(file) {
            multi_save_plot(file, plot =  iter_plot(),
                            type = input$plot_format,
                            width = input$widthGraph,
                            height = input$heightGraph,
                            units = input$unitsGraph,
                            dpi = input$dpiGraph)
        } 
    )

# Main Plots --------------------------------------------------------------


    vk_plot <- reactive({
        req(tidy_data())
        
        data <- filter_iter(tidy_data())
        
        
        data %>%
            group_by(.,H_C, O_C, Iter) %>%
            mutate(n = n()) %>%
            ggplot(aes(O_C, H_C, fill = as.factor(Iter), size = n)) +
            geom_point(shape = 21) +
            scale_fill_viridis_d(end = .8) +
            labs(
                y = "H/C",
                x = "O/C",
                fill = "Iteration"
            ) +
            {if(input$set_limits) expand_limits(x = input$set_limits_x, y = input$set_limits_y)} +
            theme(text = element_text(size = input$textSize))
    })
    
    output$vk_plot <- renderPlotly({
        ggplotly(vk_plot(), source = "plotly_vk")
        })
    
    output$download_vk_plot <- downloadHandler(
        filename =  function() {
            paste("vanKrevelen_plot", input$plot_format, sep=".")
        },
        # content is a function with argument file. content writes the plot to the device
        content = function(file) {
            multi_save_plot(file, plot =  vk_plot(),
                            type = input$plot_format,
                            width = input$widthGraph,
                            height = input$heightGraph,
                            units = input$unitsGraph,
                            dpi = input$dpiGraph)
        } 
    )
    
    # store all formula at a certain point
    
    #Display formula present at that point
    vk_datatable<- reactive({
        req(tidy_data())
        s <- event_data("plotly_click",source = "plotly_vk")
        if(is.null(s))return(NULL)
        filter_iter(tidy_data()) %>%
            dplyr::filter(O_C == s$x & H_C == s$y) 
    })
    
    output$vk_datatable <- DT::renderDataTable({
        req(vk_datatable())
        vk_datatable()
    })
        
   #OC plot
    OC_plot <- reactive({
        req(tidy_data())
        
        filter_iter(tidy_data()) %>%
            ggplot(aes(x = O/C, colour = as.factor(Iter))) +
            geom_density(size =2) +
            scale_colour_viridis_d(end = .8) +
            labs(
                y = "Density of molecules",
                x = "O/C",
                colour = "Iteration"
            )+
            {if(input$set_limits) expand_limits(x = input$set_limits_x)} +
            theme(text = element_text(size = input$textSize))
    })
    
    output$OC_plot <- renderPlot({
        OC_plot()
    })
    
    output$download_OC_plot <- downloadHandler(
        filename =  function() {
            paste("OC_plot", input$plot_format, sep=".")
        },
        # content is a function with argument file. content writes the plot to the device
        content = function(file) {
            multi_save_plot(file, plot =  OC_plot(),
                            type = input$plot_format,
                            width = input$widthGraph,
                            height = input$heightGraph,
                            units = input$unitsGraph,
                            dpi = input$dpiGraph)
        } 
    )
    
    #HC plot
    HC_plot <- reactive({
        req(tidy_data())
        
        filter_iter(tidy_data()) %>%
            ggplot(aes(x = H/C, colour = as.factor(Iter))) +
            geom_density(size =2) +
            scale_colour_viridis_d(end = .8) +
            labs(
                y = "Density of molecules",
                x = "H/C",
                colour = "Iteration"
            )+
            {if(input$set_limits) expand_limits(x = input$set_limits_y)}+
            theme(text = element_text(size = input$textSize))
    })
    
    output$HC_plot <- renderPlot({
        HC_plot()
    })
    
    output$download_HC_plot <- downloadHandler(
        filename =  function() {
            paste("HC_plot", input$plot_format, sep=".")
        },
        # content is a function with argument file. content writes the plot to the device
        content = function(file) {
            multi_save_plot(file, plot =  HC_plot(),
                            type = input$plot_format,
                            width = input$widthGraph,
                            height = input$heightGraph,
                            units = input$unitsGraph,
                            dpi = input$dpiGraph)
        } 
    )
    

# Class distribution ------------------------------------------------------

       class_dist <- reactive({
               filter_iter(tidy_data()) %>%
               mutate(class = case_when(
                   O == 0 ~ "HC",
                   O > 0 ~ paste0("O",O)
               )) %>%
               group_by(Iter, class) %>%
               summarise(n = n(), .groups = "drop") %>%
               complete(class,Iter, fill = list(n = 0)) %>%
               ggplot(aes(
                   x = factor(class, levels = str_sort(unique(class), numeric = TRUE), ordered = TRUE),
                   y = n,
                   fill = factor(Iter)
                   )) +
               geom_col(position = position_dodge2(preserve = "single")) +
               scale_fill_viridis_d(end = .8) +
               labs(
                   y = "Number of molecules",
                   x = "Heteroatom class",
                   fill = "Iteration"
               )+
               theme(text = element_text(size = input$textSize))
               
       }) 
    
    output$class_plot <- renderPlotly({
        ggplotly(class_dist())
    })
    
    output$download_class_plot <- downloadHandler(
        filename =  function() {
            paste("class_plot", input$plot_format, sep=".")
        },
        # content is a function with argument file. content writes the plot to the device
        content = function(file) {
            multi_save_plot(file, plot =  class_dist(),
                            type = input$plot_format,
                            width = input$widthGraph,
                            height = input$heightGraph,
                            units = input$unitsGraph,
                            dpi = input$dpiGraph)
        } 
    )
    

# Violin Plot -------------------------------------------------------------
    get_classes<-function(tidy_data){
        tidy_data %>%
            mutate(class = case_when(
                O == 0 ~ "HC",
                O > 0 ~ paste0("O",O)
            )) %>%
            pull(class) %>%
            unique()
    }
    
    observe({
        req(tidy_data())
        updateSelectInput(session = shiny::getDefaultReactiveDomain(), inputId = "violin_class_filter", label = "Select classes", choices = get_classes(tidy_data()))
    })
    
    density_data <- reactive({
        filter_iter(tidy_data()) %>%
            mutate(class = case_when(
                O == 0 ~ "HC",
                O > 0 ~ paste0("O",O)
            )) %>%
            {if(input$violin_filter_classes) filter(., class %in% input$violin_class_filter) else .} %>%
            ggplot(aes(
                x = factor(class, levels = str_sort(unique(class), numeric = TRUE), ordered = TRUE),
                y = H_C,
                fill = factor(Iter)
                )) +
            geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
            scale_fill_viridis_d(end = .8) +
            labs(
                y = "H/C",
                x = "Heteroatom class",
                fill = "Iteration"
            ) +
            {if(input$set_limits) expand_limits(y = input$set_limits_y)}+
            theme(text = element_text(size = input$textSize))
        
    })
    
    output$violin_plot <- renderPlot({
        density_data()
    })
    
    output$download_violin_plot <- downloadHandler(
        filename =  function() {
            paste("violin_plot", input$plot_format, sep=".")
        },
        # content is a function with argument file. content writes the plot to the device
        content = function(file) {
            multi_save_plot(file, plot =  density_data(),
                            type = input$plot_format,
                            width = input$widthGraph,
                            height = input$heightGraph,
                            units = input$unitsGraph,
                            dpi = input$dpiGraph)
        } 
    )
    

# H/C material ------------------------------------------------------------

  HC_data <-reactive({
    req(tidy_data())
    filter_iter(tidy_data()) %>%
      filter(O == 0) %>%
      mutate(HC_val = case_when(
        H_C <= 0.67 ~ "Aromatic",
        H_C <= 1.67 ~ "Napthenic",
        H_C < 2 ~ "Other H/C < 2",
        H_C >= 2 ~ "Paraffinic",
        TRUE ~ "Other"
      )) %>%
      group_by(Iter, HC_val) %>%
      summarise(n = n(), .groups = "drop")
  })  
    
    output$HC_datatable <- DT::renderDataTable({
      req(HC_data())
      HC_data()
    })  
    
    HC_val_plot <- reactive({
      req(HC_data())
      HC_data() %>%
        ggplot(aes(Iter,n, colour = HC_val, linetype = HC_val))+
        geom_line(size = input$HC_val_line_size)+
        geom_point(size = input$HC_val_point_size)+
        labs(
          y = "Number of molecules",
          x = "Iteration",
          colour = "Hydrocarbon type",
          linetype = "Hydrocarbon type"
        )+
        scale_x_continuous(n.breaks = 10) +
        theme(text = element_text(size = input$textSize))
    })
    
    output$HC_val_plot <- renderPlot({
      HC_val_plot()
    })
    
    output$download_HC_val_plot <- downloadHandler(
      filename =  function() {
        paste("HC_val_plot", input$plot_format, sep=".")
      },
      # content is a function with argument file. content writes the plot to the device
      content = function(file) {
        multi_save_plot(file, plot =  HC_val_plot(),
                        type = input$plot_format,
                        width = input$widthGraph,
                        height = input$heightGraph,
                        units = input$unitsGraph,
                        dpi = input$dpiGraph)
      } 
    )
    
    output$download_HC_val_data <- downloadHandler(
      filename = function() {
        paste("HC_data", ".csv", sep = "")
      },
      content = function(file) {
        write_csv(HC_data(), file)
      }
    )

# Generate gif ------------------------------------------------------------
    output$vk_gif <- renderImage({
        # A temp file to save the output.
        # This file will be removed later by renderImage
        outfile <- tempfile(fileext='.gif')
        
        # now make the animation
        p = vk_plot() +
            transition_manual(Iter) +
            labs(subtitle = "Iteration: {current_frame}")
        
        anim_save("outfile.gif", p) # New
        
        # Return a list containing the filename
        list(src = "outfile.gif",
             contentType = 'image/gif'
             # width = 400,
             # height = 300,
             # alt = "This is alternate text"
        )}, deleteFile = TRUE)
    
    

# Comparison Section ------------------------------------------------------

    stored_data<-reactive({
        req(!is.null(reactive_values$sim_data))
      # 
      # reaction <- reactive_values$reaction_data %>%
      #   map_chr(~stringr::str_c(paste0(.x$reaction,"(",.x$rep,")")))

        tibble(ID = reactive_values$sim_names)
    })

    output$stored_data <- renderTable({
        req(stored_data())
      stored_data()
    })
    
    observe({
      req(!is.null(reactive_values$sim_data))
      updatePickerInput(session = shiny::getDefaultReactiveDomain(),
                        inputId = "selected_entries",
                        "Select data to use",
                        choices = reactive_values$sim_names
      )
      updateSelectInput(session = shiny::getDefaultReactiveDomain(),
                        inputId = "combined_violin_class_filter",
                        label = "Select classes",
                        choices = get_classes(map_df(reactive_values$sim_data,~filter(.x, Iter == max(Iter))))
      )
    })
    
    
    combined_data<-reactive({
      req(length(reactive_values$sim_data) > 0)
      req(length(reactive_values$sim_names) > 0)
      
      
      
      data <- reactive_values$sim_data %>%
        set_names(.,reactive_values$sim_names) %>%
        map_df(~filter(.x, Iter == max(Iter)), .id = "id") %>%
        filter(id %in% input$selected_entries)
      
      if(nrow(data) >0) return(data) else return(NULL)
    })
    

## Combined VK -------------------------------------------------------------

    
    combined_vk_plot <- reactive({
      req(combined_data())
      
      combined_data() %>%
        group_by(.,H_C, O_C, Iter) %>%
        mutate(n = n()) %>%
        ggplot(aes(O_C, H_C, fill = id, size = n)) +
        geom_point(shape = 21) +
        scale_fill_viridis_d(end = .8) +
        labs(
          y = "H/C",
          x = "O/C",
          fill = "Sample"
        ) +
        {if(input$set_limits) expand_limits(x = input$set_limits_x, y = input$set_limits_y)} +
        theme(text = element_text(size = input$textSize))
        
    })
    
    output$combined_vk_plot <- renderPlot({
      combined_vk_plot()
    })
    
    output$download_combined_vk_plot <- downloadHandler(
      filename =  function() {
        paste("combined_vk_plot", input$plot_format, sep=".")
      },
      # content is a function with argument file. content writes the plot to the device
      content = function(file) {
        multi_save_plot(file, plot =  combined_vk_plot(),
                        type = input$plot_format,
                        width = input$widthGraph,
                        height = input$heightGraph,
                        units = input$unitsGraph,
                        dpi = input$dpiGraph)
      } 
    )
    

## Combined class plot -----------------------------------------------------

    combined_class_plot<- reactive({
      req(combined_data())
      
      combined_data() %>%
        mutate(class = case_when(
          O == 0 ~ "HC",
          O > 0 ~ paste0("O",O)
        )) %>%
        group_by(id, class) %>%
        summarise(n = n(), .groups = "drop") %>%
        complete(class,id, fill = list(n = 0)) %>%
        ggplot(aes(
          x = factor(class, levels = str_sort(unique(class), numeric = TRUE), ordered = TRUE),
          y = n,
          fill = id
        )) +
        geom_col(position = position_dodge2(preserve = "single")) +
        scale_fill_viridis_d(end = .8) +
        labs(
          y = "Number of molecules",
          x = "Heteroatom class",
          fill = "Sample"
        )+
        theme(text = element_text(size = input$textSize))
      
    }) 
    
    output$combined_class_plot <- renderPlotly({
      ggplotly(combined_class_plot())
    })
    
    output$download_combined_class_plot <- downloadHandler(
      filename =  function() {
        paste("combined_class_plot", input$plot_format, sep=".")
      },
      # content is a function with argument file. content writes the plot to the device
      content = function(file) {
        multi_save_plot(file, plot =  combined_class_plot(),
                        type = input$plot_format,
                        width = input$widthGraph,
                        height = input$heightGraph,
                        units = input$unitsGraph,
                        dpi = input$dpiGraph)
      } 
    )

## Combined violin ---------------------------------------------------------

    
    
    
    
    combined_violin_plot <- reactive({
      req(combined_data())
      
      combined_data() %>%
        mutate(class = case_when(
          O == 0 ~ "HC",
          O > 0 ~ paste0("O",O)
        )) %>%
        {if(input$combined_violin_filter_classes) filter(., class %in% input$combined_violin_class_filter) else .} %>%
        ggplot(aes(
          x = factor(class, levels = str_sort(unique(class), numeric = TRUE), ordered = TRUE),
          y = H_C,
          fill = id
        )) +
        geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
        scale_fill_viridis_d(end = .8) +
        labs(
          y = "H/C",
          x = "Heteroatom class",
          fill = "Sample"
        ) +
        {if(input$set_limits) expand_limits(y = input$set_limits_y)}+
        theme(text = element_text(size = input$textSize))
      
    })
    
    output$combined_violin_plot <- renderPlot({
      combined_violin_plot()
    })
    
    output$download_combined_violin_plot <- downloadHandler(
      filename =  function() {
        paste("combined_violin_plot", input$plot_format, sep=".")
      },
      # content is a function with argument file. content writes the plot to the device
      content = function(file) {
        multi_save_plot(file, plot =  combined_violin_plot(),
                        type = input$plot_format,
                        width = input$widthGraph,
                        height = input$heightGraph,
                        units = input$unitsGraph,
                        dpi = input$dpiGraph)
      } 
    )
    
    ## Combined H/C material ------------------------------------------------------------
    
    combined_HC_data <-reactive({
      req(combined_data())
      
      combined_data() %>%
        filter(O == 0) %>%
        mutate(HC_val = case_when(
          H_C <= 0.67 ~ "Aromatic",
          H_C <= 1.57 ~ "Napthenic",
          H_C < 2 ~ "Other H/C < 2",
          H_C >= 2 ~ "Paraffinic",
          TRUE ~ "Other"
        )) %>%
        group_by(id, HC_val) %>%
        summarise(n = n(), .groups = "drop")
    })  
    
    output$combined_HC_datatable <- DT::renderDataTable({
      req(combined_HC_data())
      combined_HC_data()
    })  
    
    combined_HC_val_plot <- reactive({
      req(combined_HC_data())
      combined_HC_data() %>%
        ggplot(aes(id, n, colour = HC_val))+
        geom_point(size = input$HC_val_point_size)+
        labs(
          y = "Number of molecules",
          x = "sample",
          colour = "Hydrocarbon type"
        )+
        theme(text = element_text(size = input$textSize))
    })
    
    output$combined_HC_val_plot <- renderPlot({
      combined_HC_val_plot()
    })
    
    output$download_combined_HC_val_plot <- downloadHandler(
      filename =  function() {
        paste("combined_HC_val_plot", input$plot_format, sep=".")
      },
      # content is a function with argument file. content writes the plot to the device
      content = function(file) {
        multi_save_plot(file, plot =  combined_HC_val_plot(),
                        type = input$plot_format,
                        width = input$widthGraph,
                        height = input$heightGraph,
                        units = input$unitsGraph,
                        dpi = input$dpiGraph)
      } 
    )
        
# End ---------------------------------------------------------------------
    
    
}  

# Create APP --------------------------------------------------------------

shinyApp(ui, server)  