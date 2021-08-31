

#helper functions

string_to_list<-function(string_vec){
    if(any(is.na(as.numeric(string_vec[1:3])))) return(NULL)
    
    rep_value <- if(is.na(as.numeric(string_vec[4]))) 1 else as.numeric(string_vec[4])
    
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


atom_check<-function(formula){
    str_extract_all(formula,'([A-Z][a-z]?)') %>%
        map_lgl(~all(.x %in% c("C","H","N","O","S")))
}

ion_polarity_finder <- function(data){
    output <- data %>%
        mutate(electron_status = ifelse (dbe != as.integer(dbe),'even','odd')) %>%
        filter(electron_status == 'even' & atom_check(formula)) %>%
        mutate(calc = C-H/2+N/2+1) %>%
        mutate(polarity = ifelse(calc-dbe==-0.5,'neg','pos')) %>%
        pull(polarity)
    
    if(length(unique(output))!=1){stop("Something went wrong")}
    return(unique(output))
}

adjust_H<-function(H, electron_status, ion_polarity){
    case_when(
        electron_status == "even" & ion_polarity == "neg" ~ H + 1,
        electron_status == "even" & ion_polarity == "pos" ~ H - 1,
        electron_status == "odd" ~ H + 0
    )
}

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
    radioButtons(inputId = inputId,label = label, choices = list("png", "pdf", "eps", "tiff","svg"))
}

# Packages to Load --------------------------------------------------------
suppressPackageStartupMessages(c(
    library(shiny),
    library(shinydashboard),
    library(shinyWidgets),
    library(magrittr),
    library(tidyverse),
    library(glue),
    library(plotly)
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
                menuItem("Options",tabName = "Options"),
                hr(),
                    materialSwitch("set_limits","Set plot limits"),
                    conditionalPanel(condition = 'input.set_limits',
                                 numericRangeInput("set_limits_x", "Set O/C axis limits", value = c(0,3)),
                                 numericRangeInput("set_limits_y", "Set H/C axis limits", value = c(0,3))
                )
                
                
                
                
    ))


# Body --------------------------------------------------------------------

body <- dashboardBody(
    tabItems(
        
        # Main Tab ----------------------------------------------------------------
        tabItem(tabName = "Input",
                fluidPage(
                    box(width = NULL,
                        fileInput(inputId = "file_input",
                                  label = "Choose data",
                                  multiple= TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv"))
                        ),
                    box(width = NULL,
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
                                actionButton("simulate","Simulate Reaction")
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
                    )
                )),
        tabItem(tabName = "Options",
                fluidPage(
                    box(width = 6,
                        h4("Options for downloaded graphics"),
                           radioButtons(inputId = "unitsGraph", label = "Select the units", choices = list("mm", "cm","in"), selected = "in"),
                           numericInput(inputId = "dpiGraph",label = "DPI of the graphics:", value = 600, max = 2000),
                           numericInput(inputId = "widthGraph",label = "Width of the graphics:", value = 12),
                           numericInput(inputId = "heightGraph",label = "Height of the graphics:", value = 8),
                           select_plot_file_type(inputId = "plot_format",label = "Select the file type")
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
   
    raw_masslist<-reactive({
        req(!is.null(input$file_input))
        validate(
            need(!is.na(KairosMSfunctions::determine_file_format(input$file_input$datapath)), "Ensure file is the correct type")
        )
        KairosMSfunctions::Master_ReadeR(file_to_read = input$file_input$datapath)
    })
    
    
    
    
    
    output$debug <- renderText({
        "debug"
        })
    
    
    

# Reaction details input --------------------------------------------------

    
    raw_vector_input <-reactive({
        req(!is.null(input$list_input))
        
        str_split(input$list_input,"\n") %>%
            unlist() %>% 
            str_extract_all("(-?)[:digit:]+") %>%
            keep(~length(.x) >= 3) %>%
            map(~string_to_list(.x)) %>%
            keep(~!is_null(.x))
        
    })
    
    reaction_details <- reactive({
        req(length(raw_vector_input()) > 0)
        map_df(raw_vector_input(), ~extract_inputs(.x)) %>%
            mutate(reaction = paste0("C",C,"H",H,"O",O))
    })
    
    output$reaction_vectors<-DT::renderDataTable({
        req(reaction_details())
        reaction_details()
        })
    
    

# Prepare data ------------------------------------------------------------

    
    filtered_data <- reactive({
        req(raw_masslist())
        ion_polarity <- ion_polarity_finder(raw_masslist())
        raw_masslist() %>%
            filter(isotope == FALSE) %>%
            filter(str_detect(class,"N|S|B|Na", negate = TRUE)) %>%
            mutate(electron_status = ifelse(dbe != as.integer(dbe),'even','odd'))
            #mutate(nH = adjust_H(H, electron_status, ion_polarity)) 
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
            )
            
    })
    
    output$iter_plot <- renderPlotly({
        ggplotly(iter_plot())
    })
    
    output$download_iter_plot <- downloadHandler(
        filename =  function() {
            paste("Iteration_plot", input$dl_VKSim_format, sep=".")
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
            {if(input$set_limits) expand_limits(x = input$set_limits_x, y = input$set_limits_y)}
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
            ggplot(aes(x = O/C, y = after_stat(count), colour = as.factor(Iter))) +
            geom_density(size =2) +
            scale_colour_viridis_d(end = .8) +
            labs(
                y = "Number of molecules",
                x = "O/C",
                colour = "Iteration"
            )+
            {if(input$set_limits) expand_limits(x = input$set_limits_x)}
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
            ggplot(aes(x = H/C, y = after_stat(count), colour = as.factor(Iter))) +
            geom_density(size =2) +
            scale_colour_viridis_d(end = .8) +
            labs(
                y = "Number of molecules",
                x = "H/C",
                colour = "Iteration"
            )+
            {if(input$set_limits) expand_limits(x = input$set_limits_y)}
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
               )
               
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
            {if(input$set_limits) expand_limits(y = input$set_limits_y)}
        
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
        
# End ---------------------------------------------------------------------
    
    
}  

# Create APP --------------------------------------------------------------

shinyApp(ui, server)  