# AVAILABLE METHODS
norm_options <- c(
    "Raw",
    "Illumina",
    "Funnorm",
    "Noob",
    "SWAN",
    "Quantile",
    "Noob+Quantile"
)

hclust_methods <- c(
    "single",
    "complete",
    "average",
    "mcquitty",
    "median",
    "centroid"
)

controlNames <- c(
    "BISULFITE CONVERSION I",
    "BISULFITE CONVERSION II",
    "HYBRIDIZATION",
    "SPECIFICITY I",
    "SPECIFICITY II",
    "TARGET REMOVAL",
    "BISULFITE CONVERSION I",
    "EXTENSION",
    "STAINING",
    "NON-POLYMORPHIC"
)

convertMenuItem <- function(mi,tabName) {
    mi$children[[1]]$attribs['data-toggle']="tab"
    mi$children[[1]]$attribs['data-value'] = tabName
    if(length(mi$attribs$class)>0 && mi$attribs$class=="treeview"){
        mi$attribs$class=NULL
    }
    mi
}


library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyWidgets)
library(shinycssloaders)
library(dplyr)

shinyUI(
    dashboardPage( 

    dashboardHeader(title = "VisualMethyl"),

    dashboardSidebar(
        
        sidebarMenu(id = "menu",
                    menuItem("Data", tabName = "data"),
                    convertMenuItem(menuItem("Analysis", tabName = "analysis", startExpanded = TRUE,
                    menuSubItem("QC", tabName = "qc"),
                    menuSubItem("Exploratory Analysis", tabName = "exploratory_analysis"),
                    menuSubItem("DMPs/DMRs", tabName = "dmp_dmr"),
                    menuSubItem("Functional Enrichment", tabName = "functional_enrichment"),
                    menuSubItem("Survival", tabName = "survival"),
                    menuSubItem("Predicted Models", tabName = "predicted_models"),
                    menuSubItem("External Sources", tabName = "external_sources"),
                    menuSubItem("Genome Browser", tabName = "genome_browser")), "analysis"),
                    menuItem("Export", tabName = "export"),
                    menuItem("Help", tabName = "help")
        )
    ),
    
    dashboardBody(includeCSS("www/style.css"), 
        
        tabItems(
            
            tabItem(
                tabName = "data",
                verticalLayout( 
                    box(title = "INPUT DATA", width = 12, collapsible = FALSE, status = "primary",
                            fileInput("input_data", "Upload input data (.zip)", multiple = FALSE, accept = ".zip"),
                            uiOutput("ui_input_data")
                    ),
                    fluidRow(style = "margin-left: 2px; margin-right: 2px",
                    conditionalPanel(
                        "input.b_input_data > 0",
                        box(title = "SELECTION OPTIONS", id = "selection_options",
                                collapsible = FALSE, status = "primary",
                            width = 3,
                            selectInput("select_input_samplenamevar", "", c()),
                            selectInput("select_input_groupingvar", "", c()),
                            selectInput("select_input_donorvar", "", c()),
                            pickerInput(
                                inputId = "selected_samples",
                                label = "",
                                choices = c(),
                                options = list(
                                    `actions-box` = TRUE,
                                    size = 10,
                                    `selected-text-format` = "count > 3"
                                ),
                                multiple = TRUE
                            ),
                            selectInput("select_input_sex", "", c()),
                            selectInput("select_input_age", "", c()),
                            actionButton("button_input_next", "Continue to Analysis")
                        )
                    ),
                    
                    # Box 2: Table
                    conditionalPanel(
                        "input.b_input_data > 0", box(title = "SAMPLES TABLE",
                                                          collapsible = FALSE, status = "primary",
                        width = 9,
                        withSpinner(DT::DTOutput("samples_table"))
                        
                    )))
                    
                )
            ),
            tabItem(
                tabName = "analysis",
                verticalLayout(
                splitLayout(style = "margin-left: 14px; margin-right: 14px;", cellArgs = list(style = "padding: 4px;"),
                            actionButton("b_qc", width = "100%", label = "QC", class = "btn-info", style = "padding:3.55vw;", disable = TRUE),
                            actionButton("b_exploratory_analysis", width = "100%",label = HTML("Exploratory <br/> Analysis"), class = "btn-info", style = "padding:1.75vw;"),
                            actionButton("b_dmp_dmr", width = "100%", label = "DMPs/DMRs", class = "btn-info", style = "padding:3.55vw;"),
                            actionButton("b_functional_enrichment", width = "90%",label = HTML("Functional <br/> Enrichment"), class = "btn-info", style = "padding:1.75vw;")
                ),
                splitLayout(style = "margin-left: 14px; margin-right: 14px", cellArgs = list(style = "padding: 4px;"),
                            actionButton("b_survival", width = "100%", label = "Survival", class = "btn-info", style = "padding:3.55vw;"),
                            actionButton("b_predicted_models", width = "100%", label = HTML("Predicted <br/> Models"), class = "btn-info", style = "padding:1.75vw;"),
                            actionButton("b_external_sources", width = "100%", label = HTML("External <br/> Sources"), class = "btn-info", style = "padding:1.75vw;"),
                            actionButton("b_genome_browser", width = "90%", label = HTML("Genome <br/> Browser"), class = "btn-info", style = "padding:1.75vw;")
                ))
            ),
            tabItem(
                tabName = "qc",
                fluidPage(
                    # Box1
                    sidebarPanel(
                        width = 3,
                        selectInput("select_minfi_norm", "Select Normalization", norm_options),
                        div(
                            margin_left = "50px",
                            switchInput(
                                inputId = "select_minfi_dropcphs",
                                label = "Drop CpHs",
                                labelWidth = "fit",
                                value = TRUE,
                                inline = TRUE
                            ),
                            
                            switchInput(
                                inputId = "select_minfi_dropsnps",
                                label = "Drop SNPs",
                                labelWidth = "fit",
                                value = TRUE,
                                inline = TRUE
                            )
                        ),
                        
                        conditionalPanel(
                            "input.select_minfi_dropsnps",
                            sliderInput(
                                inputId = "slider_minfi_maf",
                                label = "Minimum MAF to filter",
                                min = 0,
                                max = 1,
                                step = 0.01,
                                value = 0,
                                width = "75%"
                            )
                        ),
                        
                        switchInput(
                            inputId = "select_minfi_chromosomes",
                            label = "Drop X/Y Chr.",
                            labelWidth = "fit",
                            value = FALSE
                        ),
                        
                        shinyjs::disabled(actionButton("button_minfi_select", "Select")),
                        h4(),
                        textOutput("text_minfi_probes"),
                        conditionalPanel( "input.button_minfi_select > 0",
                            #"typeof output.green_intensities_plot != 'undefined'",
                        br(),
                        downloadButton("download_html", label = "Download HTML"),
                        br(),
                        downloadButton("download_pdf", label = "Download PDF")
                        )
                    ),
                    
                    mainPanel(
                        width = 9,
                        verticalLayout(
                        box(title = "INTENSITIES BOXPLOTS", width = 12, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                h4("Green channel intensities"), 
                                withSpinner(plotOutput("green_intensities_plot")),
                                h4("Red channel intensities"), 
                                withSpinner(plotOutput("red_intensities_plot"))
                        ),
                        box(title = "FAILED PROBES", width = 12, collapsible = TRUE, collapsed = TRUE, status = "primary",
                            column(width = 9,
                                h4("Failure Rate Plot"), 
                                withSpinner(plotly::plotlyOutput("failure_rate_plot"))),
                            column(width = 3,
                                h4("Failure Rate Table"), 
                                withSpinner(DT::DTOutput("failure_rate_table")))
                        ),
                        box(title = "CONTROL TYPES", width = 12, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                selectInput("controlType", "Choose a control type:",
                                            choices = controlNames, selected = 1),
                                selectInput("select_slide", "Select slide:", choices = c()),
                                withSpinner(plotOutput("controlTypePlotGreen")),
                                withSpinner(plotOutput("controlTypePlotRed"))
                        ),
                        box(title = "PREDICTED SEX", width = 12, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                h4("X vs Y chromosomes signal plot"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_sex")),
                                withSpinner(DT::DTOutput("table_minfi_sex"))
                        ),
                        box(title = "DENSITY PLOT", width = 12, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                selectInput("probeType", "Choose a probe type for the density curves:",
                                            choices = c("I-Green","I-Red","II"),
                                            selected="I-Green"),
                                h4("Raw"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_densityplotraw")),
                                h4("Normalized"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_densityplot"))
                        ),
                        box(title = "SNP ANALYSIS", width = 12, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                h4("SNPs beta-values (Raw)"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_snps"))
                        ),
                        box(title = "BATCH EFFECTS", width = 12, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                h4("Normalized"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_corrplot")),
                                selectInput("select_minfi_typecorrplot", "Select data to plot", choices = c("p.value", "correlation value"), selected = "p.value"),
                                withSpinner(DT::DTOutput("table_minfi_corrplot"))
                        ))
                    )
                )
            ),
            tabItem(
                tabName = "exploratory_analysis",
                fluidPage(
                        verticalLayout(
                            box(title = "VIOLIN PLOT", width = 12, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                    h4("Raw"),
                                    withSpinner(plotOutput("graph_violin_raw")),
                                    h4("Normalized"),
                                    withSpinner(plotOutput("graph_violin_normalized"))
                            ),
                            box(title = "PRINCIPAL COMPONENT ANALYSIS", width = 12, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                    withSpinner(plotly::plotlyOutput("graph_minfi_pcaplot")),
                                    withSpinner(DT::DTOutput("table_minfi_pcaplot")),
                                    column(
                                        6,
                                        selectInput(
                                            inputId = "select_minfi_pcaplot_pcx",
                                            choices = c(),
                                            label = "Select x variable"
                                        ),
                                        
                                        selectInput(
                                            inputId = "select_minfi_pcaplot_color",
                                            choices = c(),
                                            label = "Select color variable"
                                        )
                                    ),
                                    column(
                                        6,
                                        selectInput(
                                            inputId = "select_minfi_pcaplot_pcy",
                                            choices = c(),
                                            label = "Select y variable"
                                        )
                                    ),
                                    actionButton("button_pca_update", "Update")
                            ),
                            box(title = "HEATMAP", width = 12, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                    h4("Random 1000 CpGs"),
                                withSpinner(plotOutput("graph_random_heatmap")),
                                h4("Top 1000 variable CpGs"),
                                withSpinner(plotOutput("graph_top_heatmap"))
                                
                            ),
                            box(title = "DECONVOLUTION", width = 12, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                    withSpinner(plotOutput("deconvolution_heatmap"))
                            ),
                            box(title = "AGE METH", width = 12, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                    withSpinner(DT::DTOutput("table_age"))
                            ),
                            #box(title = "GENOME ANALYSIS PLOT", collapsible = TRUE, collapsed = TRUE, status = "primary",
                            #        h1("genome analysis plot")
                            #),
                            box(title = "HYPO/HYPER", width = 12, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                sliderInput(
                                    inputId = "slider_beta",
                                    label = "Beta threshold",
                                    min = 0,
                                    max = 1,
                                    step = 0.01,
                                    value = 0.33,
                                    width = "75%"
                                ),
                                pickerInput(
                                    inputId = "selected_samples_h",
                                    label = "",
                                    choices = c(),
                                    options = list(
                                        `actions-box` = TRUE,
                                        size = 10,
                                        `selected-text-format` = "count > 3"
                                    ),
                                    multiple = TRUE
                                ),
                                actionButton("button_hyper_hypo_update", "Update"),
                                br(),
                                withSpinner(plotOutput("plot_chr")),
                                withSpinner(plotOutput("plot_relation_to_island")),
                                withSpinner(plotOutput("plot_group"))
                            #),
                            #box(title = "CIRCOS", collapsible = TRUE, collapsed = TRUE, status = "primary",
                            #        h1("circos")
                            )
                        )
                    
                )
            ),
            tabItem(
                tabName = "dmp_dmr",
                fluidPage(
                    tabsetPanel(
                        tabPanel("DMPs", style = "margin-top: 10px;", sidebarLayout(
                            sidebarPanel(
                                width = 3,
                                h4("Linear Model Options"),
                                
                                pickerInput(
                                    inputId = "select_limma_voi",
                                    label = "Select Variable of Interest",
                                    choices = c(),
                                    multiple = FALSE
                                ),
                                
                                pickerInput(
                                    inputId = "checkbox_limma_covariables",
                                    label = "Select linear model covariables",
                                    choices = c(),
                                    multiple = TRUE,
                                    options = list(
                                        `actions-box` = TRUE,
                                        size = 10,
                                        `selected-text-format` = "count > 3"
                                    )
                                ),
                                
                                pickerInput(
                                    inputId = "checkbox_limma_interactions",
                                    label = "Select linear model interactions",
                                    choices = c(),
                                    multiple = TRUE,
                                    options = list(
                                        `actions-box` = TRUE,
                                        size = 10,
                                        `selected-text-format` = "count > 3"
                                    )
                                ),
                                
                                switchInput(
                                    inputId = "select_limma_weights",
                                    label = "Array Weights",
                                    labelWidth = "80px",
                                    value = FALSE
                                ),
                                
                                
                                shinyjs::disabled(
                                    actionButton("button_limma_calculatemodel", "Generate Model")
                                ),
                                tags$br(),
                                uiOutput("button_limma_calculatedifs_container")
                            ),
                            mainPanel(
                                width = 9,
                                box(title = "DMP TABLE AND OPTIONS", width = 12, closable = FALSE, collapsible = FALSE, status = "primary",
                                    h4("DMP counts in each contrast"),
                                    tableOutput("table_limma_difcpgs") %>% shinycssloaders::withSpinner(),
                                    fluidRow(
                                        column(6,
                                            h4("Group options"),
                                            selectizeInput(
                                                "select_limma_groups2plot",
                                                "Groups to plot",
                                                c(),
                                                multiple = TRUE,
                                                options = list(plugins = list("remove_button", "drag_drop"))
                                            ),
                                            selectizeInput(
                                                "select_limma_contrasts2plot",
                                                "Contrasts to plot",
                                                c(),
                                                multiple = TRUE,
                                                options = list(plugins = list("remove_button", "drag_drop"))
                                            ),
                                            
                                            h4("Data options"),
                                            switchInput(
                                                inputId = "select_limma_removebatch",
                                                label = "Remove Batch Effect",
                                                labelWidth = "100px",
                                                value = FALSE,
                                                disabled = TRUE
                                            )
                                        ),
                                        
                                        column(6,
                                            h4("Filtering options"),
                                            sliderInput("slider_limma_deltab", "Min. DeltaBeta", 0, 1, 0.2),
                                            sliderInput("slider_limma_adjpvalue", "Max. FDR", 0, 1, 0.05),
                                            sliderInput("slider_limma_pvalue", "Max. p-value", 0, 1, 1)
                                        ),
                                        actionButton("button_limma_tablecalc", "Update")
                                    )
                                        
                                ),
                                box(title = "DMP HEATMAP", width = 12, closable = TRUE, collapsible = TRUE, collapsed = TRUE,  status = "primary",
                                    textOutput("text_limma_heatmapcount"),
                                    uiOutput("graph_limma_heatmapcontainer"),
                                    
                                    h4("Clustering options",
                                       align =
                                           "left"
                                    ),
                                    
                                    fluidRow(
                                        column(5,
                                            selectInput(
                                                "select_limma_clusteralg",
                                                "Clustering algorithm",
                                                c(
                                                    "single",
                                                    "complete",
                                                    "average",
                                                    "mcquitty",
                                                    "median",
                                                    "centroid"
                                                ),
                                                "average"
                                            ),
                                            
                                            selectInput(
                                                "select_limma_clusterdist",
                                                "Distance Function",
                                                c("pearson", "spearman", "kendall", "euclidean"),
                                                "pearson"
                                            ),
                                            
                                            selectInput("select_limma_scale", "Scale", c("row", "none"), "row"),
                                            tags$br()
                                        ),
                                        
                                        column(
                                            3,
                                            offset = 1,
                                            tags$br(),
                                            
                                            switchInput(
                                                inputId = "select_limma_graphstatic",
                                                label = "Static Graph",
                                                labelWidth = "100px",
                                                value = TRUE
                                            ),
                                            
                                            switchInput(
                                                inputId = "select_limma_colv",
                                                label = "Column Dendro.",
                                                labelWidth = "100px",
                                                value = TRUE
                                            ),
                                            
                                            switchInput(
                                                inputId = "select_limma_colsidecolors",
                                                label = "Column Colors",
                                                labelWidth = "100px",
                                                value = FALSE
                                            )
                                        ),
                                        
                                        column(
                                            3,
                                            
                                            tags$br(),
                                            
                                            switchInput(
                                                inputId = "select_limma_rowsidecolors",
                                                label = "Row Colors",
                                                labelWidth = "100px",
                                                value = FALSE
                                            ),
                                            
                                            conditionalPanel(
                                                "input.select_limma_rowsidecolors",
                                                numericInput(
                                                    "select_limma_knumber",
                                                    "Clusters number",
                                                    value = 2,
                                                    min = 1,
                                                    max = Inf,
                                                    step = 1
                                                )
                                            ),
                                            
                                            
                                            shinyjs::disabled(actionButton("button_limma_heatmapcalc", "Update"))
                                        )
                                    )   
                                ),
                                box(title = "DMPs ANNOTATION", width = 12, closable = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                    h4("DMP Boxplot"),
                                    plotOutput("graph_limma_indboxplot") %>% shinycssloaders::withSpinner(),
                                    h4("DMPs Annotation"),
                                    br(),
                                    DT::DTOutput("table_limma_ann") %>% shinycssloaders::withSpinner(),
                                    selectInput(inputId = "select_limma_anncontrast", label = "", choices = "", selected = ""),
                                    actionButton(inputId = "button_limma_indboxplotcalc", label = "Plot")
                                ),
                                box(title = "DMP MANHATTAN", width = 12, closable = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                    selectInput(inputId = "select_anncontrast_manhattan", label = "", choices = "", selected = ""),
                                    withSpinner(plotOutput("manhattan_plot"))
                                ),
                                box(title = "DMP VOLCANO", width = 12, closable = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                    selectInput(inputId = "select_anncontrast_volcano", label = "", choices = "", selected = ""),
                                    withSpinner(plotOutput("volcano_plot")) 
                                )
                            )
                        )),
                        tabPanel("DMRs", style = "margin-top: 10px;", sidebarLayout(
                            sidebarPanel(
                                width = 3,
                                
                                pickerInput(
                                    inputId = "select_dmrs_contrasts",
                                    label = "Contrasts to calculate",
                                    choices = c(),
                                    options = list(
                                        `actions-box` = TRUE,
                                        size = 10,
                                        `selected-text-format` = "count > 3"
                                    ),
                                    multiple = TRUE
                                ),
                                
                                pickerInput(
                                    inputId = "select_dmrs_regions",
                                    label = "Type of DMRs",
                                    choices = c("promoters", "genes", "CGI"),
                                    selected = c("promoters", "genes", "CGI"),
                                    options = list(
                                        `actions-box` = TRUE,
                                        size = 10,
                                        `selected-text-format` = "count > 3"
                                    ),
                                    multiple = TRUE
                                ),
                                
                                pickerInput(
                                    inputId = "select_dmrs_platform",
                                    label = "Array platform",
                                    choices = c("450k", "EPIC"),
                                    selected = c("EPIC"),
                                    multiple = FALSE
                                ),
                                
                                sliderInput(
                                    "slider_dmrs_cpgs",
                                    label = "Min. CpGs in DMR",
                                    min = 2,
                                    max = 50,
                                    value = 5
                                ),
                                
                                sliderInput(
                                    "slider_dmrs_permutations",
                                    label = "Number of permutations",
                                    min = 1000,
                                    max = 100000,
                                    value = 50000
                                ),
                                
                                shinyjs::disabled(actionButton("button_dmrs_calculate", "Calculate"))
                            ),
                            mainPanel(
                                width = 9,
                                box(title = "DMR TABLE AND OPTIONS", width = 12, closable = FALSE, collapsible = FALSE,
                                    h4("DMRs counts in each contrast"),
                                    tableOutput("table_dmrs_count") %>% shinycssloaders::withSpinner(),
                                    fluidRow(
                                        column(
                                            6,
                                            h4("Group options"),
                                            
                                            selectizeInput(
                                                "select_dmrs_groups2plot",
                                                "Groups to plot",
                                                c(),
                                                multiple = TRUE,
                                                options = list(plugins = list("remove_button", "drag_drop"))
                                            ),
                                            
                                            selectizeInput(
                                                "select_dmrs_contrasts2plot",
                                                "Contrasts to plot",
                                                c(),
                                                multiple = TRUE,
                                                options = list(plugins = list("remove_button", "drag_drop"))
                                            ),
                                            
                                            selectizeInput(
                                                "select_dmrs_regions2plot",
                                                "Regions to plot",
                                                c(),
                                                multiple = TRUE,
                                                options = list(plugins = list("remove_button", "drag_drop"))
                                            ),
                                            
                                            h4("Data options"),
                                            
                                            switchInput(
                                                inputId = "select_dmrs_removebatch",
                                                label = "Remove Batch Effect",
                                                labelWidth = "100px",
                                                value = FALSE,
                                                disabled = TRUE
                                            )
                                        ),
                                        
                                        column(
                                            6,
                                            h4("Filtering options"),
                                            sliderInput("slider_dmrs_deltab", "Min. DeltaBeta", 0, 1, 0),
                                            sliderInput("slider_dmrs_adjpvalue", "Max. FDR", 0, 1, 0.05),
                                            sliderInput("slider_dmrs_pvalue", "Max. p-value", 0, 1, 1)
                                        ),
                                        actionButton("button_dmrs_tablecalc", "Update")
                                    )
                                    
                                ),
                                box(title = "DMR HEATMAP", width = 12, closable = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                    h4("DMRs Heatmap"),
                                    textOutput("text_dmrs_heatmapcount"),
                                    uiOutput("graph_dmrs_heatmapcontainer"),
                                    
                                    h4("Clustering options",
                                       align =
                                           "left"
                                    ),
                                    
                                    fluidRow(
                                        column(
                                            5,
                                            selectInput(
                                                "select_dmrs_clusteralg",
                                                "Clustering algorithm",
                                                c(
                                                    "single",
                                                    "complete",
                                                    "average",
                                                    "mcquitty",
                                                    "median",
                                                    "centroid"
                                                ),
                                                "average"
                                            ),
                                            
                                            selectInput(
                                                "select_dmrs_clusterdist",
                                                "Distance Function",
                                                c("pearson", "spearman", "kendall", "euclidean"),
                                                "pearson"
                                            ),
                                            
                                            selectInput("select_dmrs_scale", "Scale", c("row", "none"), "row")#,
                                            #tags$br()
                                        ),
                                        
                                        column(
                                            3,
                                            offset = 1,
                                            #tags$br(),
                                            
                                            switchInput(
                                                inputId = "select_dmrs_graphstatic",
                                                label = "Static Graph",
                                                labelWidth = "100px",
                                                value = TRUE
                                            ),
                                            
                                            switchInput(
                                                inputId = "select_dmrs_colv",
                                                label = "Column Dendro.",
                                                labelWidth = "100px",
                                                value = TRUE
                                            ),
                                            
                                            switchInput(
                                                inputId = "select_dmrs_colsidecolors",
                                                label = "Column Colors",
                                                labelWidth = "100px",
                                                value = FALSE
                                            )
                                        ),
                                        
                                        column(
                                            3,
                                            
                                            #tags$br(),
                                            
                                            switchInput(
                                                inputId = "select_dmrs_rowsidecolors",
                                                label = "Row Colors",
                                                labelWidth = "100px",
                                                value = FALSE
                                            ),
                                            
                                            conditionalPanel(
                                                "input.select_dmrs_rowsidecolors",
                                                numericInput(
                                                    "select_dmrs_knumber",
                                                    "Clusters number",
                                                    value = 2,
                                                    min = 1,
                                                    max = Inf,
                                                    step = 1
                                                )
                                            ),
                                            
                                            shinyjs::disabled(actionButton("button_dmrs_heatmapcalc", "Update"))
                                        )
                                    )
                                ),
                                box(title = "DMR ANNOTATION", width = 12, closable = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                    h4("Genomic graph"),
                                    plotOutput("graph_dmrs_singledmr") %>% shinycssloaders::withSpinner(),
                                    # h4("GSEA graph"),
                                    # plotOutput("graph_dmrs_singlegsea") %>% shinycssloaders::withSpinner(),
                                    h4("DMRs table"),
                                    
                                    
                                    div(
                                        style = "display:inline-block",
                                        
                                        selectInput(
                                            "select_dmrs_selcont",
                                            label = "Contrast",
                                            choices = c()
                                        )
                                    ),
                                    
                                    div(
                                        style = "display:inline-block",
                                        selectInput("select_dmrs_selreg", label = "Region", choices = c())
                                    ),
                                    
                                    
                                    DT::DTOutput("table_dmrs_table") %>% shinycssloaders::withSpinner(),
                                    
                                    br(),
                                    
                                    actionButton("button_dmrs_graphsingle", "Plot")    
                                )
                            )
                        )
                                 
                        )
                    )
                )
            ),
            tabItem(
                tabName = "functional_enrichment",
                fluidPage(
                    box(title = "KEGG", width = 12, closable = FALSE, collapsible = TRUE, collapsed = TRUE, status = "primary",
                        withSpinner(plotOutput("plot_kegg"))
                    ),
                    box(title = "Gene Ontology (GO)", width = 12, closable = FALSE, collapsible = TRUE, collapsed = TRUE, status = "primary",
                        withSpinner(plotOutput("plot_go_mf")),
                        withSpinner(plotOutput("plot_go_bp")),
                        withSpinner(plotOutput("plot_go_cc"))
                    ),
                    box(title = "REACTOME", width = 12, closable = FALSE, collapsible = TRUE, collapsed = TRUE, status = "primary",
                        withSpinner(plotOutput("plot_reactome"))
                    )
                )
            ),
            tabItem(
                tabName = "survival",
                fluidPage(h1("Survival"))
            ),
            tabItem(
                tabName = "predicted_models",
                fluidPage(h1("Predicted Models"))
            ),
            tabItem(
                tabName = "external_sources",
                fluidPage(h1("External Sources"))
            ),
            tabItem(
                tabName = "genome_browser",
                fluidPage(h1("Genome Browser"))
            ),
            tabItem(
                tabName = "export",
                fluidPage(
                    shinyjs::useShinyjs(),
                    h3("Download RObjects"),
                    pickerInput(
                        inputId = "select_export_objects2download",
                        label = "Selected objects",
                        choices = c("RGSet", "GenomicRatioSet", "fit", "design", "ebayestables", "Bvalues", "Mvalues", "global_difs", "dmr_results"),
                        selected = c("RGSet", "GenomicRatioSet", "fit", "design", "ebayestables", "Bvalues", "Mvalues", "global_difs"),
                        options = list(
                            `actions-box` = TRUE,
                            size = 10,
                            `selected-text-format` = "count > 3"
                        ),
                        multiple = TRUE
                    ),
                    downloadButton("download_export_robjects"),
                    p(
                        "Press to download the R objects used for the analysis (RGSet, GenomicRatioSet, Bvalues, Mvalues, etc.)"
                    ),
                    h3("Download filtered bed files"),
                    
                    fluidPage(
                        div(
                            style = "display:inline-block",
                            selectInput(
                                "select_export_analysistype",
                                "Analysis type",
                                c("DMPs", "DMRs"),
                                selected = "by contrast"
                            )
                        ),
                        div(
                            style = "display:inline-block",
                            selectInput(
                                "select_export_bedtype",
                                "Subsetting mode",
                                c("by contrasts", "by heatmap cluster"),
                                selected = "by contrasts"
                            )
                        ),
                        div(
                            style = "display:inline-block",
                            selectInput(
                                "select_export_genometype",
                                "Genome version",
                                c("hg19", "hg38"),
                                selected = "hg19"
                            )
                        )
                    ),
                    
                    
                    downloadButton("download_export_filteredbeds"),
                    p(
                        "Press to download the created filtered lists of contrasts, or heatmap clusters,
        with the chosen criteria, in BED format."
                    ),
                    h3("Download Workflow Report"),
                    downloadButton("download_export_markdown"),
                    p(
                        "Press to download the report of all the steps follow and selected in the pipeline, and the results."
                    ),
                    h3("Download Custom R Script"),
                    downloadButton("download_export_script"),
                    p(
                        "Press to download an R script with the main parameters and steps follow in the DMP/DMR pipeline, to reproduce the results later outside the shiny application."
                    ),
                    h3("Download Heatmap"),
                    selectInput(
                        "select_export_heatmaptype",
                        label = "Heatmap type",
                        choices = c("DMPs", "DMRs"),
                        selected = "DMPs"
                    ),
                    downloadButton("download_export_heatmaps"),
                    p("Press to download the custom heatmap in the gplots::heatmap.2 version.")
                )
            ),
            tabItem(
                tabName = "help",
                fluidPage(h1("Help"))
            )
        )
    )
))
