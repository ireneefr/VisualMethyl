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

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)

shinyUI(
    dashboardPage(

    dashboardHeader(title = "APP"),

    dashboardSidebar(
        
        sidebarMenu(id = "menu",
                    menuItem("Home", tabName = "home"),
                    menuItem("QC", tabName = "qc"),
                    menuItem("Exploratory Analysis", tabName = "exploratory_analysis"),
                    menuItem("DMPs/DMRs", tabName = "dmp_dmr"),
                    menuItem("Functional Enrichment", tabName = "functional_enrichment"),
                    menuItem("Survival", tabName = "survival"),
                    menuItem("Predicted Models", tabName = "predicted_models"),
                    menuItem("External Sources", tabName = "external_sources"),
                    menuItem("Genome Browser", tabName = "genome_browser")
        )
    ),
    
    dashboardBody(
        
        tabItems(
            
            tabItem(
                tabName = "home",
                verticalLayout( 
                    boxPlus(title = "INPUT DATA", width = 12, closable = FALSE, collapsible = TRUE, 
                            fileInput("input_data", "Upload input data (.zip)", multiple = FALSE, accept = ".zip"),
                            uiOutput("ui_input_data")
                    ),
                    fluidRow(style = "margin-left: 2px; margin-right: 2px",
                    conditionalPanel(
                        "typeof output.samples_table != 'undefined'",
                        boxPlus(title = "SELECTION OPTIONS", id = "selection_options",
                                closable = FALSE, collapsible = TRUE,
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
                            actionButton("button_input_next", "Continue to Analysis")
                        )
                    ),
                    
                    # Box 2: Table
                    conditionalPanel(
                        "input.b_input_data > 0", boxPlus(title = "SAMPLES TABLE",
                                                          closable = FALSE, collapsible = TRUE,
                        width = 9,
                        withSpinner(DT::DTOutput("samples_table"))
                        
                    ))),
                    conditionalPanel("input.button_input_next > 0",
                    splitLayout(style = "margin-left: 14px; margin-right: 14px;", cellArgs = list(style = "padding: 4px;"),
                        actionButton("b_qc", width = "100%", label = strong("QC"), class = "btn-primary", style = "color:#fff; padding:50px; font-size:300%", disable = TRUE),
                        actionButton("b_exploratory_analysis", width = "100%",label = strong(HTML("Exploratory <br/> Analysis")), class = "btn-primary", style = "color:#fff; padding:20px; font-size:300%"),
                        actionButton("b_dmp_dmr", width = "100%", label = strong("DMPs/DMRs"), class = "btn-primary", style = "color:#fff; padding:50px; font-size:300%"),
                        actionButton("b_functional_enrichment", width = "90%",label = strong(HTML("Functional <br/> Enrichment")), class = "btn-primary", style = "color:#fff; padding:20px; font-size:300%")
                    ),
                    splitLayout(style = "margin-left: 14px; margin-right: 14px", cellArgs = list(style = "padding: 4px;"),
                        actionButton("b_survival", width = "100%", label = strong("Survival"), class = "btn-primary", style = "color:#fff; padding:50px; font-size:300%"),
                        actionButton("b_predicted_models", width = "100%", label = strong(HTML("Predicted <br/> Models")), class = "btn-primary", style = "color:#fff; padding:20px; font-size:300%"),
                        actionButton("b_external_sources", width = "100%", label = strong(HTML("External <br/> Sources")), class = "btn-primary", style = "color:#fff; padding:20px; font-size:300%"),
                        actionButton("b_genome_browser", width = "90%", label = strong(HTML("Genome <br/> Browser")), class = "btn-primary", style = "color:#fff; padding:20px; font-size:300%")
                    ))
                    
                )
            ),
            tabItem(
                tabName = "qc",
                verticalLayout(
                    boxPlus(title = "NORMALIZATION OPTIONS",
                            collapsible = FALSE, closable = FALSE, width = 12,
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
                            h4()
                    ),
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
                        textOutput("text_minfi_probes")
                    ),
                    
                    
                    mainPanel(
                        width = 9,
                        tabsetPanel(
                            
                            ###################################################################################
                            
                            tabPanel(
                                "Control Type",
                                selectInput("controlType", "Choose a control type:",
                                            choices = c(
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
                                            ), selected = 1),
                                selectInput("select_slide", "Select slide:", choices = c()),
                                withSpinner(plotOutput("controlTypePlotGreen")),
                                withSpinner(plotOutput("controlTypePlotRed"))
                            ),
                            
                            ###################################################################################
                            
                            tabPanel(
                                "Quality Control",
                                h4("Overall Signal"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_qcraw")),
                                h4("Bisulfite Conversion II"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_bisulfiterawII"))
                            ),
                            
                            tabPanel(
                                "Density plot",
                                selectInput("probeType", "Choose a probe type for the density curves:",
                                            choices = c("I-Green","I-Red","II"),
                                            selected="I-Green"),
                                textOutput("prueba1"),
                                textOutput("prueba2"),
                                textOutput("prueba3"),
                                textOutput("prueba4"),
                                h4("Raw"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_densityplotraw")),
                                h4("Processed"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_densityplot"))
                            ),
                            tabPanel(
                                "Failed probes",
                                h4("Failure Rate Plot"),
                                textOutput("text"),
                                withSpinner(plotly::plotlyOutput("failure_rate_plot"))
                            ),
                            tabPanel(
                                "Sex prediction",
                                h4("X vs Y chromosomes signal plot"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_sex")),
                                withSpinner(DT::DTOutput("table_minfi_sex"))
                            ),
                            tabPanel(
                                "Correlations",
                                h4("Processed"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_corrplot")),
                                selectInput("select_minfi_typecorrplot", "Select data to plot", choices = c("p.value", "correlation value"), selected = "correlation value"),
                                withSpinner(DT::DTOutput("table_minfi_corrplot"))
                            )
                        )
                    )
                )
            ),
            tabItem(
                tabName = "exploratory_analysis",
                fluidPage(h1("Exploratory Analysis"))
            ),
            tabItem(
                tabName = "dmp_dmr",
                fluidPage(h1("DMPs/DMRs"))
            ),
            tabItem(
                tabName = "functional_enrichment",
                fluidPage(h1("Functional Enrichment"))
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
            )
        )
    )
))
