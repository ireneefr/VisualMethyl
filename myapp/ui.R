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

shinyUI(
    dashboardPage(

    dashboardHeader(title = "APP"),

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
    
    dashboardBody(
        
        tabItems(
            
            tabItem(
                tabName = "data",
                verticalLayout( 
                    boxPlus(title = "INPUT DATA", width = 12, closable = FALSE, collapsible = TRUE, 
                            fileInput("input_data", "Upload input data (.zip)", multiple = FALSE, accept = ".zip"),
                            uiOutput("ui_input_data")
                    ),
                    fluidRow(style = "margin-left: 2px; margin-right: 2px",
                    conditionalPanel(
                        "input.b_input_data > 0",
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
                        
                    )))
                    
                )
            ),
            tabItem(
                tabName = "analysis",
                verticalLayout(
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
                        downloadButton("download_export_markdown")
                    ),
                    
                    
                    mainPanel(
                        width = 9,
                        verticalLayout(
                        boxPlus(title = "FAILED PROBES", width = "100%", closable = FALSE, collapsible = TRUE, collapsed = TRUE,
                                h4("Failure Rate Plot"), 
                                withSpinner(plotly::plotlyOutput("failure_rate_plot"))
                        ),
                        boxPlus(title = "CONTROL TYPES", width = "100%", closable = FALSE, collapsible = TRUE, collapsed = TRUE,
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
                        boxPlus(title = "PREDICTED SEX", width = "100%", closable = FALSE, collapsible = TRUE, collapsed = TRUE,
                                h4("X vs Y chromosomes signal plot"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_sex")),
                                withSpinner(DT::DTOutput("table_minfi_sex"))
                        ),
                        boxPlus(title = "DENSITY PLOT", width = "100%", closable = FALSE, collapsible = TRUE, collapsed = TRUE,
                                selectInput("probeType", "Choose a probe type for the density curves:",
                                            choices = c("I-Green","I-Red","II"),
                                            selected="I-Green"),
                                h4("Raw"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_densityplotraw")),
                                h4("Processed"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_densityplot"))
                        ),
                        boxPlus(title = "SNP ANALYSIS", width = "100%", closable = FALSE, collapsible = TRUE, collapsed = TRUE,
                                h4("SNPs beta-values (Raw)"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_snps"))
                        ),
                        boxPlus(title = "BATCH EFFECTS", width = "100%", closable = FALSE, collapsible = TRUE, collapsed = TRUE,
                                h4("Processed"),
                                withSpinner(plotly::plotlyOutput("graph_minfi_corrplot")),
                                selectInput("select_minfi_typecorrplot", "Select data to plot", choices = c("p.value", "correlation value"), selected = "correlation value"),
                                withSpinner(DT::DTOutput("table_minfi_corrplot"))
                        ))
                    )
                )
            ),
            tabItem(
                tabName = "exploratory_analysis",
                fluidPage(
                    sidebarPanel(
                        width = 3,
                        h1("OPTIONS")
                    ),
                    mainPanel(
                        widh = 9,
                        verticalLayout(
                            boxPlus(title = "DENSITY PLOT", width = "100%", closable = FALSE, collapsible = TRUE, collapsed = TRUE
                            ),
                            boxPlus(title = "VIOLIN PLOT", width = "100%", closable = FALSE, collapsible = TRUE, collapsed = TRUE,
                                    withSpinner(plotOutput("graph_violin"))
                            ),
                            boxPlus(title = "PRINCIPAL COMPONENT ANALYSIS", width = "100%", closable = FALSE, collapsible = TRUE, collapsed = TRUE,
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
                            boxPlus(title = "HEATMAP", width = "100%", closable = FALSE, collapsible = TRUE, collapsed = TRUE,
                                    h1("heatmap")
                            ),
                            boxPlus(title = "DECONVOLUTION", width = "100%", closable = FALSE, collapsible = TRUE, collapsed = TRUE,
                                    h1("deconvolution")
                            ),
                            boxPlus(title = "AGE METH", width = "100%", closable = FALSE, collapsible = TRUE, collapsed = TRUE,
                                    h1("age meth")
                            ),
                            boxPlus(title = "GENOME ANALYSIS PLOT", width = "100%", closable = FALSE, collapsible = TRUE, collapsed = TRUE,
                                    h1("genome analysis plot")
                            ),
                            boxPlus(title = "HYPO/HYPER", width = "100%", closable = FALSE, collapsible = TRUE, collapsed = TRUE,
                                    h1("hypo/hyper based on betas")
                            ),
                            boxPlus(title = "CIRCOS", width = "100%", closable = FALSE, collapsible = TRUE, collapsed = TRUE,
                                    h1("circos")
                            )
                        )
                    )
                )
            ),
            tabItem(
                tabName = "dmp_dmr",
                fluidPage(
                    tabsetPanel(
                        tabPanel("DMPs", sidebarLayout(
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
                                boxPlus(title = "DMP TABLE", width = "100%", closable = TRUE, collapsible = TRUE, collapsed = TRUE
                                        
                                )
                            )
                        )),
                        tabPanel("DMRs", fluidPage())
                    )
                )
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
            ),
            tabItem(
                tabName = "export",
                fluidPage(h1("Export"))
            ),
            tabItem(
                tabName = "help",
                fluidPage(h1("Help"))
            )
        )
    )
))
