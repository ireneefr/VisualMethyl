
library(shiny)
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
                    splitLayout(style = "margin-left: 14px; margin-right: 14px;", cellArgs = list(style = "padding: 4px;"),
                        actionButton("b_qc", width = "100%", label = strong("QC"), class = "btn-primary", style = "color:#fff; padding:80px; font-size:300%"),
                        actionButton("b_exploratory_analysis", width = "100%",label = strong(HTML("Exploratory <br/> Analysis")), class = "btn-primary", style = "color:#fff; padding:50px; font-size:300%"),
                        actionButton("b_dmp_dmr", width = "100%", label = strong("DMPs/DMRs"), class = "btn-primary", style = "color:#fff; padding:80px; font-size:300%"),
                        actionButton("b_functional_enrichment", width = "90%",label = strong(HTML("Functional <br/> Enrichment")), class = "btn-primary", style = "color:#fff; padding:50px; font-size:300%")
                    ),
                    splitLayout(style = "margin-left: 14px; margin-right: 14px", cellArgs = list(style = "padding: 4px;"),
                        actionButton("b_survival", width = "100%", label = strong("Survival"), class = "btn-primary", style = "color:#fff; padding:80px; font-size:300%"),
                        actionButton("b_predicted_models", width = "100%", label = strong(HTML("Predicted <br/> Models")), class = "btn-primary", style = "color:#fff; padding:50px; font-size:300%"),
                        actionButton("b_external_sources", width = "100%", label = strong(HTML("External <br/> Sources")), class = "btn-primary", style = "color:#fff; padding:50px; font-size:300%"),
                        actionButton("b_genome_browser", width = "90%", label = strong(HTML("Genome <br/> Browser")), class = "btn-primary", style = "color:#fff; padding:50px; font-size:300%")
                    )
                    
                )
            ),
            tabItem(
                tabName = "qc",
                fluidPage(h1("QC"))
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
