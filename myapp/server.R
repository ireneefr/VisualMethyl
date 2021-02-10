
library(shiny)

shinyServer(function(input, output, session) {
    
    # Max size
    options(shiny.maxRequestSize = 8000 * 1024^2) # 5MB getShinyOption("shiny.maxRequestSize") | 30*1024^2 = 30MB

    
    # Reaction of home action buttons
    observeEvent(input$b_qc, {
        newtab <- switch(input$menu,
                         "home" = "qc",
                         "qc" = "home")
        updateTabItems(session, "menu", newtab)
        })
    observeEvent(input$b_exploratory_analysis, {
        newtab <- switch(input$menu,
                         "home" = "exploratory_analysis",
                         "exploratory_analysis" = "home")
        updateTabItems(session, "menu", newtab)
    })
    observeEvent(input$b_dmp_dmr, {
        newtab <- switch(input$menu,
                         "home" = "dmp_dmr",
                         "dmp_dmr" = "home")
        updateTabItems(session, "menu", newtab)
    })
    observeEvent(input$b_functional_enrichment, {
        newtab <- switch(input$menu,
                         "home" = "functional_enrichment",
                         "functional_enrichment" = "home")
        updateTabItems(session, "menu", newtab)
    })
    observeEvent(input$b_survival, {
        newtab <- switch(input$menu,
                         "home" = "survival",
                         "survival" = "home")
        updateTabItems(session, "menu", newtab)
    })
    observeEvent(input$b_predicted_models, {
        newtab <- switch(input$menu,
                         "home" = "predicted_models",
                         "predicted_models" = "home")
        updateTabItems(session, "menu", newtab)
    })
    observeEvent(input$b_external_sources, {
        newtab <- switch(input$menu,
                         "home" = "external_sources",
                         "external_sources" = "home")
        updateTabItems(session, "menu", newtab)
    })
    observeEvent(input$b_genome_browser, {
        newtab <- switch(input$menu,
                         "home" = "genome_browser",
                         "genome_browser" = "home")
        updateTabItems(session, "menu", newtab)
    })

    
    # Input button
    output$ui_input_data <- renderUI({
        if (!is.null(input$input_data$datapath)) {
            return(actionButton("b_input_data", "Load Data"))
        } else {
            return()
        }
    })
    
    # Enable load button every time file input is updated
    observeEvent(input$input_data, shinyjs::enable("b_input_data"))
    
})
