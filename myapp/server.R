
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
    
    # input_data: enable b_input_data
    observeEvent(input$input_data, shinyjs::enable("b_input_data"))
    
    ################################
    ########## rval_sheet() METHARRAY SHEET ##########
    
    # When you press button_input_load, the data is unzipped and the metharray sheet is loaded
    rval_sheet <- eventReactive(input$b_input_data, {
        print("OK")
        # Check if updated file is .zip
        validate(need(tools::file_ext(input$input_data$datapath) == "zip", "File extension should be .zip"))
        
        shinyjs::disable("button_input_load") # disable the load button to avoid multiple clicks
        
        if (dir.exists(paste0(tempdir(), "/experiment_data"))) {
            unlink(paste0(tempdir(), "/experiment_data"), recursive = TRUE) # remove current files in target directory
        }
        
        zip::unzip(input$input_data$datapath,
                   exdir = paste0(tempdir(), "/experiment_data")
        ) # extracting zip
        
        sheet <- minfi::read.metharray.sheet(paste0(tempdir(), "/experiment_data"))
        
        # We check if sheet is correct
        # This is to prevent app crashes when zip updated is not correct.
        validate(
            need(
                is.data.frame(sheet) &
                    any(colnames(sheet) %in% "Slide") &
                    any(colnames(sheet) %in% "Array"),
                "SampleSheet is not correct. Please, check your samplesheet and your zip file."
            )
        )
        
        validate(
            need(
                anyDuplicated(colnames(sheet)) == 0,
                "Repeated variable names are not allowed. Please, modify your sample sheet."
            )
        )
        
        colnames(sheet) <- make.names(colnames(sheet)) # fix possible not-valid colnames
        
        sheet
    })
    
    
    rval_sheet_target <- eventReactive(
        input$button_input_next,
        rval_sheet()[rval_sheet()[[input$select_input_samplenamevar]] %in% input$selected_samples, ]
    )
    
    
    rval_clean_sheet_target <- eventReactive(rval_gset(), {
        generate_clean_samplesheet(target_samplesheet = minfi::pData(rval_gset()),
                                   donorvar = input$select_input_donorvar)
        
    })
    
    
    
    
    # When you press button_input_load, the form options are updated
    observeEvent(input$b_input_data, {
        updateSelectInput(
            session,
            "select_input_samplenamevar",
            label = "Select Sample Names Column:",
            choices = colnames(rval_sheet())
        )
        
        updateSelectInput(
            session,
            "select_input_groupingvar",
            label = "Select Variable of Interest:",
            choices = colnames(rval_sheet())
        )
        updateSelectInput(
            session,
            "select_input_donorvar",
            label = "Select Donor Variable:",
            choices = colnames(rval_sheet())
        )
        updateSelectInput(
            session,
            "select_input_sex",
            label = "Select Sex Column",
            choices = c("None", colnames(rval_sheet()))
        )
        
        
        shinyjs::enable("button_input_next") # Enable button continue
    })
    
    
    # The checkbox of samples to process is updated when samplenamevar changes
    observeEvent({input$select_input_samplenamevar
        input$select_input_groupingvar},
        updatePickerInput(
            session,
            "selected_samples",
            label = "Select Samples to Process:",
            selected = rval_sheet()[, input$select_input_samplenamevar],
            choices = rval_sheet()[, input$select_input_samplenamevar],
            choicesOpt = list(subtext = paste("Group: ", rval_sheet()[, input$select_input_groupingvar]))
        )
    )
    
    # when samples selected are changed, continue button is enabled again
    observeEvent(input$selected_samples, shinyjs::enable("button_input_next"))
    
    # The dataframe is rendered
    output$samples_table <- DT::renderDT(
        rval_sheet(),
        rownames = FALSE,
        selection = "single",
        style = "bootstrap",
        options = list(
            pageLength = 10,
            autoWidth = TRUE,
            scrollX = TRUE,
            columnDefs = list(list(
                targets = match("Basename", colnames(rval_sheet())) - 1, visible = FALSE
            ))
        )
    )
    
    
    
    ########## rval_rgset() RGCHANNEL ##########
    
    # rval_rgset loads RGSet using read.metharray.exp and the sample sheet (rval_sheet())
    rval_rgset <- eventReactive(input$button_input_next, ignoreNULL = FALSE, {
        validate(need(input$fileinput_input != "", "Data has not been uploaded yet"))
        
        # Prior check to test variable selection
        if (anyDuplicated(rval_sheet_target()[, input$select_input_samplenamevar]) > 0 |
            anyDuplicated(rval_sheet_target()[, input$select_input_groupingvar]) == 0) {
            showModal(
                modalDialog(
                    title = "Variable error",
                    "Check if selected variables are correct. Sample Name Variable should not have duplicated values
          and the variable of interest should have groups greater than 1.",
                    easyClose = TRUE,
                    footer = NULL
                )
            )
        }
        
        # Check prior conditions to read data
        validate(need(
            anyDuplicated(rval_sheet_target()[, input$select_input_samplenamevar]) == 0,
            "Sample Name Variable should not have duplicated values"
        ))
        validate(need(
            anyDuplicated(rval_sheet_target()[, input$select_input_groupingvar]) > 0,
            "Grouping variable should have groups greater than 1"
        ))
        
        # disable button to avoid multiple clicks
        shinyjs::disable("button_input_next")
        
        
        # We need to check if this step works
        withProgress(
            message = "Reading array data...",
            value = 2,
            max = 5,
            {
                try({
                    RGSet <- read_idats(
                        targets = rval_sheet_target())
                })
                
                if (!exists("RGSet", inherits = FALSE)) {
                    showModal(
                        modalDialog(
                            title = "reading error",
                            "Minfi can't read arrays specified in your samplesheet. Please, check your zipfile and your sampleSheet",
                            easyClose = TRUE,
                            footer = NULL
                        )
                    )
                    shinyjs::disable("button_minfi_select")
                }
                
                validate(
                    need(
                        exists("RGSet", inherits = FALSE),
                        "Minfi can't read arrays specified in your samplesheet. Please, check your zipfile and your sampleSheet"
                    )
                )
                
                #Checking array type and annotation
                nProbes = length(minfi::featureNames(RGSet))
                
                if(nProbes >= 622000 & nProbes <= 623000){
                    
                    if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE) |
                        !requireNamespace("IlluminaHumanMethylation450kmanifest", quietly = TRUE))
                    {
                        showModal(
                            modalDialog(
                                title = "Missing package(s)",
                                "450k annotation or manifest packages are not available. Please, install IlluminaHumanMethylation450kmanifest and IlluminaHumanMethylation450kanno.ilmn12.hg19 packages and restart the application.",
                                easyClose = TRUE,
                                footer = NULL
                            )
                        )
                    }
                    
                    validate(
                        need(
                            requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE) &
                                requireNamespace("IlluminaHumanMethylation450kmanifest", quietly = TRUE),
                            "450k annotation or manifest packages are not available. Please, install IlluminaHumanMethylation450kmanifest and IlluminaHumanMethylation450kanno.ilmn12.hg19 packages."
                        )
                    )
                }
                else if (nProbes >= 1032000 & nProbes <= 1053000){
                    
                    if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE) |
                        !requireNamespace("IlluminaHumanMethylationEPICmanifest", quietly = TRUE))
                    {
                        showModal(
                            modalDialog(
                                title = "Missing package(s)",
                                "EPIC annotation or manifest packages are not available. Please, install IlluminaHumanMethylationEPICmanifest and IlluminaHumanMethylationEPICanno.ilm10b4.hg19 packages and restart the application.",
                                easyClose = TRUE,
                                footer = NULL
                            )
                        )
                    }
                    
                    validate(
                        need(
                            requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE) &
                                requireNamespace("IlluminaHumanMethylationEPICmanifest", quietly = TRUE),
                            "EPIC annotation or manifest packages are not available. Please, install IlluminaHumanMethylationEPICmanifest and IlluminaHumanMethylationEPICanno.ilm10b4.hg19 packages."
                        )
                    )
                }
                
                # analysis restarted
                rval_analysis_finished(FALSE)
                rval_dmrs_finished(FALSE)
                
                # we return RGSet
                RGSet
            }
        )
    })
    
    # We change the page to the next one
    observeEvent(input$button_input_next, {
        # check if rgset is loaded
        req(rval_rgset())
        
        updateSelectInput(
            session,
            "select_slide",
            choices = rval_sheet_target()$Slide
        )
        
        updateSelectInput(
            session,
            "select_control",
            selected = "BISULFITE CONVERSION I"
        )
        
        # update PCA parameters
        updateSelectInput(
            session,
            "select_minfi_pcaplot_pcx",
            choices = paste0("PC", seq_len(nrow(rval_sheet_target()))),
            selected = "PC1"
        )
        
        updateSelectInput(
            session,
            "select_minfi_pcaplot_pcy",
            choices = paste0("PC", seq_len(nrow(rval_sheet_target()))),
            selected = "PC2"
        )
        
        updateSelectInput(
            session,
            "select_minfi_pcaplot_color",
            choices = c(
                colnames(rval_sheet_target()),
                "xMed",
                "yMed",
                "predictedSex"
            ),
            
            selected = input$select_input_groupingvar
        )
        
        shinyjs::enable("button_minfi_select")
        #updateTabItems(session, "menu", "normalization")
    })
    
    
    
    
    
})
