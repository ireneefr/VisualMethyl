
library(shiny)
source("utils_analysis.R")
source("utils_graphs.R")
source("utils_download.R")
library(ggplot2)


shinyServer(function(input, output, session) {
    # INITIALIZE REACTIVE VARIABLES
    rval_gset_done <- reactiveVal(value = FALSE)
    rval_generated_limma_model <- reactiveVal(value = FALSE)
    rval_analysis_finished <- reactiveVal(value = FALSE)
    rval_gset_getBeta_done <- reactiveVal(value = FALSE)
    rval_filteredlist2heatmap_valid <- reactiveVal(value = FALSE)
    rval_filteredmcsea2heatmap_valid <- reactiveVal(value = FALSE)
    rval_dmrs_finished <- reactiveVal(value = FALSE)
    rval_dmrs_ready2heatmap <- reactiveVal(value = FALSE)
    rval_dmrs_ready2mcsea <- reactiveVal(value = FALSE)
    
    
    # Max size
    options(shiny.maxRequestSize = 8000 * 1024^2) # 5MB getShinyOption("shiny.maxRequestSize") | 30*1024^2 = 30MB
    n_cores <- parallel::detectCores() / 2
    
    
    output$arrow1 <- renderImage(list(src = "img/arrow3.png", contentType = "image/png"), deleteFile = FALSE)
    output$arrow2 <- renderImage(list(src = "img/arrow3.png", contentType = "image/png"), deleteFile = FALSE)
    output$arrow3 <- renderImage(list(src = "img/arrow4.png", contentType = "image/png"), deleteFile = FALSE)
    output$arrow4 <- renderImage(list(src = "img/arrow4.png", contentType = "image/png"), deleteFile = FALSE)
    
    output$arrow1_ <- renderImage(list(src = "img/arrow3.png", contentType = "image/png"), deleteFile = FALSE)
    output$arrow3_ <- renderImage(list(src = "img/arrow4.png", contentType = "image/png"), deleteFile = FALSE)
    output$arrow4_ <- renderImage(list(src = "img/arrow4.png", contentType = "image/png"), deleteFile = FALSE)
    
    
    
    observe({
        shinyjs::disable("check_qc")
        shinyjs::disable("check_group_qc")
        shinyjs::disable("check_exploratory_analysis")
        shinyjs::disable("check_group_exploratory_analysis")
        shinyjs::disable("check_dmps")
        shinyjs::disable("check_group_dmps")
        shinyjs::disable("check_dmrs")
        shinyjs::disable("check_group_dmrs")
        shinyjs::disable("check_functional_enrichment")
        shinyjs::disable("check_group_functional_enrichment")
        })
    observeEvent(rval_gset(), {
        shinyjs::enable("check_qc")
        shinyjs::enable("check_exploratory_analysis")
        updateCheckboxInput(session, "check_qc", value = TRUE)
        updateCheckboxInput(session, "check_exploratory_analysis", value = TRUE)
        })
    observeEvent(rval_filteredlist(), {
        shinyjs::enable("check_dmps")
        shinyjs::enable("check_functional_enrichment")
        updateCheckboxInput(session, "check_dmps", value = TRUE)
        updateCheckboxInput(session, "check_functional_enrichment", value = TRUE)
    })
    observeEvent(rval_mcsea(), {
        shinyjs::enable("check_dmrs")
        updateCheckboxInput(session, "check_dmrs", value = TRUE)
    })
    observe({
        req(rval_gset())
        if(input$check_qc){
            updateCheckboxGroupInput(session, "check_group_qc",choices = list("Intensities boxplots" = 1,
                                                                              "Failed probes" = 2,
                                                                              "Density plots" = 3,
                                                                              "SNPs heatmap" = 4,
                                                                              "Sex prediction" = 5,
                                                                              "Batch effects" = 6),  selected = c(1:6))
            shinyjs::enable("check_group_qc")
        } else{
            shinyjs::disable("check_group_qc")
            updateCheckboxGroupInput(session, "check_group_qc", choices = c(), selected = 0)
        }
        if(input$check_exploratory_analysis){
            updateCheckboxGroupInput(session, "check_group_exploratory_analysis",choices = list("Violin plots" = 1,
                                                                                                "Principal Component Analysis" = 2,
                                                                                                "Heatmaps" = 3,
                                                                                                "Deconvolution" = 4,
                                                                                                "Age methylation" = 5,
                                                                                                "Hypo/Hyper" = 6), selected = c(1:6))
            shinyjs::enable("check_group_exploratory_analysis")
        } else{
            updateCheckboxGroupInput(session, "check_group_exploratory_analysis", choices = c(), selected = 0)
            shinyjs::disable("check_group_exploratory_analysis")

        }
    })
    observe({
        req(rval_filteredlist())
        if(input$check_dmps){
            updateCheckboxGroupInput(session, "check_group_dmps", choices = list("Table" = 1,
                                                                           "Heatmap" = 2,
                                                                           "Annotation" = 3,
                                                                           "Manhattan plot" = 4,
                                                                           "Volcano plot" = 5), selected = c(1:5))
            shinyjs::enable("check_group_dmps")
        } else{
            updateCheckboxGroupInput(session, "check_group_dmps", choices = c(), selected = 0)
            shinyjs::disable("check_group_dmps")
            
        }
        if(input$check_functional_enrichment){
            updateCheckboxGroupInput(session, "check_group_functional_enrichment", choices = list("Kegg" = 1,
                                                                                                  "Gene Ontology (GO)" = 2,
                                                                                                  "Reactome" = 3), selected = c(1:3))
            shinyjs::enable("check_group_functional_enrichment")
        } else{
            updateCheckboxGroupInput(session, "check_group_functional_enrichment", choices = c(), selected = 0)
            shinyjs::disable("check_group_functional_enrichment")
        }
    })
    observe({
        req(rval_mcsea())
        if(input$check_dmrs){
            updateCheckboxGroupInput(session, "check_group_dmrs", choices = list("Table" = 1,
                                                                                 "Heatmap" = 2,
                                                                                 "Annotation" = 3), selected = c(1:3))
            shinyjs::enable("check_group_dmrs")
        } else{
            updateCheckboxGroupInput(session, "check_group_dmrs", c(), selected = 0)
            shinyjs::disable("check_group_dmrs")
        }
    })
    
    
    # Reaction of data action buttons
    observeEvent(input$b_qc, {
        newtab <- switch(input$menu,
                         "analysis" = "qc",
                         "qc" = "analysis")
        updateTabItems(session, "menu", newtab)
        })
    observeEvent(input$b_exploratory_analysis, {
        newtab <- switch(input$menu,
                         "analysis" = "exploratory_analysis",
                         "exploratory_analysis" = "analysis")
        updateTabItems(session, "menu", newtab)
    })
    observeEvent(input$b_dmp_dmr, {
        newtab <- switch(input$menu,
                         "analysis" = "dmp_dmr",
                         "dmp_dmr" = "analysis")
        updateTabItems(session, "menu", newtab)
    })
    observeEvent(input$b_functional_enrichment, {
        newtab <- switch(input$menu,
                         "analysis" = "functional_enrichment",
                         "functional_enrichment" = "analysis")
        updateTabItems(session, "menu", newtab)
    })
    observeEvent(input$b_survival, {
        newtab <- switch(input$menu,
                         "analysis" = "survival",
                         "survival" = "analysis")
        updateTabItems(session, "menu", newtab)
    })
    observeEvent(input$b_predicted_models, {
        newtab <- switch(input$menu,
                         "analysis" = "predicted_models",
                         "predicted_models" = "analysis")
        updateTabItems(session, "menu", newtab)
    })
    observeEvent(input$b_external_sources, {
        newtab <- switch(input$menu,
                         "analysis" = "external_sources",
                         "external_sources" = "analysis")
        updateTabItems(session, "menu", newtab)
    })
    observeEvent(input$b_genome_browser, {
        newtab <- switch(input$menu,
                         "analysis" = "genome_browser",
                         "genome_browser" = "analysis")
        updateTabItems(session, "menu", newtab)
    })
    observeEvent(input$button_input_next, {
        newtab <- switch(input$menu,
                         "data" = "analysis",
                         "analysis" = "data")
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
        updateSelectInput(
            session,
            "select_input_age",
            label = "Select Age Column",
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
        validate(need(input$input_data != "", "Data has not been uploaded yet"))

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
        
        showModal(modalDialog(
            title = NULL, footer = NULL,
            div(
                img(src="https://upload.wikimedia.org/wikipedia/commons/7/7d/Pedro_luis_romani_ruiz.gif"),
                p("Reading array data..."),
                style = "margin: auto; text-align: center"
            )
        ))
        
        # We need to check if this step works
        withProgress(
            message = "Reading array data...",
            value = 2,
            max = 5,
            {
                try({
                    RGSet <- minfi::read.metharray.exp(targets = rval_sheet_target(), verbose = TRUE, force = TRUE) #Read idats
                })
                
                if (!exists("RGSet", inherits = FALSE)) {
                    removeModal()
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
                        removeModal()
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
                        removeModal()
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
                removeModal()
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
        
        updatePickerInput(
            session,
            "selected_samples_h",
            label = "Select Samples to Plot:",
            selected = rval_sheet_target()[, input$select_input_samplenamevar],
            choices = rval_sheet_target()[, input$select_input_samplenamevar],
            choicesOpt = list(subtext = paste("Group: ", rval_sheet_target()[, input$select_input_groupingvar]))
        )
        print("ok")
        shinyjs::enable("button_minfi_select")
    })
    
    
    ########## rval_gset() NORMALIZATION ##########
    
    rval_gset <- eventReactive(input$button_minfi_select, {
        validate(need(
            !is.null(rval_rgset()),
            "Raw data has not been loaded yet."
        ))
        
        shinyjs::disable("button_minfi_select") # disable button to avoid repeat clicking

        showModal(modalDialog(
            title = NULL, footer = NULL,
            div(
                img(src="https://upload.wikimedia.org/wikipedia/commons/7/7d/Pedro_luis_romani_ruiz.gif"),
                p("Normalization in progress..."),
                style = "margin: auto; text-align: center"
            )
        ))
        
        withProgress(
            message = "Normalization in progress...",
            value = 1,
            max = 4,
            {
                try({
                    gset <- normalize_rgset(
                        rgset = rval_rgset(), normalization_mode = input$select_minfi_norm,
                        detectP = 0.01, dropSNPs = input$select_minfi_dropsnps, maf = input$slider_minfi_maf,
                        dropCpHs = input$select_minfi_dropcphs, dropSex = input$select_minfi_chromosomes
                    )
                })
                
                # check if normalization has worked
                
                if (!exists("gset", inherits = FALSE)) {
                    removeModal()
                    showModal(
                        modalDialog(
                            title = "Normalization failure",
                            "An unexpected error has occurred during minfi normalization. Please, notify the error to the package maintainer.",
                            easyClose = TRUE,
                            footer = NULL
                        )
                    )
                    
                    shinyjs::enable("button_minfi_select")
                }
                
                validate(
                    need(
                        exists("gset", inherits = FALSE),
                        "An unexpected error has occurred during minfi normalization. Please, notify the error to the package maintainer."
                    )
                )
                # enable button
                rval_gset_done(TRUE)
                shinyjs::enable("button_minfi_select")
                removeModal()
                # return gset
                gset
            }
        )
        
        
    })
    
    
    
    # filtered probes info
    
    rval_gsetprobes <- eventReactive(input$button_minfi_select, {
        req(rval_gset())
        length(minfi::featureNames(rval_gset()))
    })
    
    output$text_minfi_probes <- renderText(paste(rval_gsetprobes(), "positions after normalization"))
    
    ########## B & M VALUES ######################
    
    rval_rgset_getBeta <- eventReactive(rval_rgset(), {
        bvalues <- as.data.frame(minfi::getBeta(rval_rgset()))
        colnames(bvalues) <- rval_sheet_target()[[input$select_input_samplenamevar]]
        bvalues
    })
    
    rval_gset_getBeta <- eventReactive(rval_gset(), {
        bvalues <- as.data.frame(minfi::getBeta(rval_gset()))
        colnames(bvalues) <- rval_sheet_target()[[input$select_input_samplenamevar]]
        rval_gset_getBeta_done(TRUE)
        print("BETA")
        print(bvalues)
        print(class(bvalues))
        bvalues
    })
    
    rval_gset_getM <- reactive({
        mvalues <- minfi::getM(rval_gset())
        colnames(mvalues) <- rval_sheet_target()[[input$select_input_samplenamevar]]
        mvalues
    }) 
    
    
    ##### INTENSITIES BOXPLOTS #####

    boxplot_intensities_green <- reactive(create_boxplot_intensities_green(rval_rgset()))
    boxplot_intensities_red <- reactive(create_boxplot_intensities_red(rval_rgset()))
    output$green_intensities_plot <- renderPlot(boxplot_intensities_green())
    output$red_intensities_plot <- renderPlot(boxplot_intensities_red())
    
    ########## FAILURE RATE PLOT ##########
    
    failure_plot <- reactive(create_failure(rval_rgset(), rval_rgset_getBeta()))
    output$failure_rate_plot <- plotly::renderPlotly(failure_plot()[["graph"]])
    output$failure_rate_table <- DT::renderDT(failure_plot()[["info"]])
    
    ########## CONTROL TYPE PLOTS ##########
    
    control_type <- reactive(create_control_type(rval_rgset(), rval_sheet_target(), input$controlType, input$select_slide))
    output$controlTypePlotGreen <- renderPlot(control_type()[["green"]])
    output$controlTypePlotRed <- renderPlot(control_type()[["red"]])
    
    
    ###### SEX PREDICTION #####
    
    # Sex prediction
    
    rval_plot_sexprediction <- reactive({
        req(rval_gset())
        create_pred_sexplot(rval_gset(), rval_sheet_target()[, input$select_input_samplenamevar])
    })
    
    rval_plot_sextable <- reactive({
        req(rval_gset())
        if(input$select_input_sex == "None"){
            data.frame(name = rval_sheet_target()[[input$select_input_samplenamevar]], predictedSex = as.data.frame(minfi::pData(rval_gset()))[["predictedSex"]])
        }
        else{
        data.frame(name = rval_sheet_target()[[input$select_input_samplenamevar]], sex = rval_sheet_target()[[input$select_input_sex]], predictedSex = as.data.frame(minfi::pData(rval_gset()))[["predictedSex"]])
        }   
    })
    
    output$graph_minfi_sex <- plotly::renderPlotly(rval_plot_sexprediction())
    
    output$table_minfi_sex <- DT::renderDT(
        rval_plot_sextable(),
        rownames = FALSE,
        selection = "single",
        style = "bootstrap",
        caption = "Predicted sex:",
        options = list(
            pageLength = 10,
            scrollX = TRUE,
            autoWidth = TRUE
        )
    )
    
    ########## DENSITY PLOT #####################
    
    channel <- reactive(getProbeInfo(rval_rgset(), type = input$probeType)[, "Name"])
    
    ###
    channel_green <- reactive(getProbeInfo(rval_rgset(), type = "I-Green")[, "Name"])
    channel_red <- reactive(getProbeInfo(rval_rgset(), type = "I-Red")[, "Name"])
    channel_II <- reactive(getProbeInfo(rval_rgset(), type = "II")[, "Name"])
    #
    beta_raw_green <- reactive(subset(rval_rgset_getBeta(), rownames(rval_rgset_getBeta()) %in% channel_green()))
    beta_raw_red <- reactive(subset(rval_rgset_getBeta(), rownames(rval_rgset_getBeta()) %in% channel_red()))
    beta_raw_II <- reactive(subset(rval_rgset_getBeta(), rownames(rval_rgset_getBeta()) %in% channel_II()))
    n_raw_green <- reactive(ifelse(nrow(beta_raw_green()) < 20000, nrow(beta_raw_green()), 20000))
    n_raw_red <- reactive(ifelse(nrow(beta_raw_red()) < 20000, nrow(beta_raw_red()), 20000))
    n_raw_II <- reactive(ifelse(nrow(beta_raw_II()) < 20000, nrow(beta_raw_II()), 20000))
    rval_plot_densityplotraw_green <- reactive(create_densityplot(beta_raw_green(), n_raw_green()))
    rval_plot_densityplotraw_red <- reactive(create_densityplot(beta_raw_red(), n_raw_red()))
    rval_plot_densityplotraw_II <- reactive(create_densityplot(beta_raw_II(), n_raw_II()))
    #
    beta_normalized_green <- reactive(rval_gset_getBeta()[rownames(rval_gset_getBeta()) %in% channel_green(),])
    beta_normalized_red <- reactive(rval_gset_getBeta()[rownames(rval_gset_getBeta()) %in% channel_red(),])
    beta_normalized_II <- reactive(rval_gset_getBeta()[rownames(rval_gset_getBeta()) %in% channel_II(),])
    n_normalized_green <- reactive(ifelse(nrow(beta_normalized_green()) < 20000, nrow(beta_normalized_green()), 20000))
    n_normalized_red <- reactive(ifelse(nrow(beta_normalized_red()) < 20000, nrow(beta_normalized_red()), 20000))
    n_normalized_II <- reactive(ifelse(nrow(beta_normalized_II()) < 20000, nrow(beta_normalized_II()), 20000))
    rval_plot_densityplot_green <- reactive(create_densityplot(beta_normalized_green(), n_normalized_green()))
    rval_plot_densityplot_red <- reactive(create_densityplot(beta_normalized_red(), n_normalized_red()))
    rval_plot_densityplot_II <- reactive(create_densityplot(beta_normalized_II(), n_normalized_II()))
    ###
    n_raw_all <- reactive(ifelse(nrow(rval_rgset_getBeta()) < 20000, nrow(rval_rgset_getBeta()), 20000))
    rval_plot_densityplotraw_all <- reactive(create_densityplot(rval_rgset_getBeta(), n_raw_all()))
    n_normalized_all <- reactive(ifelse(nrow(rval_gset_getBeta()) < 20000, nrow(rval_gset_getBeta()), 20000))
    rval_plot_densityplot_all <- reactive(create_densityplot(rval_gset_getBeta(), n_normalized()))
    ###
    
    
    beta_raw <- reactive(subset(rval_rgset_getBeta(), rownames(rval_rgset_getBeta()) %in% channel()))
    n_raw <- reactive(ifelse(nrow(beta_raw()) < 20000, nrow(beta_raw()), 20000))
    rval_plot_densityplotraw <- reactive(create_densityplot(beta_raw(), n_raw()))
    
    beta_normalized <- reactive(rval_gset_getBeta()[rownames(rval_gset_getBeta()) %in% channel(),])
    n_normalized <- reactive(ifelse(nrow(beta_normalized()) < 20000, nrow(beta_normalized()), 20000))
    rval_plot_densityplot <- reactive(create_densityplot(beta_normalized(), n_normalized()))
    
    # Density plots
    #rval_plot_densityplotraw <- reactive(create_densityplot(rval_rgset_getBeta(), 200000))
    #rval_plot_densityplot <- reactive(create_densityplot(rval_gset_getBeta(), 200000))
    
    output$graph_minfi_densityplotraw <- plotly::renderPlotly(rval_plot_densityplotraw())
    output$graph_minfi_densityplot <- plotly::renderPlotly(rval_plot_densityplot())
    ##### SNP HEATMAP #####
    rval_plot_snpheatmap <- reactive(
        create_snpheatmap(
            minfi::getSnpBeta(rval_rgset()),
            rval_sheet_target()[, input$select_input_samplenamevar],
            rval_sheet_target()[, input$select_input_donorvar]
        )
    )
    
    output$graph_minfi_snps <- plotly::renderPlotly(rval_plot_snpheatmap())
    
    ##### BATCH EFECTS #####
    rval_plot_corrplot <- reactive(
        create_corrplot(
            rval_gset_getBeta(),
            rval_clean_sheet_target(),
            rval_sheet_target(),
            p.value = input$select_minfi_typecorrplot == "p.value"
        )
    )
    
    output$graph_minfi_corrplot <- plotly::renderPlotly(rval_plot_corrplot()[["graph"]])
    output$table_minfi_corrplot <- DT::renderDT(rval_plot_corrplot()[["info"]],
                                                rownames = FALSE,
                                                selection = "single",
                                                style = "bootstrap",
                                                caption = "Autodetected variable types:",
                                                options = list(pageLength = 10, autoWidth = TRUE)
    )
    
    ########## VIOLIN PLOT #####################
    
    create_violinplot <- function(Bvalues, n = 200000) {
        # Creating density plot using a sample of n CpGs
        
        Bvalues[sample(seq_len(nrow(Bvalues)), n), ] %>%
            tidyr::pivot_longer(
                cols = seq_len(ncol(Bvalues)),
                names_to = "sample",
                values_to = "Bvalues"
            ) %>% 
            ggplot2::ggplot(ggplot2::aes(
                x = .data$Bvalues, y = .data$sample, color = .data$sample
            )) +
            ggplot2::geom_violin() +
            ggplot2::geom_vline(xintercept = 0.5, 
                                linetype = "dashed", 
                                alpha = 0.5) +
            ggplot2::stat_summary(fun = mean, 
                                  geom = "crossbar") +
            ggplot2::xlab("Beta") +
            ggplot2::ylab("") +
            ggplot2::theme_bw() +
            ggplot2::theme(legend.position = "none",
                           panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank())
    }
    
    rval_plot_violin_raw <- reactive(create_violinplot(rval_rgset_getBeta(), nrow(rval_rgset_getBeta())))
    output$graph_violin_raw <- renderPlot(rval_plot_violin_raw())
    
    rval_plot_violin_normalized <- reactive(create_violinplot(rval_gset_getBeta(), nrow(rval_gset_getBeta())))
    output$graph_violin_normalized <- renderPlot(rval_plot_violin_normalized())
    
    
    ########## PCA PLOT ##########
    
    rval_plot_pca <- eventReactive(
        list(input$button_pca_update, input$button_minfi_select),
        create_pca(
            Bvalues = rval_gset_getBeta(),
            pheno_info = as.data.frame(minfi::pData(rval_gset())),
            group = input$select_input_samplenamevar,
            pc_x = input$select_minfi_pcaplot_pcx,
            pc_y = input$select_minfi_pcaplot_pcy,
            color = input$select_minfi_pcaplot_color
        )
    )
    
    output$graph_minfi_pcaplot <- plotly::renderPlotly(rval_plot_pca()[["graph"]])
    
    output$table_minfi_pcaplot <- DT::renderDT(
        format(rval_plot_pca()[["info"]], digits = 2),
        rownames = TRUE,
        selection = "single",
        style = "bootstrap",
        options = list(
            autoWidth = TRUE,
            paging = FALSE,
            scrollX = TRUE,
            lengthChange = FALSE,
            searching = FALSE,
            info = FALSE
        )
    )
    
    ##### DECONVOLUTION #####
    deconvolution <- reactive(estimateCellCounts(rval_rgset()))
    
    graph_deconvolution <- reactive(pheatmap::pheatmap(deconvolution()))
    
    output$deconvolution_heatmap <- renderPlot(graph_deconvolution())
    
    
    ##### HEATMAP #####
    
    plot_random_heatmap <- reactive(create_random_heatmap(rval_rgset_getBeta()))
    plot_top_heatmap <- reactive(create_top_heatmap(rval_rgset_getBeta()))
    output$graph_random_heatmap <- renderPlot(plot_random_heatmap())
    output$graph_top_heatmap <- renderPlot(plot_top_heatmap())

    
    ##### AGE #####
    
    age_data <- reactive(create_age(rval_rgset_getBeta(), rval_sheet_target(), input$select_input_age))
    
    output$table_age <- DT::renderDT(age_data())
    
    
    ##### HYPER/HYPO PLOTS #####
    
    graph_hyper_hypo <- eventReactive(list(input$button_hyper_hypo_update, input$button_input_next), create_hyper_hypo(rval_rgset(), rval_rgset_getBeta(), input$slider_beta, input$selected_samples_h))
    output$plot_chr <- renderPlot(graph_hyper_hypo()[["chr"]])
    output$plot_relation_to_island <- renderPlot(graph_hyper_hypo()[["relation_to_island"]])
    output$plot_group <- renderPlot(graph_hyper_hypo()[["group"]])
    
    
    
    
    # Update of next form and move to Limma 
    observeEvent(input$button_minfi_select, {
        # check if normalization has been performed
        req(rval_gset())
        
        updatePickerInput(
            session,
            "select_limma_voi",
            label = "Select Variable of Interest",
            choices = colnames(rval_clean_sheet_target())[vapply(rval_clean_sheet_target(), is.factor, logical(1))],
            selected = input$select_input_groupingvar
        )
        
        updatePickerInput(
            session,
            "checkbox_limma_covariables",
            label = "Select linear model covariables",
            choices = colnames(rval_clean_sheet_target()),
            selected = input$select_input_donorvar
        )
        
        # possible interactions of variables:
        if (length(colnames(rval_clean_sheet_target()) > 2)) {
            interactions <- utils::combn(colnames(rval_clean_sheet_target()), 2)
            interactions <- sprintf("%s:%s", interactions[1, ], interactions[2, ])
        }
        else {
            interactions <- c()
        }
        
        
        updatePickerInput(
            session,
            "checkbox_limma_interactions",
            label = "Select linear model interactions",
            choices = interactions,
            selected = input$select_input_donorvar
        )
        
        shinyjs::enable("button_limma_calculatemodel")
    })
    
    
    ###### LIMMA #####
    
    # Variable of interest
    rval_voi <- reactive(factor(make.names(minfi::pData(rval_gset())[, input$select_limma_voi])))
    
    # Calculation of contrasts
    rval_contrasts <- reactive({
        generate_contrasts(rval_voi())
    })
    
    # Design calculation
    rval_fit <- eventReactive(input$button_limma_calculatemodel, {
        print("design")
        req(input$select_limma_voi) # a variable of interest is required
        print("design2")
        # Bulding the design matrix
        try({
            design <- generate_design(
                voi = input$select_limma_voi, sample_name = input$select_input_samplenamevar,
                covariables = input$checkbox_limma_covariables, interactions = input$checkbox_limma_interactions,
                sample_sheet = rval_clean_sheet_target()
            )
        })
        
        validate(
            need(
                exists("design", inherits = FALSE),
                "Design matrix has failed. Please, check your options and try again."
            ),
            need(
                nrow(design) == length(rval_sheet_target()[[input$select_input_samplenamevar]]),
                "The design matrix contains missing values. Please, check the selected variable and covariables."
            )
        )
        
        # checking colinearity
        if (qr(design)$rank < ncol(design)) {
            showModal(
                modalDialog(
                    title = "Colinearity warning",
                    "The design matrix presents colinear columns. Even though it is possible to try to continue the analysis, checking the selected covariables is recommended.",
                    easyClose = TRUE,
                    footer = NULL
                )
            )
        }
        print(design)
        
        
        print("fit")
        
        validate(
            need(input$input_data != "", "DMP calculation has not been performed or data has not been uploaded.")
        )
        print("ok")
        req(design)
        
        print("fit2")
        
        shinyjs::disable("button_limma_calculatemodel") # disable button to avoid repeat clicking

        showModal(modalDialog(
            title = NULL, footer = NULL,
            div(
                img(src="https://upload.wikimedia.org/wikipedia/commons/7/7d/Pedro_luis_romani_ruiz.gif"),
                p("Generating linear model..."),
                style = "margin: auto; text-align: center"
            )
        ))
        
        withProgress(
            message = "Generating linear model...",
            value = 3,
            max = 6,
            {
                try({
                    fit <- generate_limma_fit(
                        Mvalues = rval_gset_getM(), design = design,
                        weighting = as.logical(input$select_limma_weights)
                    )
                })
                
                if (!exists("fit", inherits = FALSE)) {
                    rval_generated_limma_model(FALSE) # disable contrast button
                    rval_analysis_finished(FALSE) # disable download buttons
                    
                    removeModal()
                    showModal(
                        modalDialog(
                            title = "lmFit error",
                            "lmFit model has failed. Please, check your options and try again.",
                            easyClose = TRUE,
                            footer = NULL
                        )
                    )
                }
            }
        )
        
        
        shinyjs::enable("button_limma_calculatemodel") # enable button again to allow repeting calculation
        
        removeModal()
        
        validate(need(
            exists("fit", inherits = FALSE),
            "lmFit model has failed. Please, check your options and try again."
        ))
        print(fit)
        
        fit
    })

    # rval_fit() has NAs, we remove the option to trend or robust in eBayes to prevent failure
    observeEvent(input$button_limma_calculatemodel, {
        print("next")
        if (any(vapply(rval_fit(), function(x) {
            any(is.na(unlist(x)) |
                unlist(x) == Inf | unlist(x) == -Inf)
        }, logical(1)))) {
            message("NAs or Inf values detected, trend and robust options are disabled.")
            
            updateSwitchInput(session,
                              "select_limma_trend",
                              value = FALSE,
                              disabled = TRUE
            )
            
            updateSwitchInput(session,
                              "select_limma_robust",
                              value = FALSE,
                              disabled = TRUE
            )
        }
        
        else {
            message("NAs or Inf values not detected, trend and robust options are enabled")
            
            updateSwitchInput(session,
                              "select_limma_trend",
                              value = FALSE,
                              disabled = FALSE
            )
            
            updateSwitchInput(session,
                              "select_limma_robust",
                              value = FALSE,
                              disabled = FALSE
            )
        }
        print("final")
        # Creating calculate differences button
        rval_generated_limma_model(TRUE)
    })
    
    output$button_limma_calculatedifs_container <- renderUI({
        if (rval_generated_limma_model()) {
            return(tagList(
                br(),
                
                h4("Contrasts options"),
                
                switchInput(
                    inputId = "select_limma_trend",
                    label = "eBayes Trend",
                    labelWidth = "80px",
                    value = FALSE
                ),
                
                switchInput(
                    inputId = "select_limma_robust",
                    label = "eBayes Robust",
                    labelWidth = "80px",
                    value = FALSE
                ),
                
                actionButton("button_limma_calculatedifs", "Calc. Contrasts")
            ))
        } else {
            return()
        }
    })
    
    
    # render of plots and tables
    
    rval_plot_plotSA <- reactive(create_plotSA(rval_fit()))
    output$graph_limma_plotSA <- renderPlot(rval_plot_plotSA())
    output$table_limma_design <- DT::renderDT(
        rval_fit()$design,
        rownames = TRUE,
        selection = "single",
        style = "bootstrap",
        options = list(
            autoWidth = TRUE,
            scrollX = TRUE,
            lengthChange = FALSE,
            searching = FALSE
        )
    )
    
    
    # Calculation of global difs
    rval_globaldifs <- eventReactive(input$button_limma_calculatedifs, {
        print("enter globaldifs")
        try({
            globaldifs <- calculate_global_difs(rval_gset_getBeta(), rval_voi(), rval_contrasts(),
                                                cores = n_cores
            )
        })
        
        if (!exists("globaldifs", inherits = FALSE)) {
            showModal(
                modalDialog(
                    title = "Global difs calculation error",
                    "An unexpected error has ocurred during global diffs and means calculation. Please, check your selected samples.",
                    easyClose = TRUE,
                    footer = NULL
                )
            )
        }
        
        validate(
            need(
                exists("globaldifs", inherits = FALSE),
                "An unexpected error has ocurred during global diffs and means calculation."
            )
        )
        print("globaldifs")
        print(globaldifs)
        globaldifs
    })
    
    # Calculation of differences (eBayes)
    rval_finddifcpgs <- eventReactive(input$button_limma_calculatedifs, {
        try({
            dif_cpgs <- find_dif_cpgs(
                design = rval_fit()$design,
                fit = rval_fit(),
                contrasts = rval_contrasts(),
                trend = as.logical(input$select_limma_trend),
                robust = as.logical(input$select_limma_robust),
                cores = n_cores
            )
        })
        
        if (!exists("dif_cpgs", inherits = FALSE)) {
            showModal(
                modalDialog(
                    title = "Contrasts Calculation Error",
                    "An unexpected error has ocurred during contrasts calculation. Please, generate the model again and check if it is correct.
        If the problem persists, report the error to the maintainer",
                    easyn_coresose = TRUE,
                    footer = NULL
                )
            )
            
            rval_generated_limma_model(FALSE)
            rval_analysis_finished(FALSE)
            shinyjs::disable("button_limma_heatmapcalc")
        }
        
        validate(
            need(
                exists("dif_cpgs", inherits = FALSE),
                "An unexpected error has ocurred during contrasts calculation. Please, generate the model again and check if it is correct.
        If the problem persists, report the error to the maintainer"
            )
        )
        
        rval_analysis_finished(TRUE)
        print("dif_cpgs")
        print(dif_cpgs)
        dif_cpgs
    })
    
    # Update of heatmap controls
    observeEvent(input$button_limma_calculatedifs, {
        updateTabItems(session, "tabset_limma", "differential_cpgs")
        
        updateSelectInput(
            session,
            "select_limma_groups2plot",
            label = "Groups to plot",
            choices = levels(rval_voi()),
            selected = levels(rval_voi())
        )
        
        updateSelectInput(
            session,
            "select_limma_contrasts2plot",
            label = "Contrasts to plot",
            choices = rval_contrasts(),
            selected = rval_contrasts()
        )
        
        updateSelectInput(
            session,
            "select_limma_anncontrast",
            label = "Contrast",
            choices = rval_contrasts()
        )
        
        updateSelectInput(
            session,
            "select_anncontrast_manhattan",
            label = "Contrast",
            choices = rval_contrasts()
        )
        
        updateSelectInput(
            session,
            "select_anncontrast_volcano",
            label = "Contrast",
            choices = rval_contrasts()
        )
        
        # disable button to avoid repeat n_coresicking
        shinyjs::disable("button_limma_calculatedifs")
        
        # force rval_filteredlist
        rval_filteredlist()
        
        # enable or disable removebatch option
        covariables_design <- as.matrix(rval_fit()$design[, -seq_len(length(unique(rval_voi())))])
        
        if (ncol(covariables_design) > 0) {
            updateSwitchInput(session,
                              "select_limma_removebatch",
                              value = FALSE,
                              disabled = FALSE
            )
            updateSwitchInput(session,
                              "select_dmrs_removebatch",
                              value = FALSE,
                              disabled = FALSE
            )
        }
        else {
            updateSwitchInput(session,
                              "select_limma_removebatch",
                              value = FALSE,
                              disabled = TRUE
            )
            updateSwitchInput(session,
                              "select_dmrs_removebatch",
                              value = FALSE,
                              disabled = FALSE
            )
        }
        
        # enable and n_coresick heatmap button to get default graph
        shinyjs::enable("button_limma_heatmapcalc")
        #shinyjs::n_coresick("button_limma_heatmapcalc")
        
        # enable again the button to allow repeat calculation
        shinyjs::enable("button_limma_calculatedifs")
    })
    
    # Calculation of filtered list
    rval_filteredlist <- eventReactive(rval_finddifcpgs(), {

        showModal(modalDialog(
            title = NULL, footer = NULL,
            div(
                img(src="https://upload.wikimedia.org/wikipedia/commons/7/7d/Pedro_luis_romani_ruiz.gif"),
                p("Performing contrasts calculations..."),
                style = "margin: auto; text-align: center"
            )
        ))
        
        withProgress(
            message = "Performing contrasts calculations...",
            value = 1,
            max = 6,
            {
                setProgress(message = "Calculating global difs...", value = 1)
                rval_globaldifs()
                setProgress(message = "Calculating eBayes...", value = 4)
                rval_finddifcpgs()
                setProgress(message = "Calculating filtered list...", value = 5)
                removeModal()
                create_filtered_list(
                    rval_finddifcpgs(),
                    rval_globaldifs(),
                    deltaB = input$slider_limma_deltab,
                    adjp_max = input$slider_limma_adjpvalue,
                    p.value = input$slider_limma_pvalue,
                    cores = n_cores
                )
               
            }
        )
        
    })
    
    rval_list <- reactive({
        showModal(modalDialog(
            title = NULL, footer = NULL,
            div(
                p("content"),
                img(src="https://upload.wikimedia.org/wikipedia/commons/7/7d/Pedro_luis_romani_ruiz.gif"),
                style = "margin: auto; text-align: center"
            )
        ))
        withProgress(
            message = "Performing contrasts calculations...",
            value = 1,
            max = 6,
            {
                setProgress(message = "Calculating global difs...", value = 1)
                rval_globaldifs()
                setProgress(message = "Calculating eBayes...", value = 4)
                rval_finddifcpgs()
                setProgress(message = "Calculating filtered list...", value = 5)
                removeModal()
                create_list(
                    rval_finddifcpgs(),
                    rval_globaldifs(),
                    cores = n_cores
                )
                
            }
        )
        
    })
    
    rval_filteredlist2heatmap <- reactive({
        join_table <- create_dmps_heatdata(
            rval_filteredlist(),
            input$select_limma_contrasts2plot,
            input$select_limma_removebatch,
            rval_fit()$design,
            rval_voi(),
            rval_gset_getBeta()
        )
        print("join table")
        print(join_table)
        # If the number of CpGs is not in the plotting range, return NULL to avoid errors in plot_heatmap and disable download
        if (is.null(join_table) |
            nrow(join_table) < 2 | nrow(join_table) > 12000) {
            rval_filteredlist2heatmap_valid(FALSE)
            NULL
        }
        else {
            rval_filteredlist2heatmap_valid(TRUE)
            join_table
        }
    })
    
    rval_cpgcount_heatmap <- eventReactive(list(input$button_limma_calculatedifs, input$button_limma_heatmapcalc), nrow(rval_filteredlist2heatmap()))
    
    rval_dendrogram <- eventReactive(list(input$button_limma_calculatedifs, input$button_limma_heatmapcalc), {
        if (input$select_limma_rowsidecolors) {
            print("if dendogram")
            # check if dendrogram cutting works (k should be minor than heatmap rows)
            try({
                dendrogram <- create_dendrogram(
                    rval_filteredlist2heatmap(),
                    factorgroups = factor(rval_voi()[rval_voi() %in% input$select_limma_groups2plot],
                                          levels = input$select_limma_groups2plot
                    ),
                    groups2plot = rval_voi() %in% input$select_limma_groups2plot,
                    clusteralg = input$select_limma_clusteralg,
                    distance = input$select_limma_clusterdist,
                    scale_selection = input$select_limma_scale,
                    k_number = input$select_limma_knumber
                )
            })
            print("final if dendogram")
        } else {
            print("else dendogram")
            dendrogram <- NULL
        }
        
        
        if (!exists("dendrogram", inherits = FALSE)) {
            print("if exists")
            dendrogram <- NULL
            showModal(
                modalDialog(
                    title = "Row clustering error",
                    "An error has ocurred during cluster cutting from row dendrogram. Maybe the number of clusters selected is too high.",
                    easyClose = TRUE,
                    footer = NULL
                )
            )
            print("final if exists")
        }
        print("dendogram")
        dendrogram # returning the dendrogram classification
    })
    
    plot_heatmap <- eventReactive(list(input$button_limma_calculatedifs, input$button_limma_heatmapcalc), ignoreNULL = FALSE, {
        validate(
            need(
               !is.null(rval_filteredlist2heatmap()),
                "Differences are not in the plotting range (<12000, >1)"
            ),
            need(
                !is.null(input$select_limma_groups2plot) &
                    input$select_limma_groups2plot != "",
                "Select at least one group to plot."
            )
        )
        
        create_heatmap(
            rval_filteredlist2heatmap(),
            factorgroups = factor(rval_voi()[rval_voi() %in% input$select_limma_groups2plot],
                                  levels = input$select_limma_groups2plot
            ),
            groups2plot = rval_voi() %in% input$select_limma_groups2plot,
            Colv = as.logical(input$select_limma_colv),
            ColSideColors = input$select_limma_colsidecolors,
            RowSideColors = rval_dendrogram(),
            clusteralg = input$select_limma_clusteralg,
            distance = input$select_limma_clusterdist,
            scale = input$select_limma_scale,
            static = as.logical(input$select_limma_graphstatic)
        )
    })
    
    make_table <- eventReactive(list(input$button_limma_calculatedifs, input$button_limma_tablecalc), {
        default_df <- data.frame(
            contrast = rval_contrasts(),
            Hypermethylated = 0,
            Hypomethylated = 0,
            total = 0
        )
        
        count_df <- data.table::rbindlist(rval_filteredlist(), idcol = "contrast")
        
        if (nrow(count_df) > 0 & !is.null(count_df)) {
            count_df <- count_df %>%
                dplyr::mutate(type = factor(
                    ifelse(
                        .data$dif_current < 0,
                        "Hypermethylated",
                        "Hypomethylated"
                    ),
                    levels = c("Hypermethylated", "Hypomethylated")
                )) %>%
                dplyr::group_by(.data$contrast, .data$type) %>%
                dplyr::summarise(CpGs = dplyr::n()) %>%
                tidyr::complete(.data$contrast, .data$type, fill = list(
                    CpGs =
                        0
                )) %>%
                tidyr::pivot_wider(
                    names_from = .data$type,
                    values_from = .data$CpGs
                ) %>%
                dplyr::mutate(total = .data$Hypermethylated + .data$Hypomethylated) %>%
                dplyr::mutate(dplyr::across(
                    c("Hypermethylated", "Hypomethylated", "total"),
                    ~ format(., scientific = FALSE, digits = 0)
                ))
        }
        
        rbind(as.data.frame(count_df), default_df[!(default_df[["contrast"]] %in% count_df[["contrast"]]), ])
    })
    
    observe({
        # Render the correct plot depending on the selected
        output$graph_limma_heatmapcontainer <- renderUI({
            if (!as.logical(input$select_limma_graphstatic)) {
                return(
                    plotly::plotlyOutput(
                        "graph_limma_heatmap_interactive",
                        width = "600px",
                        height = "500px"
                    ) %>% shinycssloaders::withSpinner()
                )
            } else {
                return(
                    div(
                        min_width = 400,
                        plotOutput(
                            "graph_limma_heatmap_static",
                            width = "600px",
                            height = "600px"
                        ) %>% shinycssloaders::withSpinner()
                    )
                )
            }
        })
        
        if (!as.logical(input$select_limma_graphstatic)) {
            output$graph_limma_heatmap_interactive <- plotly::renderPlotly(plot_heatmap())
        } else {
            output$graph_limma_heatmap_static <- renderPlot(plot_heatmap())
        }
    })
    
    output$text_limma_heatmapcount <- renderText(paste("DMPs in heatmap:", rval_cpgcount_heatmap()))
    
    output$table_limma_difcpgs <- renderTable(make_table(), digits = 0)
    
    table_annotation <- eventReactive(list(input$button_limma_calculatedifs, input$button_limma_heatmapcalc, input$select_limma_anncontrast), {
        req(rval_filteredlist())
        print("TABLE ANNOTATION")
        dif_target <- paste("dif",
                            limma::strsplit2(input$select_limma_anncontrast, "-")[1],
                            limma::strsplit2(input$select_limma_anncontrast, "-")[2],
                            sep = "_"
        )
        
        temp <- rval_annotation()[row.names(rval_annotation()) %in% rval_filteredlist()[[input$select_limma_anncontrast]]$cpg, ]

        temp$dif_beta <- rval_globaldifs()[[dif_target]][rval_globaldifs()[["cpg"]] %in% row.names(temp)]

        temp$fdr <- rval_filteredlist()[[input$select_limma_anncontrast]][["adj.P.Val"]][rval_filteredlist()[[input$select_limma_anncontrast]][["cpg"]] %in% row.names(temp)]

        temp$pvalue <- rval_filteredlist()[[input$select_limma_anncontrast]][["P.Value"]][rval_filteredlist()[[input$select_limma_anncontrast]][["cpg"]] %in% row.names(temp)]

        temp$chr <- as.numeric(as.character(gsub("chr", "", temp$chr)))
        gene <- vapply(strsplit(temp$UCSC_RefGene_Name,";"), `[`, 1, FUN.VALUE=character(1))
        gene[is.na(gene)]<-""
        temp$gene <- gene

        print("temp")
        
        temp
    })
    
    output$table_limma_ann <- DT::renderDT(
        table_annotation(),
        extensions = "Buttons",
        rownames = FALSE,
        selection = "single",
        style = "bootstrap",
        caption = "DMPs Annotation:",
        options = list(
            dom = "Blfrtip",
            lengthMenu = list(c(10, 25, 50, 100, 25000), c(10, 25, 50, 100, 25000)),
            pageLength = 10,
            autoWidth = TRUE,
            scrollX = TRUE,
            buttons = c("csv", "excel", "print")
        )
    )
    
    ind_boxplot <- eventReactive(input$button_limma_indboxplotcalc, {
        validate(need(!is.null(input$table_limma_ann_rows_selected), "A DMP should be selected."))
        cpg_sel <- table_annotation()[["Name"]][input$table_limma_ann_rows_selected]
        
        create_individual_boxplot(rval_gset_getBeta(), cpg_sel, rval_voi())
    })
    
    output$graph_limma_indboxplot <- renderPlot(ind_boxplot())
    
    
    ########## MANHATTAN & VOLCANO PLOT ##########
    
    table_annotation_manhattan <- eventReactive(input$select_anncontrast_manhattan, {
        req(rval_list())
        
        dif_target <- paste("dif",
                            limma::strsplit2(input$select_anncontrast_manhattan, "-")[1],
                            limma::strsplit2(input$select_anncontrast_manhattan, "-")[2],
                            sep = "_"
        )
        
        temp <- rval_annotation()[row.names(rval_annotation()) %in% rval_list()[[input$select_anncontrast_manhattan]]$cpg, ]
        temp$dif_beta <- rval_globaldifs()[[dif_target]][rval_globaldifs()[["cpg"]] %in% row.names(temp)]
        
        temp$fdr <- rval_list()[[input$select_anncontrast_manhattan]][["adj.P.Val"]][rval_list()[[input$select_anncontrast_manhattan]][["cpg"]] %in% row.names(temp)]
        
        temp$pvalue <- rval_list()[[input$select_anncontrast_manhattan]][["P.Value"]][rval_list()[[input$select_anncontrast_manhattan]][["cpg"]] %in% row.names(temp)]
        
        temp$chr[temp$chr == "chrX"] <- 23
        temp$chr[temp$chr == "chrY"] <- 24
        #temp$chr[temp$chr == "chrM"] <- 25
        
        
        
        print(unique(temp$chr))
        
        
        temp$chr <- as.numeric(as.character(gsub("chr", "", temp$chr)))
        gene <- vapply(strsplit(temp$UCSC_RefGene_Name,";"), `[`, 1, FUN.VALUE=character(1))
        gene[is.na(gene)]<-""
        temp$gene <- gene
        
        temp
    })
    
    table_annotation_volcano <- eventReactive(input$select_anncontrast_volcano, {
        req(rval_list())
        
        dif_target <- paste("dif",
                            limma::strsplit2(input$select_anncontrast_volcano, "-")[1],
                            limma::strsplit2(input$select_anncontrast_volcano, "-")[2],
                            sep = "_"
        )
        
        temp <- rval_annotation()[row.names(rval_annotation()) %in% rval_list()[[input$select_anncontrast_volcano]]$cpg, ]
        temp$dif_beta <- rval_globaldifs()[[dif_target]][rval_globaldifs()[["cpg"]] %in% row.names(temp)]
        
        temp$fdr <- rval_list()[[input$select_anncontrast_volcano]][["adj.P.Val"]][rval_list()[[input$select_anncontrast_volcano]][["cpg"]] %in% row.names(temp)]
        
        temp$pvalue <- rval_list()[[input$select_anncontrast_volcano]][["P.Value"]][rval_list()[[input$select_anncontrast_volcano]][["cpg"]] %in% row.names(temp)]
        
        temp$chr[temp$chr == "chrX"] <- 23
        temp$chr[temp$chr == "chrY"] <- 24
        #temp$chr[temp$chr == "chrM"] <- 25
        
        
        
        print(unique(temp$chr))
        
        
        temp$chr <- as.numeric(as.character(gsub("chr", "", temp$chr)))
        gene <- vapply(strsplit(temp$UCSC_RefGene_Name,";"), `[`, 1, FUN.VALUE=character(1))
        gene[is.na(gene)]<-""
        temp$gene <- gene
        
        temp
    })
    
    
    volcano_data <- eventReactive(table_annotation_volcano(), {
        pval <- table_annotation_volcano()$pvalue
        fc <- table_annotation_volcano()$dif_beta
        names <- table_annotation_volcano()$gene
        tFC <- 0.2
        show.labels <- T
        dta <- data.frame(P.Value = pval, FC = fc, names, clr = "gray87", alp = 0.5, stringsAsFactors = FALSE)
        head(dta)
        dta$PV <- -log10(dta$P.Value)
        dta$feature <- rownames(dta)
        dta$clr[abs(dta$FC) >= tFC] <- "olivedrab"
        dta$alp[abs(dta$FC) >= tFC] <- 0.7
        tPV <- -log10(0.001)
        dta$clr[dta$PV >= tPV] <- "tan3"
        dta$alp[dta$PV >= tPV] <- 0.7
        dta$clr[dta$PV >= tPV & abs(dta$FC) >= tFC] <- "lightskyblue"
        dta$alp[dta$PV >= tPV & abs(dta$FC) >= tFC] <- 0.9
        clrvalues <- c("gray87", "tan3", "olivedrab", "lightskyblue")
        names(clrvalues) <- c("gray87", "tan3", "olivedrab", "lightskyblue")
        print(head(dta))
        dta
    })
    
    manhattan_graph <- reactive(qqman::manhattan(table_annotation_manhattan(), chr = "chr", bp = "pos", snp = "gene", p = "pvalue",
                                                 annotatePval = 1, suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08), annotateTop = T))
    volcano_graph <- reactive(MultiDataSet::volcano_plot(pval = table_annotation_volcano()$pvalue, fc = table_annotation_volcano()$dif_beta,
                                                         table_annotation_volcano()$gene, tFC = 0.2, show.labels = T))
    
    
    output$manhattan_plot <- renderPlot(manhattan_graph())
    output$volcano_plot <- renderPlot(volcano_graph())
    
    

    
    ###### DMRs #####
    
    rval_mcsea <- eventReactive(input$button_dmrs_calculate, {
        validate(
            need(
                rval_analysis_finished(),
                "To calculate DMRs, you have to finish first the DMP calculation."
            )
        )
        
        validate(
            need(
                requireNamespace("mCSEA", quietly = TRUE),
                "mCSEA is not installed. You should install the package to calculate DMRs"
            ),
            need(
                input$select_dmrs_contrasts != "",
                "You should select at least one contrast."
            ),
            need(
                input$select_dmrs_regions != "",
                "You should select at least one DMR type."
            )
        )
        
        updateSelectInput(session,
                          "select_dmrs_selcont",
                          choices = input$select_dmrs_contrasts
        )
        updateSelectInput(session,
                          "select_dmrs_selreg",
                          choices = input$select_dmrs_regions
        )
        
        updateSelectInput(
            session,
            "select_dmrs_groups2plot",
            label = "Groups to plot",
            choices = levels(rval_voi()),
            selected = levels(rval_voi())
        )
        updateSelectInput(
            session,
            "select_dmrs_contrasts2plot",
            label = "Contrasts to plot",
            choices = input$select_dmrs_contrasts,
            selected = input$select_dmrs_contrasts
        )
        updateSelectInput(
            session,
            "select_dmrs_regions2plot",
            label = "Regions to plot",
            choices = input$select_dmrs_regions,
            selected = input$select_dmrs_regions
        )
        
        shinyjs::disable("button_dmrs_calculate")

        showModal(modalDialog(
            title = NULL, footer = NULL,
            div(
                img(src="https://upload.wikimedia.org/wikipedia/commons/7/7d/Pedro_luis_romani_ruiz.gif"),
                p("Calculating DMRs..."),
                style = "margin: auto; text-align: center"
            )
        ))
        
        withProgress(
            message = "Calculating DMRs...",
            max = 3,
            value = 1,
            {
                try({
                    dmrs_result <- find_dmrs(
                        rval_finddifcpgs(),
                        minCpGs = input$slider_dmrs_cpgs,
                        platform = "EPIC",
                        voi = rval_voi(),
                        regionsTypes = input$select_dmrs_regions,
                        contrasts = input$select_dmrs_contrasts,
                        bvalues = rval_gset_getBeta(),
                        permutations = input$slider_dmrs_permutations,
                        ncores = n_cores
                    )
                    
                    setProgress(value = 2, message = "Calculating differential of betas...")
                    
                    dmrs_result <- add_dmrs_globaldifs(
                        mcsea_result = dmrs_result,
                        cpg_globaldifs = rval_globaldifs(),
                        regionsTypes = input$select_dmrs_regions
                    )
                })
            }
        )
        
        shinyjs::enable("button_dmrs_calculate")
        
        removeModal()
        
        if (!exists("dmrs_result", inherits = FALSE)) {
            rval_dmrs_finished(FALSE)
            showModal(
                modalDialog(
                    title = "DMR calculation error",
                    "An unexpected error has occurred during DMRs calculation.",
                    easyClose = TRUE,
                    footer = NULL
                )
            )
            shinyjs::disable("button_dmrs_heatmapcalc")
        }
        
        validate(need(
            exists("dmrs_result", inherits = FALSE),
            "An unexpected error has occurred during DMRs calculation."
        ))
        
        # enable heatmap button
        shinyjs::enable("button_dmrs_heatmapcalc")
        rval_dmrs_finished(TRUE)
        rval_dmrs_ready2heatmap(TRUE)
        
        dmrs_result
    })
    
    rval_filteredmcsea <- eventReactive(list(input$button_dmrs_calculate, input$button_dmrs_tablecalc), {
        req(rval_mcsea())
        
        filter_dmrs(
            mcsea_list = rval_mcsea(),
            fdr = input$slider_dmrs_adjpvalue,
            pval = input$slider_dmrs_pvalue,
            dif_beta = input$slider_dmrs_deltab,
            regionsTypes = input$select_dmrs_regions,
            contrasts = input$select_dmrs_contrasts
        )
    })
    
    rval_filteredmcsea2heatmap <- reactive({
        req(input$select_dmrs_regions2plot)
        req(input$select_dmrs_contrasts2plot)
        
        # filtering contrasts
        join_table <- create_dmrs_heatdata(
            mcsea_result = rval_filteredmcsea(),
            bvalues = rval_gset_getBeta(),
            regions = input$select_dmrs_regions2plot,
            contrasts = input$select_dmrs_contrasts2plot,
            removebatch = input$select_dmrs_removebatch,
            design = rval_fit()$design,
            voi = rval_voi()
        )
        
        # If the number of CpGs is not in the plotting range, return NULL to avoid errors in plot_dmrsheatmap
        if (is.null(join_table) |
            nrow(join_table) < 2 | nrow(join_table) > 12000) {
            rval_filteredmcsea2heatmap_valid(FALSE)
            NULL
        }
        else {
            rval_filteredmcsea2heatmap_valid(TRUE)
            join_table
        }
    })
    
    observeEvent(rval_dmrs_ready2heatmap(), {
        if (rval_dmrs_ready2heatmap()) {
            rval_dmrs_ready2heatmap(FALSE)
            
            shinyjs::click("button_dmrs_heatmapcalc")
        }
    })
    
    # Heatmap DMRs
    
    plot_dmrsheatmap <- eventReactive(input$button_dmrs_heatmapcalc, ignoreNULL = FALSE, {
        validate(
            need(
                requireNamespace("mCSEA", quietly = TRUE),
                "mCSEA is not installed. You should install the package to calculate DMRs."
            )
        )
        
        validate(
            need(
                input$input_data != "",
                "DMR calculation has not been performed or data has not been uploaded."
            ),
            need(
                !is.null(input$select_dmrs_groups2plot) &
                    input$select_dmrs_groups2plot != "",
                "Select at least one group to plot."
            ),
            need(
                !is.null(rval_filteredmcsea2heatmap()),
                "Differences are not in the plotting range (<12000, >1)"
            )
        )
        
        create_heatmap(
            rval_filteredmcsea2heatmap(),
            factorgroups = factor(rval_voi()[rval_voi() %in% input$select_dmrs_groups2plot],
                                  levels = input$select_dmrs_groups2plot
            ),
            groups2plot = rval_voi() %in% input$select_dmrs_groups2plot,
            Colv = as.logical(input$select_dmrs_colv),
            ColSideColors = input$select_dmrs_colsidecolors,
            RowSideColors = rval_dendrogram_dmrs(),
            clusteralg = input$select_dmrs_clusteralg,
            distance = input$select_dmrs_clusterdist,
            scale = input$select_dmrs_scale,
            static = as.logical(input$select_dmrs_graphstatic)
        )
    })
    
    rval_dendrogram_dmrs <- eventReactive(input$button_dmrs_heatmapcalc, ignoreNULL = FALSE, {
        if (input$select_dmrs_rowsidecolors) {
            
            # check if dendrogram cutting works (k should be minor than heatmap rows)
            try({
                dendrogram <- create_dendrogram(
                    rval_filteredmcsea2heatmap(),
                    factorgroups = factor(rval_voi()[rval_voi() %in% input$select_dmrs_groups2plot],
                                          levels = input$select_dmrs_groups2plot
                    ),
                    groups2plot = rval_voi() %in% input$select_dmrs_groups2plot,
                    clusteralg = input$select_dmrs_clusteralg,
                    distance = input$select_dmrs_clusterdist,
                    scale_selection = input$select_dmrs_scale,
                    k_number = input$select_dmrs_knumber
                )
            })
        } else {
            dendrogram <- NULL
        }
        
        
        if (!exists("dendrogram", inherits = FALSE)) {
            dendrogram <- NULL
            showModal(
                modalDialog(
                    title = "Row clustering error",
                    "An error has ocurred during cluster cutting from row dendrogram. Maybe the number of clusters selected is too high.",
                    easyClose = TRUE,
                    footer = NULL
                )
            )
        }
        
        dendrogram # returning the dendrogram classification
    })
    rval_cpgcount_dmrs_heatmap <- eventReactive(input$button_dmrs_heatmapcalc, nrow(rval_filteredmcsea2heatmap()))
    
    observe({
        # Render the correct plot depending on the selected
        output$graph_dmrs_heatmapcontainer <- renderUI({
            if (!as.logical(input$select_dmrs_graphstatic)) {
                return(
                    plotly::plotlyOutput(
                        "graph_dmrs_heatmap_interactive",
                        width = "600px",
                        height = "500px"
                    ) %>% shinycssloaders::withSpinner()
                )
            } else {
                return(
                    div(
                        min_width = 400,
                        plotOutput(
                            "graph_dmrs_heatmap_static",
                            width = "600px",
                            height = "600px"
                        ) %>% shinycssloaders::withSpinner()
                    )
                )
            }
        })
        
        if (!as.logical(input$select_dmrs_graphstatic)) {
            output$graph_dmrs_heatmap_interactive <- plotly::renderPlotly(plot_dmrsheatmap())
        } else {
            output$graph_dmrs_heatmap_static <- renderPlot(plot_dmrsheatmap())
        }
    })
    
    output$text_dmrs_heatmapcount <- renderText(paste("DMRs in heatmap:", rval_cpgcount_dmrs_heatmap()))
    
    
    # DMRs count table
    
    make_table_dmrscount <- eventReactive(rval_filteredmcsea(), {
        result <- data.frame(
            contrast = rep(
                names(rval_filteredmcsea()),
                each = length(input$select_dmrs_regions)
            ),
            region_type = input$select_dmrs_regions
        )
        
        
        result$Hypermethylated <- apply(result, 1, function(x) {
            nrow(rval_filteredmcsea()[[x[1]]][[x[2]]][rval_filteredmcsea()[[x[1]]][[x[2]]]$NES < 0, ])
        })
        
        result$Hypomethylated <- apply(result, 1, function(x) {
            nrow(rval_filteredmcsea()[[x[1]]][[x[2]]][rval_filteredmcsea()[[x[1]]][[x[2]]]$NES > 0, ])
        })
        
        result$total <- result$Hypermethylated + result$Hypomethylated
        
        result
    })
    
    output$table_dmrs_count <- renderTable(make_table_dmrscount(), digits = 0)
    
    # DMRs single plot
    
    
    plot_singledmr <- eventReactive(input$button_dmrs_graphsingle, {
        validate(need(!is.null(input$table_dmrs_table_rows_selected), "A DMR should be selected."))
        
        selected_dmr <- row.names(rval_filteredmcsea()[[input$select_dmrs_selcont]][[input$select_dmrs_selreg]])[input$table_dmrs_table_rows_selected]
        col <- grDevices::rainbow(length(levels(rval_voi())))
        
        mCSEA::mCSEAPlot(
            rval_mcsea()[[input$select_dmrs_selcont]],
            regionType = input$select_dmrs_selreg,
            dmrName = selected_dmr,
            makePDF = FALSE,
            transcriptAnnotation = "symbol",
            col = col
        )
    })
    
    
    
    output$graph_dmrs_singledmr <- renderPlot(plot_singledmr())
    
    
    rval_table_sigdmrs <- reactive({
        rval_filteredmcsea()[[input$select_dmrs_selcont]][[input$select_dmrs_selreg]]
    })
    
    output$table_dmrs_table <- DT::renderDT(
        rval_table_sigdmrs(),
        rownames = TRUE,
        extensions = "Buttons",
        selection = "single",
        style = "bootstrap",
        caption = "Select DMR to plot:",
        options = list(
            pageLength = 10,
            dom = "Blfrtip",
            lengthMenu = list(c(10, 25, 50, 100, 25000), c(10, 25, 50, 100, 25000)),
            autoWidth = TRUE,
            scrollX = TRUE,
            buttons = c("csv", "excel", "print"),
            columnDefs = list(
                list(
                    targets = match("leadingEdge", colnames(rval_filteredmcsea()[[input$select_dmrs_selcont]][[input$select_dmrs_selreg]])),
                    visible = FALSE
                )
            )
        )
    )
    
    
    
    
    
    
    
    # Disable or enable buttons depending on software state
    observeEvent(
        rval_analysis_finished(),
        if (rval_analysis_finished()) {
            shinyjs::enable("download_export_robjects")
            shinyjs::enable("download_export_filteredbeds")
            shinyjs::enable("download_export_markdown")
            shinyjs::enable("download_export_script")
            shinyjs::enable("button_dmrs_calculate")
            
            updatePickerInput(
                session,
                "select_dmrs_contrasts",
                selected = rval_contrasts(),
                choices = rval_contrasts()
            )
            updatePickerInput(
                session,
                "select_dmrs_platform",
                selected = if (nrow(rval_finddifcpgs()[[1]]) > 500000) {
                    "EPIC"
                } else {
                    "450k"
                },
                choices = c("450k", "EPIC")
            )
        }
        
        else {
            shinyjs::disable("download_export_robjects")
            shinyjs::disable("download_export_filteredbeds")
            shinyjs::disable("download_export_markdown")
            shinyjs::disable("download_export_script")
            shinyjs::disable("download_export_heatmaps")
            shinyjs::disable("button_dmrs_calculate")
        }
    )
    
    
    
    
    
    
    
    
    
    
    
    rval_annotation <- eventReactive(rval_gset(), {
        int_cols <- c(
            "Name",
            "Relation_to_Island",
            "UCSC_RefGene_Name",
            "UCSC_RefGene_Group",
            "chr",
            "pos",
            "strand"
        )
        print("int")
        print(int_cols)
        #if (input$select_export_genometype == "hg38" &
        #    (!requireNamespace("rtracklayer", quietly = TRUE)) |
        #    !requireNamespace("GenomicRanges", quietly = TRUE)) {
        #    showModal(
        #        modalDialog(
        #            title = "rtracklayer::liftOver not available",
        #            "Rtracklayer is not installed and it is needed to liftOver annotation from hg19 to hg38 genome. Please, install the package and restart the R session.",
        #            easyClose = TRUE,
        #            footer = NULL
        #        )
        #    )
        #    updateSelectInput(session,
        #                      "select_export_genometype",
        #                      choices = "hg19",
        #                      selected = "hg19"
        #    )
        #}

        showModal(modalDialog(
            title = NULL, footer = NULL,
            div(
                img(src="https://upload.wikimedia.org/wikipedia/commons/7/7d/Pedro_luis_romani_ruiz.gif"),
                p("Generating annotation..."),
                style = "margin: auto; text-align: center"
            )
        ))
        
        withProgress(
            message = "Generating annotation...",
            max = 3,
            value = 1,
            {
                annotation <- as.data.frame(minfi::getAnnotation(rval_gset()))
                annotation <- annotation[, int_cols]
                annotation$genome <- "hg19"
                print("annotation")
                removeModal()
                print(annotation)
                
                #if (input$select_export_genometype == "hg19") {
                #    annotation
                #}
                #else {
                #    chain <- rtracklayer::import.chain(system.file("hg19ToHg38.over.chain", package = "shinyepico"))
                #    
                #    ann_granges <- data.frame(
                #        chr = annotation$chr,
                #        start = annotation$pos - 1,
                #        end = annotation$pos,
                #        name = row.names(annotation)
                #    )
                #    ann_granges <- GenomicRanges::makeGRangesFromDataFrame(
                #        ann_granges,
                #        starts.in.df.are.0based = TRUE,
                #        keep.extra.columns = TRUE
                #    )
                #    ann_granges <- unlist(rtracklayer::liftOver(ann_granges, chain = chain))
                #    
                #    hg38 <- data.table::data.table(
                #        Name = GenomicRanges::mcols(ann_granges)[[1]],
                #        chr = as.character(GenomicRanges::seqnames(ann_granges)),
                #        pos = GenomicRanges::start(ann_granges),
                #        genome = "hg38"
                #    )
                #    
                #    annotation <- as.data.table(annotation[, !(colnames(annotation) %in% c("chr", "pos", "genome"))])
                #    
                #    hg38 <- as.data.frame(data.table::merge.data.table(
                #        x = annotation,
                #        y = hg38,
                #        by = "Name",
                #        all.x = TRUE
                #    ))
                #    row.names(hg38) <- hg38$Name
                #    
                #    hg38
                #}
            }
        )
       
    })
    
    
    
    
    
    
    ##### FUNCTIONAL ENRICHMENT #####
    
    cluster_profiler <- eventReactive(table_annotation(), {
        testgenes <- unique(as.character(table_annotation()$gene))
        print("testgenes")
        testgenes <- testgenes[!testgenes == ""]
        print(testgenes)
        eg <- clusterProfiler::bitr(testgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
        print("eg")
        print(eg)
        eg
    })
    
    kegg <- eventReactive(cluster_profiler(), {
        kegg <- clusterProfiler::enrichKEGG(gene = cluster_profiler()$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
        kegg
    })
    reactome <- eventReactive(cluster_profiler(), {
        reactome <- ReactomePA::enrichPathway(cluster_profiler()$ENTREZID)
        reactome
    })
    go_mf <- eventReactive(cluster_profiler(), {
        mf <- clusterProfiler::enrichGO(gene = cluster_profiler()$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
        mf
    })
    go_bp <- eventReactive(cluster_profiler(), {
        bp <- clusterProfiler::enrichGO(gene = cluster_profiler()$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
        bp
    })
    go_cc <- eventReactive(cluster_profiler(), {
        cc <- clusterProfiler::enrichGO(gene = cluster_profiler()$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
        cc
    })
    gmt_kegg <- eventReactive(cluster_profiler(), {
        print("GMT")
        gmtfile <- "gsea/c2.cp.kegg.v7.1.entrez.gmt"
        print("file")
        print(gmtfile)
        c5 <- suppressWarnings(qusage::read.gmt(gmtfile))
        print("file")
        print(head(c5))
        egmtkegg <- clusterProfiler::enricher(cluster_profiler()$ENTREZID, TERM2GENE=c5)
        print("file")
        print(egmtkegg)
        #write.table(as.data.frame(egmtkegg),"msigdb_kegg.csv",col.names=T,row.names=F,sep="\t",quote=F)
    })
    gmt_go_mf <- eventReactive(cluster_profiler(), {
        gmtfile <- "gsea/c5.mf.v7.1.entrez.gmt"
        c5 <- qusage::read.gmt(gmtfile)
        egmtmf <- clusterProfiler::enricher(cluster_profiler()$ENTREZID, TERM2GENE=c5)
        egmtmf
    })
    gmt_go_bp <- eventReactive(cluster_profiler(), {
        gmtfile <- "gsea/c5.bp.v7.1.entrez.gmt"
        c5 <- qusage::read.gmt(gmtfile)
        egmtbp <- clusterProfiler::enricher(cluster_profiler()$ENTREZID, TERM2GENE=c5)
        egmtbp
    })
    gmt_go_cc <- eventReactive(cluster_profiler(), {
        gmtfile <- "gsea/c5.cc.v7.1.entrez.gmt"
        c5 <- qusage::read.gmt(gmtfile)
        egmtcc <- clusterProfiler::enricher(cluster_profiler()$ENTREZID, TERM2GENE=c5)
        egmtcc
    })
    
    dotplot_kegg <- reactive(enrichplot::dotplot(kegg(), font.size = 12, title = "KEGG"))
    dotplot_go_mf <- reactive(enrichplot::dotplot(go_mf(), font.size = 12, title = "GO - Molecular Function"))
    dotplot_go_bp <- reactive(enrichplot::dotplot(go_bp(), font.size = 12, title = "GO - Biological Process"))
    dotplot_go_cc <- reactive(enrichplot::dotplot(go_cc(), font.size = 12, title = "GO - Cellular Component"))
    dotplot_reactome <- reactive(enrichplot::dotplot(reactome(), font.size = 12, title = "Reactome"))
    dotplot_gmt_kegg <- reactive(enrichplot::dotplot(gmt_kegg(), font.size = 12, title = "MSigDB KEGG"))
    dotplot_gmt_go_mf <- reactive(enrichplot::dotplot(gmt_go_mf(), font.size = 12, title = "MSigDB GO - Molecular Function"))
    dotplot_gmt_go_bp <- reactive(enrichplot::dotplot(gmt_go_bp(), font.size = 12, title = "MSigDB GO - Biological Process"))
    dotplot_gmt_go_cc <- reactive(enrichplot::dotplot(gmt_go_cc(), font.size = 12, title = "MSigDB GO - Cellular Component"))
    
    output$plot_kegg <- renderPlot(dotplot_kegg())
    output$plot_go_mf <- renderPlot(dotplot_go_mf())
    output$plot_go_bp <- renderPlot(dotplot_go_bp())
    output$plot_go_cc <- renderPlot(dotplot_go_cc())
    output$plot_reactome <- renderPlot(dotplot_reactome())
    output$plot_gmt_kegg <- renderPlot(dotplot_gmt_kegg())
    output$plot_gmt_go_mf <- renderPlot(dotplot_gmt_go_mf())
    output$plot_gmt_go_bp <- renderPlot(dotplot_gmt_go_bp())
    output$plot_gmt_go_cc <- renderPlot(dotplot_gmt_go_cc())
    
    
    
    
    
    
    
    
    
    
    
    ###### SURVIVAL #####
    output$clinical_template <- downloadHandler(
        filename = paste("clinical_template", ".csv", sep = ""),
        content = function(file) {write.csv(data.frame(Sample = rval_sheet_target()[[input$select_input_samplenamevar]], Time = NA, Status = NA), file, row.names = FALSE)}
    )
    
    # Input button
    output$ui_clinical_data <- renderUI({
        if (!is.null(input$input_clinical$datapath)) {
            print(input$input_clinical$datapath)
            print(read.csv(input$input_clinical$datapath))
            return(actionButton("b_clinical_data", "Load Clinical Data"))
        } else {
            return()
        }
    })
    
    # input_clinical: enable b_clinical_data
    observeEvent(input$input_clinical, shinyjs::enable("b_clinical_data"))
    
    clinical_sheet <- eventReactive(input$b_clinical_data, {
        print("hello")
        file <- read.csv(input$input_clinical$datapath)
        #validate(need(tools::file_ext(input$input_data$datapath) == "csv", "File extension should be .csv"))
        validate(
            need(
                is.data.frame(file),
                "Clinical Data is not correct. Please, check it."
            )
        )
        #validate(
        #    need(
        #        colnames(file) == rval_sheet_target()[[input$select_input_samplenamevar]],
        #        "Sample names are not correct."
        #    )
        #)
        print(file)
        print(colnames(file))
        file
    })
   
    # When you press b_clinical_data, the form options are updated
    observeEvent(input$b_clinical_data, {
        updateSelectInput(
            session,
            "select_clinical_samplenamevar",
            label = "Select Sample Names Column:",
            choices = colnames(clinical_sheet())
        )
        
        updateSelectInput(
            session,
            "select_clinical_timevar",
            label = "Select Time Column:",
            choices = colnames(clinical_sheet())
        )
        updateSelectInput(
            session,
            "select_clinical_statusvar",
            label = "Select Status Column:",
            choices = colnames(clinical_sheet())
        )
        shinyjs::enable("b_clinical_next") # Enable button continue
    })
    
    output$ui_clinical_different <- renderUI({
        if(input$select_clinical_samplenamevar==input$select_clinical_timevar | input$select_clinical_samplenamevar==input$select_clinical_statusvar | input$select_clinical_timevar==input$select_clinical_statusvar){
            shinyjs::disable("b_clinical_next")
            return(helpText("Columns of Sample, Time and Status must be different"))
        }
        else{
            shinyjs::enable("b_clinical_next")
            return()
        }
    })
    
    # When you press b_clinical_next, the options are updated
    observeEvent(input$b_clinical_next, {
        
        showModal(modalDialog(
            title = NULL, footer = NULL,
            div(
                img(src="https://upload.wikimedia.org/wikipedia/commons/7/7d/Pedro_luis_romani_ruiz.gif"),
                p("Generating survival options..."),
                style = "margin: auto; text-align: center"
            )
        ))
        
        withProgress(
            message = "Generating survival options...",
            max = 3,
            value = 1,
            {
                updateSelectInput(
                    session,
                    "select_clinical_infovar",
                    label = "Select Clinical Information Column:",
                    choices = c("None", colnames(clinical_sheet())[!(colnames(clinical_sheet()) %in% c(input$select_clinical_samplenamevar, input$select_clinical_timevar, input$select_clinical_statusvar))])
                )

                removeModal()
                shinyjs::enable("b_run_survival") # Enable button continue
         }
        )
    })
    
    observeEvent(input$select_clinical_infovar, ignoreNULL = TRUE, {
        if(input$select_clinical_infovar != "None"){
            output$ui_select_clinical_infovar2 <- renderUI(
                selectInput("select_clinical_infovar2", label = "Select Second Clinical Information Column:",
                            choices = c("None", colnames(clinical_sheet())[!(colnames(clinical_sheet()) %in% c(input$select_clinical_samplenamevar, input$select_clinical_timevar, input$select_clinical_statusvar, input$select_clinical_infovar))]))
            )
            if(is.numeric(clinical_sheet()[,input$select_clinical_infovar])){
                print("CLIN1 NUMERIC")
                print(clinical_sheet()[,input$select_clinical_infovar])
                if(length(unique(clinical_sheet()[,input$select_clinical_infovar])) > 2){
                    print("CLIN1 NUMERIC > 2")
                    print(clinical_sheet()[,input$select_clinical_infovar])
                    output$ui_slider_clin1 <- renderUI(sliderInput("slider_clin1",
                                                                   label = paste(input$select_clinical_infovar, "threshold"),
                                                                   min = min(clinical_sheet()[,input$select_clinical_infovar]),
                                                                   max = max(clinical_sheet()[,input$select_clinical_infovar]),
                                                                   value = (min(clinical_sheet()[,input$select_clinical_infovar]) + max(clinical_sheet()[,input$select_clinical_infovar]))/2
                    )
                    )
                }
            }
            else{
                print("CLIN1 NO NUMERIC")
                print(clinical_sheet()[,input$select_clinical_infovar])
                output$ui_slider_clin1 <- renderUI(NULL)
            }
        }
        else{
            output$ui_select_clinical_infovar2 <- renderUI(NULL)
        }
    })
    
    
    
    output$ui_meth_data_disabled <- renderUI({
        if(rval_gset_done() == TRUE){
            updateSwitchInput(session, "select_meth_data", value = TRUE, disabled = FALSE)
            updateSelectizeInput(
                session,
                "select_gene",
                label = "Select Gene:",
                choices = sort(unique(unlist(strsplit(rval_annotation()$UCSC_RefGene_Name, ";"))))
            )
            return()
        }
        else{
            updateSwitchInput(session, "select_meth_data", value = FALSE, disabled = TRUE)
            return(helpText("Methylation options are disabled because normalization is not done.", style = "font-size:11px"))
        }
    })
    
    output$ui_group_island <- renderUI({
        if (input$select_group_island == "Genomic Region") {
            return(selectInput("select_region", label = "Select Genomic Region:", 
                               choices = c("None", sort(unique(unlist(strsplit(rval_annotation()$UCSC_RefGene_Group, ";"))))), 
                   selected = "None"))
        } 
        else if (input$select_group_island == "Relation to Island") {
            return(selectInput("select_island", label = "Select Relation to Island:",
                               choices = c("None", sort(unique(unlist(strsplit(rval_annotation()$Relation_to_Island, ";"))))),
                               selected = "None"))
        }
        else {
            return()
        }
    })
    
   observeEvent(surv_ann_gene(), {
       if(input$select_group_island == "Genomic Region"){
       updateSelectInput(
           session,
           "select_region",
           choices = c("None", sort(unique(unlist(strsplit(surv_annotation()$UCSC_RefGene_Group, ";")))))
       )
       }
       if(input$select_group_island == "Relation to Island"){
       updateSelectInput(
           session,
           "select_island",
           choices = c("None", sort(unique(unlist(strsplit(surv_annotation()$Relation_to_Island, ";")))))
       )
       }
   })
    
   surv_ann_gene <- reactive({
       a_gene <- rval_annotation()[grep(input$select_gene, rval_annotation()$UCSC_RefGene_Name), ]
       print("a_gene")
       print(head(a_gene))
       a_gene
   })
   
   
    surv_annotation <- reactive({
        a1 <- rval_annotation()[grep(input$select_gene, rval_annotation()$UCSC_RefGene_Name), ]
        
        if(input$select_group_island == "Genomic Region"){
        if(input$select_region == "None"){
        a2 <- a1
        }
        else{
            a2 <- a1[grep(input$select_region, a1$UCSC_RefGene_Group), ]
        }
        }
        else{
            a2 <- a1
        }
        
        if(input$select_group_island == "Relation to Island"){
        if(input$select_island == "None"){
            a3 <- a2
        }
        else{
            a3 <- a2[grep(input$select_island, a2$Relation_to_Island), ]
        }
        }
        else{
            a3 <- a2
        }
        a3
    })
    
    graph_survival <- eventReactive(input$b_run_survival, {
        
        surv_time <- input$select_clinical_timevar
        surv_status <- input$select_clinical_statusvar
        surv_clin <- input$select_clinical_infovar
        surv_clin2 <- input$select_clinical_infovar2
        
        new_clinical_sheet <- clinical_sheet()
#        if(!is.null(input$slider_clin1)){
#            print("INSIDE")
#            print(input$slider_clin1)
#            
#            nn <- new_clinical_sheet[,input$select_clinical_infovar]
#            print(nn)
#            nn[nn >= input$slider_clin1] <- paste0(">=", input$slider_clin1)
#            nn[nn < input$slider_clin1] <- paste0("<", input$slider_clin1)
#print(nn)
#print("HELLO")
#            new_clinical_sheet <- cbind(new_clinical_sheet, other = nn)
#            surv_clin <- "other"
#            print(new_clinical_sheet$other)
#            print(new_clinical_sheet$other[new_clinical_sheet$other >= input$slider_clin1])
#            new_clinical_sheet$other[new_clinical_sheet$other >= input$slider_clin1] <- paste0(">=", input$slider_clin1)
#            new_clinical_sheet$other[new_clinical_sheet$other < input$slider_clin1] <- paste0("<", input$slider_clin1)
#            print(new_clinical_sheet)
#            
#         }
#        print(new_clinical_sheet)
        
        
        if(input$select_meth_data == FALSE){
            
            if(surv_clin == "None"){
                showModal(
                    modalDialog(
                        title = "Variable error",
                        "Some variable is needed. Please, select Methylation Data or some Clinical Information.",
                        easyClose = TRUE,
                        footer = NULL
                    )
                )
                return()
            }
            
            else{
                
                surv <- new_clinical_sheet
                print("surv")
                print(surv)
                
                if(surv_clin2 == "None"){
                    log.rank.test <- survival::survdiff(survival::Surv(eval(parse(text = surv_time)), eval(parse(text = surv_status))) ~ eval(parse(text = surv_clin)), data = surv)
                    cox.fit <- summary(survival::coxph(survival::Surv(eval(parse(text = surv_time)), eval(parse(text = surv_status))) ~ eval(parse(text = surv_clin)), data = surv))
                    s_fit <- survminer::surv_fit(survival::Surv(eval(parse(text = surv_time)), eval(parse(text = surv_status)))  ~ eval(parse(text = surv_clin)), data = surv)
                }
                
                else{
                    log.rank.test <- survival::survdiff(survival::Surv(eval(parse(text = surv_time)), eval(parse(text = surv_status))) ~ eval(parse(text = surv_clin)) + eval(parse(text = surv_clin2)), data = surv)
                    cox.fit <- summary(survival::coxph(survival::Surv(eval(parse(text = surv_time)), eval(parse(text = surv_status))) ~ eval(parse(text = surv_clin)) + eval(parse(text = surv_clin2)), data = surv))
                    s_fit <- survminer::surv_fit(survival::Surv(eval(parse(text = surv_time)), eval(parse(text = surv_status)))  ~ eval(parse(text = surv_clin)) + eval(parse(text = surv_clin2)), data = surv)
                }

            }
        }
        
        
        else{
            
            if(input$cpg_input == ""){
                surv_data <- surv_annotation()
            }
            else if(input$cpg_input %in% surv_annotation()$Name){
                surv_data <- surv_annotation()[surv_annotation()$Name == input$cpg_input, ]
            }
            else{
                showModal(
                    modalDialog(
                        title = "CpG site error",
                        "There is not CpG site with that name. Check that the CpG is correctly written.",
                        easyClose = TRUE,
                        footer = NULL
                    )
                )
                return()
            }
            
            ann <- rbind(rval_gset_getBeta()[rownames(surv_data),], beta_median = apply(rval_gset_getBeta()[rownames(surv_data),], 2, median, na.rm = TRUE))
            print("ann")
            print(ann)
            surv <- cbind(new_clinical_sheet, beta_median = t(ann["beta_median",]))
            
            print(surv$beta_median[surv$beta_median >= input$slider_meth_data_val])
            print(surv$beta_median)
            print(surv[beta_median])
            surv$beta_median[surv$beta_median >= input$slider_meth_data_val] <- "Hypermethylated"
            surv$beta_median[surv$beta_median < input$slider_meth_data_val] <- "Hypomethylated"
            #surv$beta_median[surv$beta_median >= input$slider_meth_data_val] <- 1
            #surv$beta_median[surv$beta_median < input$slider_meth_data_val] <- 0
            surv$beta_median <- as.factor(surv$beta_median)
            print("surv")
            print(surv)
            
            
            if(surv_clin == "None"){
                log.rank.test <- survival::survdiff(survival::Surv(eval(parse(text = surv_time)), eval(parse(text = surv_status))) ~ beta_median, data = surv)
                cox.fit <- summary(survival::coxph(survival::Surv(eval(parse(text = surv_time)), eval(parse(text = surv_status))) ~ beta_median, data = surv))
                s_fit <- survminer::surv_fit(survival::Surv(eval(parse(text = surv_time)), eval(parse(text = surv_status)))  ~ beta_median, data = surv)
            }
            
            else{
                
                if(surv_clin2 == "None"){
                    log.rank.test <- survival::survdiff(survival::Surv(eval(parse(text = surv_time)), eval(parse(text = surv_status))) ~ beta_median + eval(parse(text = surv_clin)), data = surv)
                    cox.fit <- summary(survival::coxph(survival::Surv(eval(parse(text = surv_time)), eval(parse(text = surv_status))) ~ beta_median + eval(parse(text = surv_clin)), data = surv))
                    s_fit <- survminer::surv_fit(survival::Surv(eval(parse(text = surv_time)), eval(parse(text = surv_status)))  ~ beta_median + eval(parse(text = surv_clin)), data = surv)
                }
                
                else{
                    log.rank.test <- survival::survdiff(survival::Surv(eval(parse(text = surv_time)), eval(parse(text = surv_status))) ~ beta_median + eval(parse(text = surv_clin)) + eval(parse(text = surv_clin2)), data = surv)
                    cox.fit <- summary(survival::coxph(survival::Surv(eval(parse(text = surv_time)), eval(parse(text = surv_status))) ~ beta_median + eval(parse(text = surv_clin)) + eval(parse(text = surv_clin2)), data = surv))
                    s_fit <- survminer::surv_fit(survival::Surv(eval(parse(text = surv_time)), eval(parse(text = surv_status)))  ~ beta_median + eval(parse(text = surv_clin)) + eval(parse(text = surv_clin2)), data = surv)
                }
            }
        }
        
        #log.rank.test
        print("log.rank.test")
        print(log.rank.test)
        
        #cox.fit
        print("cox.fit")
        print(cox.fit)
        hz <- round(cox.fit$conf.int[1], 3)
        ci.95.low <- round(cox.fit$conf.int[3], 3)
        ci.95.up <- round(cox.fit$conf.int[4], 3)
        hz.pval <- round(cox.fit$logtest[3], 3)
        my_text <- paste0("HR=", hz, "(95% CI ", ci.95.low, " - ", ci.95.up, "); p=", hz.pval)
        print("my_text")
        print(my_text)
        
        ggsurv <- survminer::ggsurvplot(
            fit = s_fit, 
            xlab = input$select_time_unit,
            ylab = "Overall survival probability",
            legend.title = "",
            #legend.labs = c("Hypermethylated", "Hypomethylated"),
            #break.x.by = 10, 
            #palette = ezfun::msk_palette("contrast"), 
            #censor = FALSE,
            risk.table = TRUE,
            risk.table.y.text = TRUE,
            pval = TRUE, 
            pval.method = TRUE)
        ggsurv$plot <- ggsurv$plot + ggplot2::annotate("text", x = 30, y = 1, label = my_text)
        ggsurv
        
    })
    
    output$plot_survival <- renderPlot(graph_survival())
    
    
    
    ##### DOWNLOADS #####
    
    # Markdown Report
    output$download_html <- downloadHandler(
        filename = "Report.html",
        content = function(file) {
            #tempReport <- file.path(tempdir(), "Report.Rmd")
            #print(tempReport)
            #file.copy("report.Rmd", getwd(), overwrite = TRUE)
            src <- normalizePath('report_html.Rmd')
            print(src)
            owd <- setwd(tempdir())
            print(owd)
            on.exit(setwd(owd))
            print("on.exit")
            file.copy(src, 'report_html.Rmd', overwrite = TRUE)
            print("file.copy")
            shinyjs::disable("download_html")
            print("start")
            
            showModal(modalDialog(
                title = NULL, footer = NULL,
                div(
                    img(src="https://upload.wikimedia.org/wikipedia/commons/7/7d/Pedro_luis_romani_ruiz.gif"),
                    p("Generating Report..."),
                    style = "margin: auto; text-align: center"
                )
            ))
            
            withProgress(
                message = "Generating Report...",
                value = 1,
                max = 3,
                {
                    params <- list(
                        rval_sheet = rval_sheet(),
                        rval_sheet_target = rval_sheet_target(),
                        name_var = input$select_input_samplenamevar,
                        grouping_var = input$select_input_groupingvar,
                        donor_var = input$select_input_donorvar,
                        normalization_mode = input$select_minfi_norm,
                        dropsnps = input$select_minfi_dropsnps,
                        dropcphs = input$select_minfi_dropcphs,
                        dropsex = input$select_minfi_chromosomes,
                        maf = input$slider_minfi_maf,
                        probes = rval_gsetprobes(),
                        
                        
                        
                        
                        #limma_voi = input$select_limma_voi,
                        #limma_covar = input$checkbox_limma_covariables,
                        #limma_inter = input$checkbox_limma_interactions,
                        #limma_arrayweights = input$select_limma_weights,
                        #limma_ebayes_trend = input$select_limma_trend,
                        #limma_ebayes_robust = input$select_limma_robust,
                        #rval_design = rval_fit()$design,
                        #rval_contrasts = rval_contrasts(),
                        #rval_voi = rval_voi(),
                        #rval_dendrogram = rval_dendrogram(),
                        #min_deltabeta = input$slider_limma_deltab,
                        #max_fdr = input$slider_limma_adjpvalue,
                        #max_pvalue = input$slider_limma_pvalue,
                        #clusteralg = input$select_limma_clusteralg,
                        #grups2plot = input$select_limma_groups2plot,
                        #contrasts2plot = input$select_limma_contrasts2plot,
                        #Colv = input$select_limma_colv,
                        #distance = input$select_limma_clusterdist,
                        #scale = input$select_limma_scale,
                        #removebatch = input$select_limma_removebatch,
                        
                        
                        
                        
                        plot_green_intensities = boxplot_intensities_green(),
                        plot_red_intensities = boxplot_intensities_red(),
                        plot_failed_probes = failure_plot()[["graph"]],
                        plot_densityplotraw = rval_plot_densityplotraw(),
                        plot_densityplotraw_green = rval_plot_densityplotraw_green(),
                        plot_densityplotraw_red = rval_plot_densityplotraw_red(),
                        plot_densityplotraw_II = rval_plot_densityplotraw_II(),
                        plot_densityplotraw_all = rval_plot_densityplotraw_all(),
                        plot_densityplot = rval_plot_densityplot(),
                        plot_densityplot_green = rval_plot_densityplot_green(),
                        plot_densityplot_red = rval_plot_densityplot_red(),
                        plot_densityplot_II = rval_plot_densityplot_II(),
                        plot_densityplot_all = rval_plot_densityplot_all(),
                        
                        
                        
                        
                        #plot_pcaplot = rval_plot_pca()[["graph"]],
                        
                        
                        
                        plot_corrplot = rval_plot_corrplot()[["graph"]],
                        
                        
                        
                        #plot_boxplotraw = rval_plot_boxplotraw(),
                        #plot_boxplot = rval_plot_boxplot(),
                        #plot_qcraw = rval_plot_qcraw(),
                        #plot_bisulfiterawII = rval_plot_bisulfiterawII(),
                        
                        
                        
                        plot_sexprediction = rval_plot_sexprediction(),
                        plot_snpheatmap = rval_plot_snpheatmap(),
                        
                        
                        
                        #plot_plotSA = rval_plot_plotSA(),
                        #table_pcaplot = rval_plot_pca()[["info"]],
                        
                        
                        
                        table_corrplot = rval_plot_corrplot()[["info"]],
                        data_sexprediction = as.data.frame(minfi::pData(rval_gset()))[["predictedSex"]]
                        
                        
                        
                        #table_dmps = make_table(),
                        #filteredlist2heatmap = rval_filteredlist2heatmap()
                    )
                    
                    print(params)
                    newenv <- new.env(parent = globalenv())
                    #newenv$create_heatmap <- create_heatmap
                    print(newenv)

                    render_file <- rmarkdown::render(
                        #tempReport,
                        input = "report_html.Rmd",
                        output_file = file,
                        run_pandoc = TRUE,
                        params = params,
                        envir = newenv
                    )
                    
                    shinyjs::enable("download_html")
                }
            )
            shinyjs::enable("download_html")
            removeModal()
        }
    )
   
    
    
    output$download_pdf <- downloadHandler(
        filename = "Report.pdf",
        content = function(file) {
            #tempReport <- file.path(tempdir(), "Report.Rmd")
            #print(tempReport)
            #file.copy("report.Rmd", getwd(), overwrite = TRUE)
            src <- normalizePath('report_pdf.Rmd')
            owd <- setwd(tempdir())
            on.exit(setwd(owd))
            file.copy(src, 'report_pdf.Rmd', overwrite = TRUE)
            shinyjs::disable("download_pdf")
            
            showModal(modalDialog(
                title = NULL, footer = NULL,
                div(
                    img(src="https://upload.wikimedia.org/wikipedia/commons/7/7d/Pedro_luis_romani_ruiz.gif"),
                    p("Generating Report..."),
                    style = "margin: auto; text-align: center"
                )
            ))
            
            withProgress(
                message = "Generating Report...",
                value = 1,
                max = 3,
                {
                    params <- list(
                        rval_sheet = rval_sheet(),
                        rval_sheet_target = rval_sheet_target(),
                        name_var = input$select_input_samplenamevar,
                        grouping_var = input$select_input_groupingvar,
                        donor_var = input$select_input_donorvar,
                        normalization_mode = input$select_minfi_norm,
                        dropsnps = input$select_minfi_dropsnps,
                        dropcphs = input$select_minfi_dropcphs,
                        dropsex = input$select_minfi_chromosomes,
                        maf = input$slider_minfi_maf,
                        probes = rval_gsetprobes(),
                        #limma_voi = input$select_limma_voi,
                        #limma_covar = input$checkbox_limma_covariables,
                        #limma_inter = input$checkbox_limma_interactions,
                        #limma_arrayweights = input$select_limma_weights,
                        #limma_ebayes_trend = input$select_limma_trend,
                        #limma_ebayes_robust = input$select_limma_robust,
                        #rval_design = rval_fit()$design,
                        #rval_contrasts = rval_contrasts(),
                        #rval_voi = rval_voi(),
                        #rval_dendrogram = rval_dendrogram(),
                        #min_deltabeta = input$slider_limma_deltab,
                        #max_fdr = input$slider_limma_adjpvalue,
                        #max_pvalue = input$slider_limma_pvalue,
                        #clusteralg = input$select_limma_clusteralg,
                        #grups2plot = input$select_limma_groups2plot,
                        #contrasts2plot = input$select_limma_contrasts2plot,
                        #Colv = input$select_limma_colv,
                        #distance = input$select_limma_clusterdist,
                        #scale = input$select_limma_scale,
                        #removebatch = input$select_limma_removebatch,
                        plot_green_intensities = green_intensities_graph(),
                        plot_red_intensities = red_intensities_graph(),
                        plot_failed_probes = failure_plot(),
                        plot_densityplotraw = rval_plot_densityplotraw(),
                        plot_densityplotraw_green = rval_plot_densityplotraw_green(),
                        plot_densityplotraw_red = rval_plot_densityplotraw_red(),
                        plot_densityplotraw_II = rval_plot_densityplotraw_II(),
                        plot_densityplotraw_all = rval_plot_densityplotraw_all(),
                        plot_densityplot = rval_plot_densityplot(),
                        plot_densityplot_green = rval_plot_densityplot_green(),
                        plot_densityplot_red = rval_plot_densityplot_red(),
                        plot_densityplot_II = rval_plot_densityplot_II(),
                        plot_densityplot_all = rval_plot_densityplot_all(),
                        #plot_pcaplot = rval_plot_pca()[["graph"]],
                        plot_corrplot = rval_plot_corrplot()[["graph"]],
                        #plot_boxplotraw = rval_plot_boxplotraw(),
                        #plot_boxplot = rval_plot_boxplot(),
                        #plot_qcraw = rval_plot_qcraw(),
                        #plot_bisulfiterawII = rval_plot_bisulfiterawII(),
                        plot_sexprediction = rval_plot_sexprediction(),
                        plot_snpheatmap = rval_plot_snpheatmap(),
                        #plot_plotSA = rval_plot_plotSA(),
                        #table_pcaplot = rval_plot_pca()[["info"]],
                        table_corrplot = rval_plot_corrplot()[["info"]],
                        data_sexprediction = as.data.frame(minfi::pData(rval_gset()))[["predictedSex"]]
                        #table_dmps = make_table(),
                        #filteredlist2heatmap = rval_filteredlist2heatmap()
                    )
                    
                    
                    newenv <- new.env(parent = globalenv())
                    #newenv$create_heatmap <- create_heatmap
                    print(newenv)
                    print(params)
                    
                    render_file <- rmarkdown::render(
                        #tempReport,
                        input = "report_pdf.Rmd",
                        output_file = file,
                        run_pandoc = TRUE,
                        params = params,
                        envir = newenv
                    )
                    
                    shinyjs::enable("download_pdf")
                }
            )
            removeModal()
        }
    )
    
    
    
    #Complete report
    output$complete_report_html <- downloadHandler(
        filename = "Report_complete.html",
        content = function(file) {
            #tempReport <- file.path(tempdir(), "Report.Rmd")
            #print(tempReport)
            #file.copy("report.Rmd", getwd(), overwrite = TRUE)
            src <- normalizePath('report_complete_html.Rmd')
            print(src)
            owd <- setwd(tempdir())
            print(owd)
            on.exit(setwd(owd))
            print("on.exit")
            file.copy(src, 'report_complete_html.Rmd', overwrite = TRUE)
            print("file.copy")
            shinyjs::disable("complete_report_html")
            print("start")

            
            showModal(modalDialog(
                title = NULL, footer = NULL,
                div(
                    img(src="https://upload.wikimedia.org/wikipedia/commons/7/7d/Pedro_luis_romani_ruiz.gif"),
                    p("Generating Report..."),
                    style = "margin: auto; text-align: center"
                )
            ))
            
            withProgress(
                message = "Generating Report...",
                value = 1,
                max = 3,
                {
                    params <- list(
                        rval_sheet = rval_sheet(),
                        rval_sheet_target = rval_sheet_target(),
                        name_var = input$select_input_samplenamevar,
                        grouping_var = input$select_input_groupingvar,
                        donor_var = input$select_input_donorvar,
                        normalization_mode = input$select_minfi_norm,
                        dropsnps = input$select_minfi_dropsnps,
                        dropcphs = input$select_minfi_dropcphs,
                        dropsex = input$select_minfi_chromosomes,
                        maf = input$slider_minfi_maf,
                        probes = rval_gsetprobes(),
                        
                        
                        
                        plot_densityplotraw = rval_plot_densityplotraw(),
                        plot_densityplotraw_green = rval_plot_densityplotraw_green(),
                        plot_densityplotraw_red = rval_plot_densityplotraw_red(),
                        plot_densityplotraw_II = rval_plot_densityplotraw_II(),
                        
                        plot_densityplot = rval_plot_densityplot(),
                        plot_densityplot_green = rval_plot_densityplot_green(),
                        plot_densityplot_red = rval_plot_densityplot_red(),
                        plot_densityplot_II = rval_plot_densityplot_II(),
                        
                      
                        #plot_pcaplot = rval_plot_pca()[["graph"]],
                        
                        
                        
                        #plot_boxplotraw = rval_plot_boxplotraw(),
                        #plot_boxplot = rval_plot_boxplot(),
                        #plot_qcraw = rval_plot_qcraw(),
                        #plot_bisulfiterawII = rval_plot_bisulfiterawII(),
                        
                        plot_sexprediction = rval_plot_sexprediction()
                        
                        
                        #plot_plotSA = rval_plot_plotSA(),
                        #table_pcaplot = rval_plot_pca()[["info"]],
                        
                        
                        
                        #exploratory analysis
                        

                    )
                    
                    if (input$check_qc){
                        params[["check_qc"]] <- TRUE
                        if("1" %in% input$check_group_qc){
                            params[["plot_green_intensities"]] <- boxplot_intensities_green()
                            params[["plot_red_intensities"]] <- boxplot_intensities_red()
                        }
                        if("2" %in% input$check_group_qc){
                            params[["plot_failed_probes"]] <- failure_plot()[["graph"]]
                        }
                        if("3" %in% input$check_group_qc){
                            params[["plot_densityplotraw_all"]] <- rval_plot_densityplotraw_all()
                            params[["plot_densityplot_all"]] <- rval_plot_densityplot_all()
                        }
                        if("4" %in% input$check_group_qc){
                            params[["plot_snpheatmap"]] <- rval_plot_snpheatmap()
                        }
                        if("5" %in% input$check_group_qc){
                            params[["data_sexprediction"]] <- as.data.frame(minfi::pData(rval_gset()))[["predictedSex"]]
                        }
                        if("6" %in% input$check_group_qc){
                            params[["plot_corrplot"]] <- rval_plot_corrplot()[["graph"]]
                            params[["table_corrplot"]] <- rval_plot_corrplot()[["info"]]
                        }
                    }
                    
                    if(input$check_exploratory_analysis){
                        params[["check_exploratory_analysis"]] <- TRUE
                        if("1" %in% input$check_group_exploratory_analysis){
                            params[["plot_violin_raw"]] <- rval_plot_violin_raw()
                            params[["plot_violin_normalized"]] <- rval_plot_violin_normalized()
                        }
                        if("2" %in% input$check_group_exploratory_analysis){
                            params[["plot_pca"]] <- rval_plot_pca()[["graph"]]
                        }
                        if("3" %in% input$check_group_exploratory_analysis){
                            params[["plot_random_heatmap"]] <- plot_random_heatmap()
                            params[["plot_top_heatmap"]] <- plot_top_heatmap()
                        }
                        if("4" %in% input$check_group_exploratory_analysis){
                            params[["plot_deconvolution"]] <- graph_deconvolution()
                        }
                        if("5" %in% input$check_group_exploratory_analysis){
                            params[["table_age"]] <- age_data()
                        }
                        if("6" %in% input$check_group_exploratory_analysis){
                            params[["plot_hyper_hypo_chr"]] <- graph_hyper_hypo()[["chr"]]
                            params[["plot_hyper_hypo_relation_to_island"]] <- graph_hyper_hypo()[["relation_to_island"]]
                            params[["plot_hyper_hypo_group"]] <- graph_hyper_hypo()[["group"]]
                        }
                    }

                    # if DMP analysis has been done, we add specific parameters
                    if (rval_analysis_finished()) {
                        params[["limma_voi"]] <- input$select_limma_voi
                        params[["limma_covar"]] <- input$checkbox_limma_covariables
                        params[["limma_inter"]] <- input$checkbox_limma_interactions
                        params[["limma_arrayweights"]] <- input$select_limma_weights
                        params[["limma_ebayes_trend"]] <- input$select_limma_trend
                        params[["limma_ebayes_robust"]] <- input$select_limma_robust
                        params[["rval_design"]] <- rval_fit()$design
                        params[["rval_contrasts"]] <- rval_contrasts()
                        params[["rval_voi"]] <- rval_voi()
                        params[["rval_dendrogram"]] <- rval_dendrogram()
                        params[["min_deltabeta"]] <- input$slider_limma_deltab
                        params[["max_fdr"]] <- input$slider_limma_adjpvalue
                        params[["max_pvalue"]] <- input$slider_limma_pvalue
                        params[["clusteralg"]] <- input$select_limma_clusteralg
                        params[["groups2plot"]] <- input$select_limma_groups2plot
                        params[["contrasts2plot"]] <- input$select_limma_contrasts2plot
                        params[["Colv"]] <- input$select_limma_colv
                        params[["distance"]] <- input$select_limma_clusterdist
                        params[["scale"]] <- input$select_limma_scale
                        params[["removebatch"]] <- input$select_limma_removebatch
                        params[["table_dmps"]] <- make_table()
                        params[["filteredlist2heatmap"]] <- rval_filteredlist2heatmap()
                        params[["table_annotation_manhattan"]] <- table_annotation_manhattan()
                        params[["plot_volcano"]] <- volcano_graph()
                        params[["table_annotation"]] <- table_annotation()
                    }
                    
                    
                    
                    # if DMR analysis has been done, we add specific parameters
                    if (rval_dmrs_finished()) {
                        params[["dmrs_contrasts"]] <- input$select_dmrs_contrasts
                        params[["dmrs_rval_dendrogram"]] <- rval_dendrogram_dmrs()
                        params[["dmrs_min_deltabeta"]] <- input$slider_dmrs_deltab
                        params[["dmrs_max_fdr"]] <- input$slider_dmrs_adjpvalue
                        params[["dmrs_max_pvalue"]] <- input$slider_dmrs_pvalue
                        params[["dmrs_clusteralg"]] <- input$select_dmrs_clusteralg
                        params[["dmrs_groups2plot"]] <- input$select_dmrs_groups2plot
                        params[["dmrs_contrasts2plot"]] <- input$select_dmrs_contrasts2plot
                        params[["dmrs_regions2plot"]] <- input$select_dmrs_regions2plot
                        params[["dmrs_Colv"]] <- input$select_dmrs_colv
                        params[["dmrs_distance"]] <- input$select_dmrs_clusterdist
                        params[["dmrs_scale"]] <- input$select_dmrs_scale
                        params[["dmrs_removebatch"]] <- input$select_dmrs_removebatch
                        params[["table_dmrs"]] <- make_table_dmrscount()
                        params[["filteredmcsea2heatmap"]] <- rval_filteredmcsea2heatmap()
                        params[["table_sigdmrs"]] <- rval_table_sigdmrs()
                    }
                    
                    # functional analysis
                    if(rval_analysis_finished()) {
                        params[["kegg"]] <- dotplot_kegg()
                        params[["go_mf"]] <- dotplot_go_mf()
                        params[["go_bp"]] <- dotplot_go_bp()
                        params[["go_cc"]] <- dotplot_go_cc()
                        params[["reactome"]] <- dotplot_reactome()
                        #params[["gmt_kegg"]] <- dotplot_gmt_kegg()
                        #params[["gmt_go_mf"]] <- dotplot_gmt_go_mf()
                        #params[["gmt_go_bp"]] <- dotplot_gmt_go_bp()
                        #params[["gmt_go_cc"]] <- dotplot_gmt_go_cc()
                    }
                    
                    
                    

                    newenv <- new.env(parent = globalenv())
                    #newenv$create_heatmap <- create_heatmap
                    
                    render_file <- rmarkdown::render(
                        #tempReport,
                        input = "report_complete_html.Rmd",
                        output_file = file,
                        run_pandoc = TRUE,
                        params = params,
                        envir = newenv
                    )
                    
                    shinyjs::enable("complete_report_html")
                }
            )
            removeModal()
        }
    )
    
    
    
})
