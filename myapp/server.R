
library(shiny)
source("utils_analysis.R")
source("utils_graphs.R")
source("utils_download.R")
library(ggplot2)


shinyServer(function(input, output, session) {
    
    # INITIALIZE REACTIVE VARIABLES
    rval_generated_limma_model <- reactiveVal(value = FALSE)
    rval_analysis_finished <- reactiveVal(value = FALSE)
    rval_filteredlist2heatmap_valid <- reactiveVal(value = FALSE)
    rval_filteredmcsea2heatmap_valid <- reactiveVal(value = FALSE)
    rval_dmrs_finished <- reactiveVal(value = FALSE)
    rval_dmrs_ready2heatmap <- reactiveVal(value = FALSE)
    rval_dmrs_ready2mcsea <- reactiveVal(value = FALSE)
    
    
    # Max size
    options(shiny.maxRequestSize = 8000 * 1024^2) # 5MB getShinyOption("shiny.maxRequestSize") | 30*1024^2 = 30MB
    
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
        
        updatePickerInput(
            session,
            "selected_samples_h",
            label = "Select Samples to Plot:",
            selected = rval_sheet_target()[, input$select_input_samplenamevar],
            choices = rval_sheet_target()[, input$select_input_samplenamevar],
            choicesOpt = list(subtext = paste("Group: ", rval_sheet_target()[, input$select_input_groupingvar]))
        )
        
        shinyjs::enable("button_minfi_select")
    })
    
    
    ########## rval_gset() NORMALIZATION ##########
    
    rval_gset <- eventReactive(input$button_minfi_select, {
        validate(need(
            !is.null(rval_rgset()),
            "Raw data has not been loaded yet."
        ))
        
        shinyjs::disable("button_minfi_select") # disable button to avoid repeat clicking
        
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
                shinyjs::enable("button_minfi_select")
                
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
        bvalues
    })
    
    rval_gset_getM <- reactive({
        mvalues <- minfi::getM(rval_gset())
        colnames(mvalues) <- rval_sheet_target()[[input$select_input_samplenamevar]]
        mvalues
    }) 
    
    
    ##### INTENSITIES BOXPLOTS #####

    output$green_intensities_plot <- renderPlot(create_boxplot_intensities_green(rval_rgset()))
    output$red_intensities_plot <- renderPlot(create_boxplot_intensities_red(rval_rgset()))
    
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
    
    output$table_age <- DT::renderDT(create_age(rval_rgset_getBeta(), rval_sheet_target(), input$select_input_age))
    
    
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
    
    
    
    
    
    
    
    
    # LIMMA
    
    # Variable of interest
    rval_voi <- reactive(factor(make.names(minfi::pData(rval_gset())[, input$select_limma_voi])))
    
    # Design calculation
    rval_design <- eventReactive(input$button_limma_calculatemodel, {
        print("HELLO")
        req(input$select_limma_voi) # a variable of interest is required
        print(input$select_limma_voi)
        print(rval_voi)
        print("LIMMA")
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
        
        design
    })
    
    # Calculation of contrasts
    rval_contrasts <- reactive({
        generate_contrasts(rval_voi())
    })
    
    
    # Calculation of limma model
    rval_fit <- eventReactive(input$button_limma_calculatemodel, ignoreNULL = FALSE, {
        validate(
            need(input$fileinput_input != "", "DMP calculation has not been performed or data has not been uploaded.")
        )
        
        req(rval_design())
        
        shinyjs::disable("button_limma_calculatemodel") # disable button to avoid repeat clicking
        
        withProgress(
            message = "Generating linear model...",
            value = 3,
            max = 6,
            {
                try({
                    fit <- generate_limma_fit(
                        Mvalues = rval_gset_getM(), design = rval_design(),
                        weighting = as.logical(input$select_limma_weights)
                    )
                })
                
                if (!exists("fit", inherits = FALSE)) {
                    rval_generated_limma_model(FALSE) # disable contrast button
                    rval_analysis_finished(FALSE) # disable download buttons
                    
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
        
        validate(need(
            exists("fit", inherits = FALSE),
            "lmFit model has failed. Please, check your options and try again."
        ))
        print(fit)
        
        fit # returning linear model
    })
    
    
    # rval_fit() has NAs, we remove the option to trend or robust in eBayes to prevent failure
    observeEvent(input$button_limma_calculatemodel, {
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
        rval_design(),
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
        
        globaldifs
    })
    
    # Calculation of differences (eBayes)
    rval_finddifcpgs <- eventReactive(input$button_limma_calculatedifs, {
        try({
            dif_cpgs <- find_dif_cpgs(
                design = rval_design(),
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
            "select_anncontrast",
            label = "Contrast",
            choices = rval_contrasts()
        )
        
        # disable button to avoid repeat n_coresicking
        shinyjs::disable("button_limma_calculatedifs")
        
        # force rval_filteredlist
        rval_filteredlist()
        
        # enable or disable removebatch option
        covariables_design <- as.matrix(rval_design()[, -seq_len(length(unique(rval_voi())))])
        
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
    rval_filteredlist <- reactive({
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
            rval_design(),
            rval_voi(),
            rval_gset_getBeta()
        )
        
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
    
    rval_cpgcount_heatmap <- eventReactive(input$button_limma_heatmapcalc, nrow(rval_filteredlist2heatmap()))
    
    rval_dendrogram <- eventReactive(input$button_limma_heatmapcalc, {
        if (input$select_limma_rowsidecolors) {
            
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
    
    plot_heatmap <- eventReactive(input$button_limma_heatmapcalc, {
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
    
    make_table <- eventReactive(input$button_limma_heatmapcalc, {
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
    
    table_annotation <- eventReactive(list(input$button_limma_heatmapcalc, input$select_limma_anncontrast), {
        req(rval_filteredlist())
        
        print("rval_filteredlist")
        print(head(rval_filteredlist()))
        
        dif_target <- paste("dif",
                            limma::strsplit2(input$select_limma_anncontrast, "-")[1],
                            limma::strsplit2(input$select_limma_anncontrast, "-")[2],
                            sep = "_"
        )
        print(dif_target)
        print("annotation")
        print(head(rval_annotation()))
        print("globaldifs")
        print(head(rval_globaldifs()))
        
        
        temp <- rval_annotation()[row.names(rval_annotation()) %in% rval_filteredlist()[[input$select_limma_anncontrast]]$cpg, ]
        print(head(temp))
        print(nrow(rval_annotation()))
        print(nrow(rval_filteredlist()[[input$select_limma_anncontrast]]$cpg))
        print(nrow(temp))
        temp$dif_beta <- rval_globaldifs()[[dif_target]][rval_globaldifs()[["cpg"]] %in% row.names(temp)]
        print(head(temp))
        print(nrow(rval_globaldifs()[[dif_target]]))
        print(nrow(rval_globaldifs()[["cpg"]]))
        print(nrow(temp))
        temp$fdr <- rval_filteredlist()[[input$select_limma_anncontrast]][["adj.P.Val"]][rval_filteredlist()[[input$select_limma_anncontrast]][["cpg"]] %in% row.names(temp)]
        print(head(temp))
        print(nrow(rval_filteredlist()[[input$select_limma_anncontrast]][["adj.P.Val"]]))
        print(nrow(rval_filteredlist()[[input$select_limma_anncontrast]][["cpg"]]))
        print(nrow(temp))
        temp$pvalue <- rval_filteredlist()[[input$select_limma_anncontrast]][["P.Value"]][rval_filteredlist()[[input$select_limma_anncontrast]][["cpg"]] %in% row.names(temp)]
        print(head(temp))
        print(nrow(rval_filteredlist()[[input$select_limma_anncontrast]][["P.Value"]]))
        print(nrow(rval_filteredlist()[[input$select_limma_anncontrast]][["cpg"]]))
        print(nrow(temp))
        temp$chr <- as.numeric(as.character(gsub("chr", "", temp$chr)))
        gene <- vapply(strsplit(temp$UCSC_RefGene_Name,";"), `[`, 1, FUN.VALUE=character(1))
        gene[is.na(gene)]<-""
        print(head(gene))
        temp$gene <- gene
        
        print("temp")
        print(head(temp))
        
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
    
    
    ########## MANHATTAN PLOT ##########
    
    table_annotation2 <- eventReactive(input$select_anncontrast, {
        req(rval_list())
        
        print("rval_list")
        print(head(rval_list()))
        
        dif_target <- paste("dif",
                            limma::strsplit2(input$select_anncontrast, "-")[1],
                            limma::strsplit2(input$select_anncontrast, "-")[2],
                            sep = "_"
        )
        print(dif_target)
        print("annotation")
        print(head(rval_annotation()))
        print("globaldifs")
        print(head(rval_globaldifs()))
        
        
        temp <- rval_annotation()[row.names(rval_annotation()) %in% rval_list()[[input$select_anncontrast]]$cpg, ]
        print(head(temp))
        print(nrow(rval_annotation()))
        print(nrow(rval_list()[[input$select_anncontrast]]$cpg))
        print(nrow(temp))
        print(which(is.na(temp$Name)))
        temp$dif_beta <- rval_globaldifs()[[dif_target]][rval_globaldifs()[["cpg"]] %in% row.names(temp)]
        print(head(temp))
        print(nrow(rval_globaldifs()[[dif_target]]))
        print(nrow(rval_globaldifs()[["cpg"]]))
        print(nrow(temp))
        print(which(is.na(temp$dif_beta)))
        
        temp$fdr <- rval_list()[[input$select_anncontrast]][["adj.P.Val"]][rval_list()[[input$select_anncontrast]][["cpg"]] %in% row.names(temp)]
        print(head(temp))
        print(nrow(rval_list()[[input$select_anncontrast]][["adj.P.Val"]]))
        print(nrow(rval_list()[[input$select_anncontrast]][["cpg"]]))
        print(nrow(temp))
        print(which(is.na(temp$fdr)))
        
        temp$pvalue <- rval_list()[[input$select_anncontrast]][["P.Value"]][rval_list()[[input$select_anncontrast]][["cpg"]] %in% row.names(temp)]
        print(head(temp))
        print(nrow(rval_list()[[input$select_anncontrast]][["P.Value"]]))
        print(nrow(rval_list()[[input$select_anncontrast]][["cpg"]]))
        print(nrow(temp))
        print(which(is.na(temp$pvalue)))
        
        print("CHR")
        
        temp$chr[temp$chr == "chrX"] <- 23
        temp$chr[temp$chr == "chrY"] <- 24
        #temp$chr[temp$chr == "chrM"] <- 25
        
        
        
        print(unique(temp$chr))
        
        
        temp$chr <- as.numeric(as.character(gsub("chr", "", temp$chr)))
        gene <- vapply(strsplit(temp$UCSC_RefGene_Name,";"), `[`, 1, FUN.VALUE=character(1))
        gene[is.na(gene)]<-""
        print(head(gene))
        temp$gene <- gene
        print(which(is.na(temp$gene)))
        print(unique(temp$chr))
        
        print("temp")
        print(head(temp))
        
        temp
    })
    
    volcano_data <- eventReactive(table_annotation2(), {
        pval <- table_annotation2()$pvalue
        fc <- table_annotation2()$dif_beta
        names <- table_annotation2()$gene
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
    
    #volcano_graph <- reactive(plotly::ggplotly(ggplot2::ggplot(volcano_data(), ggplot2::aes_string(x="FC", y="PV", color="clr", fill="clr")) +
    #                           ggplot2::theme_bw() +
    #                          ggplot2::geom_point(alpha=1) #+
    #ggplot2::scale_colour_manual(values=clrvalues) +
    #ggplot2::ylab("expression(-log[10](P-Value))") +
    #ggplot2::theme(legend.position="none") +
    #ggplot2::xlab("expression(log[2](Fold~~Change))") +
    #ggrepel::geom_text_repel(
    #  data = subset(volcano_data(), volcano_data()$PV >= tPV & abs(volcano_data()$FC) >= tFC),
    #  ggplot2::aes_string("FC", "PV", label="names"),
    #  size = 2,
    #  box.padding = ggplot2::unit(0.35, "lines"),
    #  point.padding = ggplot2::unit(0.3, "lines"),
    #  color="black"
    #) +
    #ggplot2::geom_hline(yintercept=tPV,
    #                    linetype="dotdash", color="gray69", size=0.75) +
    #ggplot2::geom_vline(xintercept=-tFC,
    #                    linetype="dotdash", color="gray69", size=0.75) +
    #ggplot2::geom_vline(xintercept=tFC,
    #                    linetype="dotdash", color="gray69", size=0.75)))
    #))
    manhattan_graph <- reactive(qqman::manhattan(table_annotation2(), chr = "chr", bp = "pos", snp = "gene", p = "pvalue",
                                                 annotatePval = 1, suggestiveline = T, genomewideline = T, annotateTop = T))
    volcano_graph <- reactive(MultiDataSet::volcano_plot(pval = table_annotation2()$pvalue, fc = table_annotation2()$dif_beta,
                                                         table_annotation2()$gene, tFC = 0.2, show.labels = T))
    
    
    output$manhattan_plot <- renderPlot(manhattan_graph())
    output$volcano_plot <- renderPlot(volcano_graph())
    #output$volcano_plot1 <- plotly::renderPlotly(volcano_graph1())
    
    
    
    
    
    
    ########## VOLCANO PLOT ##########
    
    
    
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
    
    
    # DMRs
    
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
    
    rval_filteredmcsea <- eventReactive(list(input$button_dmrs_calculate, input$button_dmrs_heatmapcalc), {
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
            design = rval_design(),
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
                input$fileinput_input != "",
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ##### DOWNLOADS #####
    
    # Markdown Report
    output$download_html <- downloadHandler(
        filename = "Report.html",
        content = function(file) {
            #tempReport <- file.path(tempdir(), "Report.Rmd")
            #print(tempReport)
            #file.copy("report.Rmd", getwd(), overwrite = TRUE)
            src <- normalizePath('report_html.Rmd')
            owd <- setwd(tempdir())
            on.exit(setwd(owd))
            file.copy(src, 'report_html.Rmd', overwrite = TRUE)
            shinyjs::disable("download_html")
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
                        #rval_design = rval_design(),
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
                        input = "report_html.Rmd",
                        output_file = file,
                        run_pandoc = TRUE,
                        params = params,
                        envir = newenv
                    )
                    
                    shinyjs::enable("download_html")
                }
            )
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
                        #rval_design = rval_design(),
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
        }
    ) 
    
    
    
})
