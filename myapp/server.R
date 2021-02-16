
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
    green_intensities <- reactive(stack(as.data.frame(getGreen(rval_rgset()))))
    red_intensities <- reactive(stack(as.data.frame(getRed(rval_rgset()))))
    
    green_intensities_graph <- reactive(ggplot2::ggplot(green_intensities(), ggplot2::aes(x = green_intensities()$ind, y = green_intensities()$values)) +
                         ggplot2::geom_boxplot(fill = "darkgreen", outlier.shape = NA, show.legend = FALSE) +
                         ggplot2::ylim(boxplot.stats(green_intensities()$values)$stats[1], boxplot.stats(green_intensities()$values)$stats[5]) +
                         ggplot2::theme_bw() +
                         ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90), 
                                        axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank()))
    
    red_intensities_graph <- reactive(ggplot2::ggplot(red_intensities(), ggplot2::aes(x = red_intensities()$ind, y = red_intensities()$values)) +
        ggplot2::geom_boxplot(fill = "red", outlier.shape = NA, show.legend = FALSE) +
        ggplot2::ylim(boxplot.stats(red_intensities()$values)$stats[1], boxplot.stats(red_intensities()$values)$stats[5]) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90), 
                       axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank()))
    
    output$green_intensities_plot <- renderPlot(green_intensities_graph())
    output$red_intensities_plot <- renderPlot(red_intensities_graph())
    
    ########## FAILURE RATE PLOT ##########
    
    beta_pvalue <- function(rgset, betas){
        pval <- as.matrix(minfi::detectionP(rgset))
        beta_pval <- as.matrix(betas)
        beta_pval[pval>=0.01] <- NA
        fail <- as.data.frame(cbind(sort((apply(beta_pval,2,function(x) sum(is.na(x)))/nrow(beta_pval)*100))))
        colnames(fail)[1] <- "probe_failure_rate"
        
        fail
    }
    
    failed_beta <- reactive(beta_pvalue(rval_rgset(), rval_rgset_getBeta()))
    
    failure_graph <- reactive(plotly::ggplotly(ggplot(failed_beta(), aes(x = rownames(failed_beta()), y = failed_beta()$probe_failure_rate)) +
                         geom_bar(stat = "identity", fill = "steelblue") +
                         scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0, 100, 5)) +  ###max(fail$probe_failure_rate) + 0.5
                         geom_hline(yintercept = 5, linetype = "dashed", color = "red", size = 0.5) + 
                         geom_hline(yintercept = 10, linetype = "dashed", color = "red", size = 0.5) +
                         xlab("Sample Name") + ylab("% Probe Failure Rate") + 
                         coord_flip() +
                         annotation_logticks(base = 2, sides = "bl") +
                         theme_bw()))
    
    failure_plot <- reactive(failure_graph())
    output$failure_rate_plot <- plotly::renderPlotly(failure_plot())
    
    ########## CONTROL TYPE PLOTS ##########
    
    output$controlTypePlotGreen <- renderPlot({
        if (!is.null(input$controlType)){
            
            groupNames <- rval_sheet_target()$Sample_Group
            sampleNames <- rval_sheet_target()$Samples_Name #sampleNames(shinyMethylSet1())
            
            if (input$controlType %in% c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II", "HYBRIDIZATION", "SPECIFICITY I",
                                         "SPECIFICITY II", "TARGET REMOVAL")){
                threshold <- 1
            } else if (input$controlType %in% c("EXTENSION", "STAINING", "NON-POLYMORPHIC")){
                threshold <- 5 # you can increase the threshold
                
            } else {threshold <- 0}
            
            green <- getGreen(rval_rgset())
            ctrlAddress <- getControlAddress(rval_rgset(), controlType = input$controlType)
            
            green_control <- log2(green[ctrlAddress, ,drop = FALSE])
            
            slide <- input$select_slide
            title <- paste("-", slide)
            
            subset <- reshape2::melt(green_control)
            subset2 <- subset[grepl(slide, subset$Var2), ]
            
            ggplot(data=as.data.frame(subset2), aes(x=Var2, y=value)) +
                geom_point(color="darkgreen", size=1.5) + scale_y_continuous(limits = c(-1, 20)) +
                theme(axis.text.x = element_text(hjust = 1, angle=45)) +
                geom_hline(yintercept =threshold, linetype="dashed") + ylab("Log2 Intensity") +
                scale_x_discrete(labels=groupNames) + xlab("Samples") + ggtitle(paste("Green Channel", title))
        }
    })
    
    output$controlTypePlotRed <- renderPlot({
        if (!is.null(input$controlType)){
            groupNames <- rval_sheet_target()$Sample_Group
            sampleNames <- rval_sheet_target()$Samples_Name
            
            if (input$controlType %in% c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II", "HYBRIDIZATION", "SPECIFICITY I",
                                         "SPECIFICITY II", "TARGET REMOVAL")){
                threshold <- 1
            } else if (input$controlType %in% c("EXTENSION", "STAINING", "NON-POLYMORPHIC")){
                threshold <- 5 # you can increase the threshold
                
            } else {threshold <- 0}
            
            red <- getRed(rval_rgset())
            ctrlAddress <- getControlAddress(rval_rgset(), controlType = input$controlType)
            
            red_control <- log2(red[ctrlAddress, ,drop = FALSE])
            
            slide <- input$select_slide
            title <- paste("-", slide)
            
            subset <- reshape2::melt(red_control)
            subset2 <- subset[grepl(slide, subset$Var2), ]
            
            ggplot(data=as.data.frame(subset2), aes(x=Var2, y=value)) +
                geom_point(color="red", size=1.5) + scale_y_continuous(limits = c(-1, 20)) +
                theme(axis.text.x = element_text(hjust = 1, angle=45)) +
                geom_hline(yintercept =threshold, linetype="dashed") + ylab("Log2 Intensity") +
                scale_x_discrete(labels=groupNames) + xlab("Samples") + ggtitle(paste("Red Channel", title))
        } })
    
    
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
    
    rval_plot_violin <- reactive(create_violinplot(rval_rgset_getBeta(), nrow(rval_rgset_getBeta())))
    output$graph_violin <- renderPlot(rval_plot_violin())
    
    
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
