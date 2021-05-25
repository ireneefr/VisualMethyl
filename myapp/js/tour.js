
var intro = introJs();

Shiny.addCustomMessageHandler("intro_steps",

  function(message){
    intro.setOptions({
      exitOnOverlayClick: false,
      showBullets: false,
      disableInteraction: true,
      steps: [
      
      {
        title: "<h1 style='text-align:center;'>Welcome to <b>VisualMethyl</b></h1><br><h4 style='text-align:center;'>A shiny web application for interactive DNA methylation analysis and visualization</h4>",
        intro: "This is ................."
      },
      {
        element: document.querySelector("#ui_input_data"),
        title: "Input Data",
        intro: "First, you have to upload the data, which must be in .zip format containing the IDAT files and the sample sheet in .csv.<br>In this case, we are going to load the example data."
      },
      {
        element: document.querySelector("#div_samples_table"),
        title: "Samples Table",
        intro: "In this box, there is the Samples Table where you can inspect the data loaded."
      },
      {
        element: document.querySelector("#div_select_options"),
        title: "Selection Options",
        intro: "Before the analysis, you have to select the correct column of the sample sheet in each option.<br>In this example, the grouping variable is <i>Sample_Group</i> and there is no Donor/Patient so is selected <i>Sample_Name</i> (Sample Names column).",
        position: "right"
      },
      {
        intro: ""
      }
      
      
      /*,
      
      {
        element: document.getElementById('startButton').onclick = 
        function() {
       introJs().onchange(function(targetElement) {
          if (this._currentStep==0) {
             $('a[data-value=\"Second tab\"]').removeClass('active');
             $('a[data-value=\"First tab\"]').addClass('active');
             $('a[data-value=\"First tab\"]').trigger('click');
          }
          if (this._currentStep==1) {
             $('a[data-value=\"First tab\"]').removeClass('active');
             $('a[data-value=\"Second tab\"]').addClass('active');
             $('a[data-value=\"Second tab\"]').trigger('click');
          }
       })
       },
        intro:
      },{
        element: document.querySelector("#button_input_next"),
        intro: "Continue to analysis"
      },
      {
        element: document.querySelector(""),
        intro:
      }
      */
      
    ]});
    
  }
);


Shiny.addCustomMessageHandler("intro_start", 
  function(message){
    
    intro.onchange(function(targetElement) {
/*
          if (intro._currentStep==1) {
            //document.querySelector(".introjs-nextbutton").style.display = "none";
            if ($("#button_input_next").prop("disabled")){
              alert("yes");
            }
            if ($("#button_input_next").prop("enabled")){
              alert("no");
            }
            document
                .querySelector(".introjs-donebutton")
                .classList.remove("introjs-skipbutton");
            
             document.querySelector(".introjs-nextbutton").style.display =
                "none";
               
          } */
          if(this._currentStep==0){
            $('a[data-value=\"data\"]').trigger('click');
            intro.setOptions({'nextLabel': 'Start tour', "disableInteraction": "true"});
          }
          if(this._currentStep==1){
            intro.setOptions({'nextLabel': 'Load Data', "disableInteraction": "true"});
            document.querySelector(".introjs-prevbutton").style.display = "none";
            $("#b_input_data").click();
            
          }
          if(this._currentStep==2){
            
            intro.setOptions({'nextLabel': 'Next', "disableInteraction": "false"});
            document.querySelector(".introjs-prevbutton").style.display = "block";
          }
          if(this._currentStep==3){
            intro.setOptions({'nextLabel': 'Continue to Analysis', "disableInteraction": "true"});
            document.querySelector(".introjs-prevbutton").style.display = "block";
            /*
            if ($("#button_input_next").prop("disabled")){
              document.querySelector(".introjs-nextbutton").style.display = "none";
            }
            if ($("#button_input_next").prop("enabled")){
              document.querySelector(".introjs-nextbutton").style.display = "block";
            }*/
          }

          if (this._currentStep==4) {
             //$('a[data-value=\"data\"]').removeClass('active');
             //$('a[data-value=\"analysis\"]').addClass('active');
             $('#button_input_next').click();
             //$('#shiny_modal').show();
             intro.exit();
             //Shiny.onInputChange("inside_tour", change);
             //$('a[data-value=\"analysis\"]').trigger('click');
          }
       }).start();
  }
);









Shiny.addCustomMessageHandler("intro_steps_continue1",

  function(message){
    intro.setOptions({
      exitOnOverlayClick: false,
      showBullets: false,
      disableInteraction: true,
      steps: [
      
      {
        element: document.querySelector("#b_qc"),
        intro: "Next Step"
      },
      {
        element: document.querySelector("#div_norm"),
        title: "Normalization Options",
        intro: "To do the normalization, you have to select the type of normalization that you want to apply and determine other parameters. <br> For the current data is selected the Illumina normalization and the other parameters with the default values. <hr> <small style='color:steelblue;'> For more information you can consult <b>Help</b> section </small>",
        position: "right"
      },
      {
        intro: ""
      }
      
    ]});
    
  }
);




Shiny.addCustomMessageHandler("intro_continue1", 
  function(message){
    intro.onchange(function(targetElement) {
      if (intro._currentStep==0){
        $('a[data-value=\"analysis\"]').trigger('click');
        intro.setOptions({'nextLabel': 'Go to Quality Control', "disableInteraction": "true"});
      }
          if(intro._currentStep==1){
            $('a[data-value=\"qc\"]').trigger('click');
            intro.setOptions({'nextLabel': 'Run Normalzation'});
            //document.querySelector(".introjs-prevbutton").style.display = "none";
          }
          if (intro._currentStep==2){
            $("#button_minfi_select").click();
            intro.exit();
          }
       }).start();
  }
);



Shiny.addCustomMessageHandler("intro_steps_continue2",

  function(message){
    intro.setOptions({
      exitOnOverlayClick: false,
      showBullets: false,
      disableInteraction: false,
      steps: [
      
      {
        element: document.querySelector("#vertical"),
        title: "Quality Control Plots",
        intro: "There are several graphics to inspect the data and compare raw data with normalized data. <hr> <b style='color:steelblue;'>Click + to expand the boxes and see the content</b>",
        position: "left"
      },
      {
        element: document.querySelector("#b_exploratory_analysis"),
        intro: "Next step"
      },
      {
        element: document.querySelector("#exploratory_plots"),
        title: "Exploratory Analysis Plots",
        intro: "Other plots to examine the data in more detail. <hr> <b style='color:steelblue;'>Click + to expand the boxes and see the content</b>",
        position: "left"
      },
      {
        element: document.querySelector("#b_dmp_dmr"),
        intro: "Next step"
      },
      {
        element: document.querySelector("#div_dmp_calculation_options"),
        title: "DMPs Calculation Options",
        intro: "Before the DMPs calculation, is needed a model and you can select some options to generate it. Then, there are other options for the contrasts calculation. <br> For this example, we use the default options. <hr> <small style='color:steelblue;'> For more information you can consult <b>Help</b> section </small>"
      },
      { 
        intro: ""
      }
      
    ]});
    
  }
);




Shiny.addCustomMessageHandler("intro_continue2", 
  function(message){
    intro.onchange(function(targetElement) {
      if (intro._currentStep==0){
        intro.setOptions({'nextLabel': 'Finish Quality Control'});
      }
          if (intro._currentStep==1) {
              $('a[data-value=\"analysis\"]').trigger('click');
              intro.setOptions({'nextLabel': 'Go to Exploratory Analysis', "disableInteraction": "true"});
              document.querySelector(".introjs-prevbutton").style.display = "none";
          }
          if (intro._currentStep==2) {
              $('a[data-value=\"exploratory_analysis\"]').trigger('click');
              intro.setOptions({'nextLabel': 'Finish Exploratory Analysis'});
          }
          if (intro._currentStep==3) {
              $('a[data-value=\"analysis\"]').trigger('click');
              intro.setOptions({'nextLabel': 'Go to DMP/DMR', "disableInteraction": "true"});
          }
          if (intro._currentStep==4) {
              $('a[data-value=\"dmp_dmr\"]').trigger('click');
              intro.setOptions({'nextLabel': 'Calculate', "disableInteraction": "true"});
          }
          if (intro._currentStep==5) {
              $("#button_limma_calculatemodel").click();
              intro.exit();
          }
       }).start();
  }
);

/*
document.getElementById('startButton').onclick = function() {
       introJs().onchange(function(targetElement) {
          if (this._currentStep==0) {
             $('a[data-value=\"analysis\"]').removeClass('active');
             $('a[data-value=\"data\"]').addClass('active');
             $('a[data-value=\"data\"]').trigger('click');
          }
          if (this._currentStep==1) {
             $('a[data-value=\"data\"]').removeClass('active');
             $('a[data-value=\"analysis\"]').addClass('active');
             $('a[data-value=\"analysis\"]').trigger('click');
          }
       }).start();
};

intro.onafterchange(function(targetElement) {  
  if(intro._currentStep == 2) { // your disabled step 2 for example
    var original_onclick = $('.introjs-nextbutton').get(0).onclick;
    $('.introjs-nextbutton').addClass('introjs-disabled');
    $('.introjs-nextbutton').get(0).onclick = null;
    $('#searchbox').on('search:done', function() {
      $('.introjs-nextbutton').removeClass('introjs-disabled');
      $('.introjs-nextbutton').get(0).onclick = original_onclick;
    }
*/






Shiny.addCustomMessageHandler("intro_steps_continue3",

  function(message){
    intro.setOptions({
      exitOnOverlayClick: false,
      showBullets: false,
      disableInteraction: false,
      steps: [
      
      {
        element: document.querySelector("#div_dmp_table_options"),
        title: "DMPs Table and Options",
        intro: "In this section, you can see the DMPs table and some options that can be modified and updated. <br> For this example, we use the default values."
      },
      {
        element: document.querySelector("#div_dmp_plots"),
        title: "DMPs Plots",
        intro: "The DMPs results are represented in different plots that you can explore. <hr> <b style='color:steelblue;'>Click + to expand the boxes and see the content</b>",
        position: "left"
      },
      {
        element: document.querySelector("#div_dmr_calculation_options"),
        title: "DMRs Calculation Options",
        intro: "For the DMRs calculation, you can select some options. <br> For this example, we select only <i>promoters</i> and use the default options for the rest. <hr> <small style='color:steelblue;'> For more information you can consult <b>Help</b> section </small>"
      },
      {
        intro: ""
      }
    ]});
    
  }
);




Shiny.addCustomMessageHandler("intro_continue3", 
  function(message){
    intro.onchange(function(targetElement) {
      if (intro._currentStep==0){
        intro.setOptions({'nextLabel': 'Next'});
      }
          if (intro._currentStep==1) {
              $("#button_limma_tablecalc").click();
              intro.setOptions({'nextLabel': 'Finish DMPs', "disableInteraction": "false"});
          }
          if (intro._currentStep==2) {
              $('a[data-value=\"DMRs\"]').trigger('click');
              intro.setOptions({'nextLabel': 'Calculate DMR'});
              document.querySelector(".introjs-prevbutton").style.display = "none";
          }
          if (intro._currentStep==3) {
              $("#button_dmrs_calculate").click();
              intro.exit();
          }
       }).start();
  }
);





Shiny.addCustomMessageHandler("intro_steps_continue4",

  function(message){
    intro.setOptions({
      exitOnOverlayClick: false,
      showBullets: false,
      disableInteraction: false,
      steps: [
      
      {
        element: document.querySelector("#div_dmr_table_options"),
        title: "DMRs Table and Options",
        intro: "In this section, you can see the DMRs table and some options that can be modified and updated. <br> For this example, we use the default values."
      },
      {
        element: document.querySelector("#div_dmr_plots"),
        title: "DMRs Plots",
        intro: "The DMRs results are represented in different plots that you can explore. <hr> <b style='color:steelblue;'>Click + to expand the boxes and see the content</b>",
        position: "left"
      },
      {
        element: document.querySelector("#b_functional_enrichment"),
        intro: "Next step"
      },
      {
        element: document.querySelector("#functional_enrichment_plots"),
        title: "Functional Enrichment",
        intro: "There are some plots of different biological ontologies to visualize the functional enrichment. <hr> <b style='color:steelblue;'>Click + to expand the boxes and see the content</b>"
      },
      {
        element: document.querySelector("#b_survival"),
        intro: "Next step"
      },
      {
        element: document.querySelector("#ui_clinical_data"),
        title: "Load Clinical Data",
        intro: "Clinical data must be in a .csv file. <br> In this case, we are going to load the example data."
      },
      {
        element: document.querySelector("#div_clinical_options"),
        title: "Clinical data options",
        intro: "Before doing the survival analysis, you have to select the correct column of the clinical data table in each option."
      },
      {
        intro: ""
      }
    ]});
    
  }
);




Shiny.addCustomMessageHandler("intro_continue4", 
  function(message){
    intro.onchange(function(targetElement) {
      if (intro._currentStep==0){
        intro.setOptions({'nextLabel': 'Next'});
      }
          if (intro._currentStep==1) {
              $("#button_dmrs_tablecalc").click();
              intro.setOptions({'nextLabel': 'Finish DMRs'});
          }
          if (intro._currentStep==2) {
              $('a[data-value=\"analysis\"]').trigger('click');
              intro.setOptions({'nextLabel': 'Go to Functional Enrichment'});
              document.querySelector(".introjs-prevbutton").style.display = "none";
          }
          if (intro._currentStep==3) {
              $('a[data-value=\"functional_enrichment\"]').trigger('click');
              intro.setOptions({'nextLabel': 'Finish Functional Enrichment'});
              document.querySelector(".introjs-prevbutton").style.display = "block";
          }
          if (intro._currentStep==4) {
              $('a[data-value=\"analysis\"]').trigger('click');
              intro.setOptions({'nextLabel': 'Go to Survival'});
          }
          if (intro._currentStep==5) {
              $('a[data-value=\"survival\"]').trigger('click');
              $("#b_clinical_data").click();
              intro.setOptions({'nextLabel': 'Load Clinical Data'});
          }
          if (intro._currentStep==6) {
              //$("#b_clinical_data").click();
              intro.setOptions({'nextLabel': 'Continue to Survival'});
          }
          if (intro._currentStep==7) {
              $("#b_clinical_next").click();
              intro.exit();
          }
       }).start();
  }
);





Shiny.addCustomMessageHandler("intro_steps_continue5",

  function(message){
    intro.setOptions({
      exitOnOverlayClick: false,
      showBullets: false,
      disableInteraction: false,
      steps: [
      
      {
        intro: "next"
      },
      {
        element: document.querySelector("#div_clin_meth_options"),
        title: "Clinical and Methylation Options",
        intro: "In this section, you have to select the variables for survival. <br> For this example, <i>Sex</i> variable from clinical data and the methylation data computed before. <hr> <small style='color:steelblue;'><i>It is also possible to perform a survival analysis without methylation data</i></small>"
      },
      {
        element: document.querySelector("#div_survival_plots"),
        title: "Survival Plots",
        intro: "The Kaplan-Meier plot is obtained and other useful information. <hr> <b style='color:steelblue;'>Click + to expand the boxes and see the content</b>",
        position: "left"
      },
      {
        element: document.querySelector("#div_export"),
        title: "Download Report",
        intro: "Finally, you can download a report with the results obtained in HTML format.",
        position: "left"
      },
      {
        intro: ""
      }
    ]});
    
  }
);




Shiny.addCustomMessageHandler("intro_continue5", 
  function(message){
    intro.onchange(function(targetElement) {
      if (intro._currentStep==0){
        intro.setOptions({'nextLabel': 'Run Survival'});
      }
      if (intro._currentStep==1){
        intro.setOptions({'nextLabel': 'Run Survival'});
        document.querySelector(".introjs-prevbutton").style.display = "none";
      }
          if (intro._currentStep==2) {
              $("#b_run_survival").click();
              intro.setOptions({'nextLabel': 'Finish Survival', "disableInteraction": "false"});
              document.querySelector(".introjs-prevbutton").style.display = "block";
          }
          if (intro._currentStep==3) {
              $('a[data-value=\"export\"]').trigger('click');
              intro.setOptions({'nextLabel': 'Download Report'});
              document.querySelector(".introjs-prevbutton").style.display = "none";
          }
          if (intro._currentStep==4) {
            $("#complete_report_html").click();
              intro.exit();
          }
       }).start();
  }
);







Shiny.addCustomMessageHandler("intro_steps_continue6",

  function(message){
    intro.setOptions({
      exitOnOverlayClick: false,
      showBullets: false,
      disableInteraction: true,
      steps: [
      {
        title: "TOUR INTROJS FINISHED",
        intro: "end"
      }
    ]});
    
  }
);




Shiny.addCustomMessageHandler("intro_continue6", 
  function(message){
    intro.start();
  }
);