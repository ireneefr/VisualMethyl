
var intro = introJs();

Shiny.addCustomMessageHandler("intro_steps",

  function(message){
    intro.setOptions({
      exitOnOverlayClick: false,
      showBullets: false,
      disableInteraction: true,
      steps: [
      
      {
        title: "TOUR INTROJS",
        intro: "This is ................."
      },
      {
        element: document.querySelector("#ui_input_data"),
        intro: "Load data"
      },
      {
        element: document.querySelector("#div_samples_table"),
        intro: "Samples table"
      },
      {
        element: document.querySelector("#div_select_options"),
        intro: "Select options",
        position: "right"
      },
      {
        element: document.querySelector("#b_qc"),
        intro: "Quality Control"
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
            
          }
          if(this._currentStep==2){
            $("#b_input_data").click();
            intro.setOptions({'nextLabel': 'Next', "disableInteraction": "false"});
            document.querySelector(".introjs-prevbutton").style.display = "block";
          }
          if(this._currentStep==3){
            intro.setOptions({'nextLabel': 'Continue to Analysis', "disableInteraction": "true"});
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
        intro: "Quality Control"
      },
      {
        element: document.querySelector("#div_norm_type"),
        intro: "Select type of normalization"//,
        //position: "right"
      },
      {
        element: document.querySelector("#div_norm_options"),
        intro: "Other options",
        position: "right"
      },
      {
        element: document.querySelector("#vertical"),
        intro: "Samples table"
      }
      
    ]});
    
  }
);




Shiny.addCustomMessageHandler("intro_continue1", 
  function(message){
    intro.onchange(function(targetElement) {
      if (intro._currentStep==0){
        intro.setOptions({'nextLabel': 'Go to Quality Control', "disableInteraction": "true"});
      }
          if (intro._currentStep==1) {
             $('#b_qc').click();
             intro.setOptions({'nextLabel': 'Next'});
             document.querySelector(".introjs-prevbutton").style.display = "none";
          }
          if(intro._currentStep==2){
            intro.setOptions({'nextLabel': 'Run Normalzation'});
            document.querySelector(".introjs-prevbutton").style.display = "block";
          }
          if (intro._currentStep==3){
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
      //disableInteraction: true,
      steps: [
      
      {
        element: document.querySelector("#vertical"),
        intro: "Quality Control Plots",
        position: "left"
      },
      {
        element: document.querySelector("#b_exploratory_analysis"),
        intro: "Next exploratory analysis"
      },
      {
        element: document.querySelector("#exploratory_plots"),
        intro: "Exploratory Analysis Plot",
        position: "left"
      },
      {
        element: document.querySelector("#b_dmp_dmr"),
        intro: "Next DMP/DMR"
      },
      {
        element: document.querySelector("#div_dmp_calculation_options"),
        intro: "Select options to generate the model and to calculate the contrasts"
      },
      { 
        intro: "Select contrast options"
      }
      
    ]});
    
  }
);




Shiny.addCustomMessageHandler("intro_continue2", 
  function(message){
    intro.onchange(function(targetElement) {
      if (intro._currentStep==0){
        intro.setOptions({'nextLabel': 'Finish Quality Control', "disableInteraction": "false"});
      }
          if (intro._currentStep==1) {
              $('a[data-value=\"analysis\"]').trigger('click');
              intro.setOptions({'nextLabel': 'Go to Exploratory Analysis', "disableInteraction": "true"});
          }
          if (intro._currentStep==2) {
              $("#b_exploratory_analysis").click();
              intro.setOptions({'nextLabel': 'Finish Exploratory Analysis', "disableInteraction": "false"});
          }
          if (intro._currentStep==3) {
              $('a[data-value=\"analysis\"]').trigger('click');
              intro.setOptions({'nextLabel': 'Go to DMP/DMR', "disableInteraction": "true"});
          }
          if (intro._currentStep==4) {
              $("#b_dmp_dmr").click();
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
      disableInteraction: true,
      steps: [
      
      {
        element: document.querySelector("#table_limma_difcpgs"),
        intro: "Table DMPs"
      },
      {
        element: document.querySelector("#div_dmp_options"),
        intro: "Options DMPs that can be updated"
      },
      {
        element: document.querySelector("#div_dmp_plots"),
        intro: "DMP Plot",
        position: "left"
      },
      {
        element: document.querySelector("#div_dmr_calculation_options"),
        intro: "Options to calculate DMRs"
      },
      {
        element: document.querySelector("#div_model_options"),
        intro: "Select options to generate the model"
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
              intro.setOptions({'nextLabel': 'Next'});
          }
          if (intro._currentStep==2) {
              $("#button_limma_tablecalc").click();
              intro.setOptions({'nextLabel': 'Finish DMPs', "disableInteraction": "false"});
          }
          if (intro._currentStep==3) {
              $('a[data-value=\"DMRs\"]').trigger('click');
              intro.setOptions({'nextLabel': 'Calculate DMR'});
          }
          if (intro._currentStep==4) {
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
      disableInteraction: true,
      steps: [
      
      {
        element: document.querySelector("#div_dmr_table"),
        intro: "Table DMRs"
      },
      {
        element: document.querySelector("#div_dmr_options"),
        intro: "Options DMRs that can be updated"
      },
      {
        element: document.querySelector("#div_dmr_plots"),
        intro: "DMR Plot",
        position: "left"
      },
      {
        element: document.querySelector("#b_functional_enrichment"),
        intro: "Next Functional Enrichment"
      },
      {
        element: document.querySelector("#div_functional_enrichment_plots"),
        intro: "Next Functional Enrichment"
      },
      {
        element: document.querySelector("#b_survival"),
        intro: "Next Survival"
      },
      {
        element: document.querySelector("#ui_clinical_data"),
        intro: "Load clinical data"
      },
      {
        element: document.querySelector("#div_clinical_options"),
        intro: "Clinical data options"
      },
      {
        element: document.querySelector("#div_clinical_options"),
        intro: "Clinical data options"
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
              intro.setOptions({'nextLabel': 'Next'});
          }
          if (intro._currentStep==2) {
              $("#button_dmrs_tablecalc").click();
              intro.setOptions({'nextLabel': 'Finish DMRs', "disableInteraction": "false"});
          }
          if (intro._currentStep==3) {
              $('a[data-value=\"analysis\"]').trigger('click');
              intro.setOptions({'nextLabel': 'Go to Functional Enrichment'});
          }
          if (intro._currentStep==4) {
              $("#b_functional_enrichment").click();
              intro.setOptions({'nextLabel': 'Finish Functional Enrichment', "disableInteraction": "false"});
          }
          if (intro._currentStep==5) {
              $('a[data-value=\"analysis\"]').trigger('click');
              intro.setOptions({'nextLabel': 'Go to Survival'});
          }
          if (intro._currentStep==6) {
              $("#b_survival").click();
              intro.setOptions({'nextLabel': 'Load Clinical Data'});
          }
          if (intro._currentStep==7) {
              $("#b_clinical_data").click();
              intro.setOptions({'nextLabel': 'Continue to Survival'});
          }
          if (intro._currentStep==8) {
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
      disableInteraction: true,
      steps: [
      
      {
        element: document.querySelector("#div_clin_meth_options"),
        intro: "Clinical and methylation options"
      },
      {
        element: document.querySelector("#div_survival_plots"),
        intro: "Survival plots"
      },
      {
        element: document.querySelector("#div_export"),
        intro: "Download Report",
        position: "left"
      },
      {
        element: document.querySelector("#b_functional_enrichment"),
        intro: "Next Functional Enrichment"
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
          if (intro._currentStep==1) {
              $("#b_run_survival").click();
              intro.setOptions({'nextLabel': 'Finish Survival', "disableInteraction": "false"});
          }
          if (intro._currentStep==2) {
              $('a[data-value=\"export\"]').trigger('click');
              intro.setOptions({'nextLabel': 'Download Report'});
          }
          if (intro._currentStep==3) {
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