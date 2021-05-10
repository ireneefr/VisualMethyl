
var intro = introJs();

Shiny.addCustomMessageHandler("intro_steps",

  function(message){
    intro.setOptions({
      showBullets: false,
      steps: [
      
      {
        title: "TOUR INTROJS",
        intro: "This is ................."
      },
      {
        element: document.querySelector("#div_upload_data"),
        intro: "Upload methylation data"
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
          
          if(intro._currentStep==4){
            intro.setOptions({'nextLabel': 'Continue to Analysis'});
            /*
            if ($("#button_input_next").prop("disabled")){
              document.querySelector(".introjs-nextbutton").style.display = "none";
            }
            if ($("#button_input_next").prop("enabled")){
              document.querySelector(".introjs-nextbutton").style.display = "block";
            }*/
          }

          if (intro._currentStep==5) {
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
    intro.setOptions({steps: [
      
      {
        element: document.querySelector("#div_btn_qc"),
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
        intro.setOptions({'nextLabel': 'Go to Quality Control'});
      }
          if (intro._currentStep==1) {
             $('#b_qc').click();
             intro.setOptions({'nextLabel': 'Next'});
          }
          if(intro._currentStep==2){
            intro.setOptions({'nextLabel': 'Run Normalzation'});
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
    intro.setOptions({steps: [
      
      {
        element: document.querySelector("#vertical"),
        intro: "Samples table",
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
        element: document.querySelector("#vertical"),
        intro: "Samples table"
      },
      {
        element: document.querySelector("#div_select_options"),
        intro: "Select options"
      },
      {
        element: document.querySelector("#b_qc"),
        intro: "Quality Control"
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
              intro.setOptions({'nextLabel': 'Go to Exploratory Analysis'});
          }
          if (intro._currentStep==2) {
              $("#b_exploratory_analysis").click();
              intro.setOptions({'nextLabel': 'Finish Exploratory Analysis'});
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

