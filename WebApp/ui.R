library(shiny)
library(ggplot2)
library(dr4pl)

### Define UI for app
ui = fluidPage(
  
  theme = "flatly",
  
  tags$head(
    tags$style(HTML("hr {border-top: 1px solid #4D4D4D;}"))
  ),
  
  #3 App title
  titlePanel("Serum Shift Factor Calculator"),
  
  ## Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    ## Sidebar panel for inputs
    sidebarPanel(
      
      # Header for dose response data
      h4("Input dose response data:"),
      h5("Data should be provided as a .csv file organized in tidy format without column names. The required columns are:"),
      p(HTML("<ol><li><strong>Drug Concentration</strong>: Numeric indicating drug concentration. May be log<sub>10</sub> transformed (see below).</li>
             <li><strong>Serum Added</strong>: Logical indicating whether serum was added.</li>
             <li><strong>Replicate</strong>: Numeric indicating experimental replicate.</li>
             <li><strong>Assay value</strong>: Numeric indicating assay value (e.g. relative viability normalized to vehicle control).</li>
             </ol>")),
      p(HTML("Exmaple data can be downloaded <a href='https://github.com/pritchardlabatpsu/DrugRepurposing/blob/main/WebAppSampleData.csv'><b>here</b></a>.")),
      
      # Input: dose response data
      fileInput(inputId = "dr_data",
                label = "Upload dose response data:",
                accept = ".csv"),
      
      # Input: Drug concentration unit
      selectInput(inputId = "conc_unit",
                  label = "Drug concentration units:",
                  choices = c("M" = "M",
                              "mM" = "mM",
                              "uM" = "uM",
                              "nM" = "nM",
                              "pM" = "pM")),
      
      # Input: Drug concentration format
      radioButtons(inputId = "conc_format",
                   label = "Drug concentration format:",
                   c("log10 Transformed" = "log10trns",
                     "Untransformed" = "untrns")
      ),
      
      # Header for Cave data
      hr(),
      h4(HTML("<em>**Optional**")),
      p(HTML("<em>Enter C<sub>ave</sub> value to be shift-corrected. C<sub>ave</sub> parameter may be measured or theoretical.</em>")),
      
      # Input: Cave values
      textInput("Cave",HTML("Enter C<sub>ave</sub> :")),
      
      # Action buttom: Submit
      hr(),
      actionButton(
        inputId = "submit_button",
        label = "Calculate"
      ),
      
      width = 4
    ),
    
    ## Main panel for displaying outputs
    mainPanel(
      
      tabsetPanel(
        type="tabs",
        tabPanel("Tutorial",
                 p(),
                 h3("Introduction"),
                 p(HTML("This app is used to calculate serum shift factors for <em>in vitro</em> studies of drug activity. 
                   Before using the app, conduct a dose-response assay in the presence and absence of human serum proteins, as in Liu et al.
                   Normalize the readout of the assay to the vehicle control for both conditions, and format the data in tidy format
                   (i.e. each row contains information for a single observation.) A description of the input file format is provided in the side panel,
                   and sample data can be downloaded <a href='https://github.com/pritchardlabatpsu/DrugRepurposing/blob/main/WebAppSampleData.csv'><b>here</b></a>.")),
                 h3("App Details"),
                 p("The app fits a 4 parameter logistic curve for each of the two conditions (with and without human serum proteins). 
                   The addition of serum should shift the dose response curve rightwards, as below."),
                 img(src = 'SerumShiftSchematicWebApp.png',width='600px',height='400px'),
                 p(HTML("The app then calculates the fold-change in IC<sub>50</sub> values. 
                   This ratio - termed the <b>serum shift factor</b> - is a single parameter that reflects the affinity of a drug for its target and for off-target serum proteins.
                   The serum shift factor can be used to determine the effective drug concentration for subsequent experiments conducted in the absence human serum proteins.")),
                 p(HTML("In the case where the pharmacokinetic properties of the drug are known (or can be reasonably estimated), 
                   the serum shift factor can be used to correct off-target binding effects to provide the equivalent <em>effective</em> pharmacokinetic parameters <em>in vitro</em>.
                   For example, a drug with a measured C<sub>ave</sub> of 500 nM and a serum shift factor &sigma; of 5 would have an effective C<sub>ave</sub> of 100 nM.")),
                 img(src="EffectiveCaveSchematicWebApp.png",width='600px',height='400px'),
                 p("A table of approved tyrosine kinase inhibitors, their pharmacokinetic parameters, 
                   and their serum shift factors as reported in Liu et al. is provided below:"),
                 img(src="TKITableWebApp.png",width='600px',height='200px'),
                 h3("App Execution"),
                 p(HTML("To calculate the serum shift factor for your drug of interest, upload your dose response data in the format indicated in th side panel. Then, note the drug concentration units and format used. 
                   Enter the C<sub>ave</sub> if known. Finally, click 'Calculate' and navigate to the 'Results' tab to visualize the curve fitting and access the calculated outputs."))
        ),
        tabPanel("Results",
                 h1("Dose Response Curve"),
                 
                 # Display plot output
                 plotOutput("dr_plot"),
                 
                 h1("Output"),
                 
                 # Display results output
                 tableOutput("shift_table"))
      )
      
      
    )
  )
)