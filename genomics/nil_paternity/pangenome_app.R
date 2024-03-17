# Dependencies 
library(shiny)
library(dplyr)
library(stringr)

# Load the passport data 
passport <- read.csv("passport_results.csv")

# Load the accessions data 
accessions <- readRDS("accessions.rds")

# Read some functions 
get_passport <- function(id) {
  # Exceptions to syntax 'prefix_empty space_number'
  if (id %in% c("GCP 105", "WR 315")) {
    id = str_replace(id, " ", "")
  }
  
  if (id %in% c("CDC 512 51", "DCP 92 3")) { 
    id = str_replace_all(id, " ", "-")
    id = str_replace(id, "-", " ")
  }
  
  
  if (id %in% passport$Genotype) {
    print(passport[which(str_detect(passport$Genotype, id)), ])
  } else if (id %in% passport$Alternate.name) {
    print(passport[which(str_detect(passport$Alternate.name, id)), ])
  } else {
    print("Apparently, there are no records for that entry.")
  }
}

get_colSNP <- function(id) {
  
  # Exceptions to syntax 'prefix_empty space_number'
  if(id %in% c("GCP 105", "WR 315")) {
    id = str_replace(id, " ", "")
    print(which(id == accessions))
  } else {
    id = str_replace(id, " ", ".")
    if(length(which(id == accessions) > 0)) {
      print(which(id == accessions))
    } else {
      id = str_replace(id, "\\.", " ")
      
      id = passport[which(str_detect(passport$Alternate.name, id)),1] #Genotype name
      id = str_replace(id, " ", ".")
      print(which(id == accessions))
    }
  }
}



# Define UI for the application
ui <- fluidPage(
  h4("Our app extracts essential information from the entry passport. 
Additionally, the app obtains the the number of column required to obtain a comprehensive 
collection of Single Nucleotide Polymorphism (SNP) data associated with each entry. "),
  titlePanel("Passport Information"),
  
  sidebarLayout(
    sidebarPanel(
      uiOutput("intro"),
      textInput("id", "Enter ID", ""),
      actionButton("submit", "Submit")
    ),
    mainPanel(
      fluidRow(
        column(12, h4("Passport Information:")),
        column(12, tableOutput("passportTable")),
        column(12, verbatimTextOutput("snpOutput"))
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  output$intro <- renderUI({
    HTML("Enter your genotype using the following syntax : prefix_empty space_number.<br>
          For example, if you need 'JG62', type `JG 62`; if 'ICC7560', type 'ICC 7560'.<br>
          If the id has more than two elements, separate them with spaces. <br>
          F. ex., DC 512 51, DCP 92 3, or Phule G 12")
  })
  
  observeEvent(input$submit, {
    id <- input$id
    passport_info <- get_passport(id)
    output$passportTable <- renderTable({
      if (is.data.frame(passport_info)) {
        passport_info
      } else {
        data.frame(Message = passport_info)
      }
    })
    
    # Display SNP information
    snp_info <- get_colSNP(id)
    output$snpOutput <- renderText({
      paste("SNPs for genotype", id, "in 'Cultivated chickpea SNPs files' column:", snp_info)
    })
  })
}

# Run the application
shinyApp(ui, server)