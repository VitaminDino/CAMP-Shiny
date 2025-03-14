library(shiny)
library(shinyWidgets)
library(seqinr)
library(dplyr)
library(VennDiagram)
library(ggplot2)
library(shinyjs)
library(shinycssloaders)
library(openxlsx)
library(DT)
library(gridExtra)

# Define amino acid monoisotopic masses
AMINOACID_MASS <- c(
  A = 71.03711, R = 156.10111, N = 114.04293, D = 115.02694, C = 103.00919,
  E = 129.04259, Q = 128.05858, G = 57.02146, H = 137.05891, I = 113.08406,
  L = 113.08406, K = 128.09496, M = 131.04049, F = 147.06841, P = 97.05276,
  S = 87.03203, T = 101.04768, W = 186.07931, Y = 163.06333, V = 99.06841
)

# Function to calculate isoelectric point (pI)
calculate_pI <- function(sequence) {
  aa_pKa <- c(
    A = 2.35, R = 9.09, N = 2.02, D = 2.09, C = 1.71, E = 2.19, Q = 2.17,
    G = 2.34, H = 2.09, I = 2.32, L = 2.36, K = 2.18, M = 2.28, F = 2.58,
    P = 1.99, S = 2.21, T = 2.09, W = 2.83, Y = 2.20, V = 2.29
  )
  return(mean(aa_pKa[unlist(strsplit(sequence, ""))], na.rm = TRUE))
}

# Function to analyze peptides with adjustable scoring parameters
analyze_peptides <- function(fasta_files, fasta_names,
                             apply_charge, charge_range,
                             apply_length, length_range,
                             apply_hydrophobic, hydrophobic_range,
                             apply_pI, pI_range,
                             apply_cationic, cationic_min,
                             apply_aromatic, aromatic_min) {
  peptide_data <- list()
  for (i in seq_along(fasta_files)) {
    file <- fasta_files[i]
    file_name <- fasta_names[i]
    sequences <- read.fasta(file, seqtype = "AA", as.string = TRUE)
    for (seq in sequences) {
      peptide <- toupper(unlist(strsplit(seq[[1]], "")))
      len <- length(peptide)
      hydrophobic_residues <- sum(peptide %in% c("A", "V", "I", "L", "M", "F", "W", "Y"))
      percent_hydrophobic <- (hydrophobic_residues / len) * 100
      charge_state <- sum(peptide %in% c("R", "K", "H")) - sum(peptide %in% c("D", "E"))
      mz <- sum(sapply(peptide, function(aa) AMINOACID_MASS[aa])) + 1.0073
      pI <- calculate_pI(paste(peptide, collapse = ""))
      cationic_residues <- sum(peptide %in% c("R", "K")) / len * 100
      aromatic_residues <- sum(peptide %in% c("W", "F", "Y")) / len * 100
      
      # Evaluate each criterion based on user settings
      charge_condition <- if (apply_charge) {
        charge_state >= charge_range[1] & charge_state <= charge_range[2]
      } else { TRUE }
      
      length_condition <- if (apply_length) {
        len >= length_range[1] & len <= length_range[2]
      } else { TRUE }
      
      hydrophobic_condition <- if (apply_hydrophobic) {
        percent_hydrophobic >= hydrophobic_range[1] & percent_hydrophobic <= hydrophobic_range[2]
      } else { TRUE }
      
      pI_condition <- if (apply_pI) {
        pI >= pI_range[1] & pI <= pI_range[2]
      } else { TRUE }
      
      cationic_condition <- if (apply_cationic) {
        cationic_residues > cationic_min
      } else { TRUE }
      
      aromatic_condition <- if (apply_aromatic) {
        aromatic_residues > aromatic_min
      } else { TRUE }
      
      # 1 point per criterion met
      score <- sum(charge_condition, length_condition, hydrophobic_condition, 
                   pI_condition, cationic_condition, aromatic_condition)
      
      peptide_data <- append(peptide_data, list(data.frame(
        Sequence = paste(peptide, collapse = ""),
        Length = len,
        Hydrophobicity = percent_hydrophobic,
        Charge = charge_state,
        m_z = mz,
        pI = pI,
        Cationic_Residue_Percentage = cationic_residues,
        Aromatic_Residue_Percentage = aromatic_residues,
        Score = score,
        Source = basename(file_name)
      )))
    }
  }
  return(do.call(rbind, peptide_data))
}

ui <- fluidPage(
  useShinyjs(),
  titlePanel("CAMP Peptide Analyzer"),
  sidebarLayout(
    sidebarPanel(
      fileInput("fasta", "Upload up to 3 FASTA Files", multiple = TRUE, accept = ".fasta"),
      actionButton("analyze", "Analyze"),
      actionButton("reset", "Reset"),
      br(), br(),
      h3("Adjust Scoring Parameters"),
      checkboxInput("apply_charge", "Apply Charge State Criterion", value = TRUE),
      sliderInput("charge_range", "Charge State Range (Recommended: 2-8)", 
                  min = -10, max = 20, value = c(2, 8)),
      checkboxInput("apply_length", "Apply Peptide Length Criterion", value = TRUE),
      sliderInput("length_range", "Peptide Length Range (Recommended: 10-50)", 
                  min = 5, max = 100, value = c(10, 50)),
      checkboxInput("apply_hydrophobic", "Apply Hydrophobicity Criterion", value = TRUE),
      sliderInput("hydrophobic_range", "Hydrophobicity % Range (Recommended: 30-50)", 
                  min = 0, max = 100, value = c(30, 50)),
      checkboxInput("apply_pI", "Apply pI Criterion", value = TRUE),
      sliderInput("pI_range", "Isoelectric Point (pI) Range (Recommended: 6-10)", 
                  min = 0, max = 14, value = c(6, 10)),
      checkboxInput("apply_cationic", "Apply Cationic Residue Criterion", value = TRUE),
      sliderInput("cationic_min", "Minimum Cationic Residue Percentage (Recommended: 10)", 
                  min = 0, max = 100, value = 10),
      checkboxInput("apply_aromatic", "Apply Aromatic Residue Criterion", value = TRUE),
      sliderInput("aromatic_min", "Minimum Aromatic Residue Percentage (Recommended: 5)", 
                  min = 0, max = 100, value = 5),
      br(),
      radioButtons("region_filter", "Filter by Region",
                   choices = list(
                     "All peptides" = "all",
                     "Unique to File 1" = "unique1",
                     "Unique to File 2" = "unique2",
                     "Unique to File 3" = "unique3",
                     "Overlap File 1 & 2" = "overlap12",
                     "Overlap File 1 & 3" = "overlap13",
                     "Overlap File 2 & 3" = "overlap23",
                     "Overlap all 3" = "all3"
                   ),
                   selected = "all"
      ),
      br(),
      downloadButton("downloadResults", "Download Filtered Results"),
      br(), br(),
      downloadButton("downloadVenn", "Download Venn Diagram")
    ),
    mainPanel(
      h3("How to Use the CAMP Peptide Analyzer"),
      p("1. Upload up to three FASTA files containing peptide sequences."),
      p("2. Click 'Analyze' to process the peptides and generate results."),
      p("3. Use the 'Filter by Region' options to view peptides in specific overlapping regions."),
      p("4. The table below displays the filtered results."),
      p("5. Download the filtered results or the Venn diagram using the buttons."),
      br(),
      h3("Scoring Criteria"),
      p("Each peptide is given a score based on the following criteria (1 point per condition met):"),
      tags$ul(
        tags$li("Charge State (z) between selected range (Recommended: 2-8)"),
        tags$li("Peptide Length between selected range (Recommended: 10-50 residues)"),
        tags$li("Percent Hydrophobic Residues between selected range (Recommended: 30-50%)"),
        tags$li("Isoelectric Point (pI) between selected range (Recommended: 6-10)"),
        tags$li("Cationic Residue Content (R & K) greater than selected minimum (Recommended: >10%)"),
        tags$li("Aromatic Residue Content (W, F, Y) greater than selected minimum (Recommended: >5%)")
      ),
      br(),
      DT::dataTableOutput("results"),
      br(),
      plotOutput("vennPlot")
    )
  )
)

server <- function(input, output, session) {
  # Reactive value to store all results
  results <- reactiveVal(data.frame())
  
  # Process FASTA files when "Analyze" is clicked using the current scoring settings
  observeEvent(input$analyze, {
    req(input$fasta)
    validate(need(length(input$fasta$datapath) > 0, "Please upload at least one FASTA file."))
    peptides <- analyze_peptides(
      input$fasta$datapath,
      input$fasta$name,
      input$apply_charge,
      input$charge_range,
      input$apply_length,
      input$length_range,
      input$apply_hydrophobic,
      input$hydrophobic_range,
      input$apply_pI,
      input$pI_range,
      input$apply_cationic,
      input$cationic_min,
      input$apply_aromatic,
      input$aromatic_min
    )
    if (nrow(peptides) == 0) {
      showNotification("No peptides found in the uploaded files.", type = "warning")
    } else {
      results(peptides)
    }
  })
  
  # Reset inputs and clear results when "Reset" is clicked
  observeEvent(input$reset, {
    reset("fasta")
    results(data.frame())
    showNotification("Reset complete. Upload new files to analyze.", type = "message")
  })
  
  # Filter results based on region selections
  filtered_results <- reactive({
    req(results())
    all_data <- results()
    if (input$region_filter == "all") return(all_data)
    
    sources <- unique(all_data$Source)
    if (length(sources) == 2) {
      if (input$region_filter == "unique1") {
        unique1 <- setdiff(all_data$Sequence[all_data$Source == sources[1]],
                           all_data$Sequence[all_data$Source == sources[2]])
        return(all_data[all_data$Sequence %in% unique1, ])
      } else if (input$region_filter == "unique2") {
        unique2 <- setdiff(all_data$Sequence[all_data$Source == sources[2]],
                           all_data$Sequence[all_data$Source == sources[1]])
        return(all_data[all_data$Sequence %in% unique2, ])
      } else if (input$region_filter == "overlap12") {
        common <- intersect(all_data$Sequence[all_data$Source == sources[1]],
                            all_data$Sequence[all_data$Source == sources[2]])
        return(all_data[all_data$Sequence %in% common, ])
      } else {
        return(all_data)
      }
    }
    
    if (length(sources) >= 3) {
      if (input$region_filter == "unique1") {
        unique1 <- setdiff(all_data$Sequence[all_data$Source == sources[1]],
                           union(all_data$Sequence[all_data$Source == sources[2]],
                                 all_data$Sequence[all_data$Source == sources[3]]))
        return(all_data[all_data$Sequence %in% unique1, ])
      } else if (input$region_filter == "unique2") {
        unique2 <- setdiff(all_data$Sequence[all_data$Source == sources[2]],
                           union(all_data$Sequence[all_data$Source == sources[1]],
                                 all_data$Sequence[all_data$Source == sources[3]]))
        return(all_data[all_data$Sequence %in% unique2, ])
      } else if (input$region_filter == "unique3") {
        unique3 <- setdiff(all_data$Sequence[all_data$Source == sources[3]],
                           union(all_data$Sequence[all_data$Source == sources[1]],
                                 all_data$Sequence[all_data$Source == sources[2]]))
        return(all_data[all_data$Sequence %in% unique3, ])
      } else if (input$region_filter == "overlap12") {
        common12 <- intersect(all_data$Sequence[all_data$Source == sources[1]],
                              all_data$Sequence[all_data$Source == sources[2]])
        common12 <- setdiff(common12, all_data$Sequence[all_data$Source == sources[3]])
        return(all_data[all_data$Sequence %in% common12, ])
      } else if (input$region_filter == "overlap13") {
        common13 <- intersect(all_data$Sequence[all_data$Source == sources[1]],
                              all_data$Sequence[all_data$Source == sources[3]])
        common13 <- setdiff(common13, all_data$Sequence[all_data$Source == sources[2]])
        return(all_data[all_data$Sequence %in% common13, ])
      } else if (input$region_filter == "overlap23") {
        common23 <- intersect(all_data$Sequence[all_data$Source == sources[2]],
                              all_data$Sequence[all_data$Source == sources[3]])
        common23 <- setdiff(common23, all_data$Sequence[all_data$Source == sources[1]])
        return(all_data[all_data$Sequence %in% common23, ])
      } else if (input$region_filter == "all3") {
        common_all <- Reduce(intersect, lapply(sources, function(s) 
          all_data$Sequence[all_data$Source == s]
        ))
        return(all_data[all_data$Sequence %in% common_all, ])
      } else {
        return(all_data)
      }
    }
    return(all_data)
  })
  
  output$results <- DT::renderDataTable({
    DT::datatable(filtered_results(), options = list(pageLength = 10, autoWidth = TRUE, filter = "top"))
  })
  
  # Render the Venn diagram based on the number of uploaded files
  output$vennPlot <- renderPlot({
    req(results())
    sources <- unique(results()$Source)
    if (length(sources) < 2) {
      showNotification("At least 2 files are needed for a Venn diagram", type = "warning")
      return(NULL)
    }
    venn_data <- lapply(sources, function(src) unique(results()$Sequence[results()$Source == src]))
    names(venn_data) <- sources
    
    if (length(sources) == 2) {
      venn.plot <- draw.pairwise.venn(
        area1 = length(venn_data[[1]]),
        area2 = length(venn_data[[2]]),
        cross.area = length(intersect(venn_data[[1]], venn_data[[2]])),
        category = names(venn_data),
        fill = c("red", "blue"),
        alpha = 0.5,
        cex = 1.5,
        cat.cex = 1.5,
        cat.col = c("red", "blue")
      )
      grid.draw(venn.plot)
    } else {
      venn.plot <- draw.triple.venn(
        area1 = length(venn_data[[1]]),
        area2 = length(venn_data[[2]]),
        area3 = length(venn_data[[3]]),
        n12 = length(intersect(venn_data[[1]], venn_data[[2]])),
        n13 = length(intersect(venn_data[[1]], venn_data[[3]])),
        n23 = length(intersect(venn_data[[2]], venn_data[[3]])),
        n123 = length(Reduce(intersect, venn_data)),
        category = names(venn_data),
        fill = c("red", "blue", "green"),
        alpha = 0.5,
        cex = 1.5,
        cat.cex = 1.5,
        cat.col = c("red", "blue", "green")
      )
      grid.draw(venn.plot)
    }
  })
  
  # Download filtered results as an Excel file
  output$downloadResults <- downloadHandler(
    filename = function() {
      paste("peptide_results_", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      req(filtered_results())
      write.xlsx(filtered_results(), file)
    }
  )
  
  # Download the Venn diagram as a PNG image
  output$downloadVenn <- downloadHandler(
    filename = function() {
      paste("venn_diagram_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      png(file, width = 800, height = 600)
      sources <- unique(results()$Source)
      if (length(sources) < 2) {
        dev.off()
        return(NULL)
      }
      venn_data <- lapply(sources, function(src) unique(results()$Sequence[results()$Source == src]))
      names(venn_data) <- sources
      
      if (length(sources) == 2) {
        venn.plot <- draw.pairwise.venn(
          area1 = length(venn_data[[1]]),
          area2 = length(venn_data[[2]]),
          cross.area = length(intersect(venn_data[[1]], venn_data[[2]])),
          category = names(venn_data),
          fill = c("red", "blue"),
          alpha = 0.5,
          cex = 1.5,
          cat.cex = 1.5,
          cat.col = c("red", "blue")
        )
        grid.draw(venn.plot)
      } else {
        venn.plot <- draw.triple.venn(
          area1 = length(venn_data[[1]]),
          area2 = length(venn_data[[2]]),
          area3 = length(venn_data[[3]]),
          n12 = length(intersect(venn_data[[1]], venn_data[[2]])),
          n13 = length(intersect(venn_data[[1]], venn_data[[3]])),
          n23 = length(intersect(venn_data[[2]], venn_data[[3]])),
          n123 = length(Reduce(intersect, venn_data)),
          category = names(venn_data),
          fill = c("red", "blue", "green"),
          alpha = 0.5,
          cex = 1.5,
          cat.cex = 1.5,
          cat.col = c("red", "blue", "green")
        )
        grid.draw(venn.plot)
      }
      dev.off()
    }
  )
}

shinyApp(ui, server)
