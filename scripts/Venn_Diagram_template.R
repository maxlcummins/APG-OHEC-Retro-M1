# Install required packages if not already installed
# devtools::install_github("yanlinlin82/ggvenn")
# devtools::install_github("gaospecial/ggVennDiagram")
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
        devtools::install_github("gaospecial/ggVennDiagram")
}
if (!requireNamespace("shiny", quietly = TRUE)) {
        install.packages("shiny")
}

library(ggVennDiagram)
library(shiny)

ui <- fluidPage(
        fluidRow(
                column(4, 
                       #sliderInput("SNP_threshold", "SNP Threshold", min = 1, max = 200, value = 100, width = "80%"),
                       numericInput("SNP_specific", "Specific SNP Value", value = 100, width = "80%"),
                       #sliderInput("allelic_dist_threshold", "Allelic Dist Threshold", min = 1, max = 100, value = 40, width = "80%"),
                       numericInput("allelic_specific", "Specific Allelic Dist Value", value = 40, width = "80%")
                ),
                column(4, plotOutput("venn_plot", height = "800px", width = "800px"))
        ),
        style = "margin-top: 20px;"  # Adjust the margin-top value as needed
)

server <- function(input, output, session) {
        observe({
                SNP_threshold <- input$SNP_specific
                allelic_dist_threshold <- input$allelic_specific
                
                dists_venn <- distsnew %>%
                        filter(ST.x %in% ST_list) %>% filter(Revised_Source_Niche.x != Revised_Source_Niche.y)
                
                SNP_overlaps <- dists_venn %>%
                        filter(SNPs <= SNP_threshold) %>%
                        mutate(name_combo = paste0(`Sample 1`, `Sample 2`, collapse = '-'))%>%
                        pull(name_combo)
                
                allelic_overlaps <- dists_venn %>%
                        filter(value <= allelic_dist_threshold) %>%
                        mutate(name_combo = paste0(`Sample 1`, `Sample 2`, collapse = '-')) %>%
                        pull(name_combo)
                
                x <- list(SNPs = SNP_overlaps, cgMLST = allelic_overlaps)#, Same_source = same_source, Different_Source = different_source)
                
                output$venn_plot <- renderPlot({
                        ggVennDiagram(x, label_alpha = 0, gridsize = 10, label_size = 12) + ggplot2::scale_fill_gradient(low = "white", high = "grey50")
                })
        })
}

shinyApp(ui, server)
