# Dataset Download ####

output$downloadDatasets1 <- downloadHandler(
    filename <- function() {
        paste("data-sets", "zip", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/data-sets.zip", file)
    },
    contentType = "application/zip"
)

output$downloadDatasets2 <- downloadHandler(
    filename <- function() {
        paste("Albatross_data", "zip", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/Albatross_data.zip", file)
    },
    contentType = "application/zip"
)

# Manual Download ####

output$downloadManual1 <- downloadHandler(
    filename <- function() {
        paste("Manual", "pdf", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/Manual.pdf", file)
    },
    contentType = "application/zip"
)

output$downloadManual2 <- downloadHandler(
    filename <- function() {
        paste("Manual(Minimal)", "pdf", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/Manual(Minimal).pdf", file)
    },
    contentType = "application/zip"
)

output$downloadManual3 <- downloadHandler(
    filename <- function() {
        paste("Manual(Korean)", "pdf", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/Manual(Korean).pdf", file)
    },
    contentType = "application/zip"
)

# Lecture Download ####

output$downloadLectureNote <- downloadHandler(
    filename <- function() {
        paste("Lecture Note", "zip", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/Lecture Note.zip", file)
    },
    contentType = "application/zip"
)

output$downloadLecture0 <- downloadHandler(
    filename <- function() {
        paste("[Chapter0] Introduction", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Chapter0] Introduction.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture1 <- downloadHandler(
    filename <- function() {
        paste("[Chapter1] Linear Models", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Chapter1] Linear Models.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture2 <- downloadHandler(
    filename <- function() {
        paste("[Chapter2] Generalized Linear Models", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Chapter2] Generalized Linear Models.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture3 <- downloadHandler(
    filename <- function() {
        paste("[Chapter3] Inference for Models with Unobservables", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Chapter3] Inference for Models with Unobservables.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture4 <- downloadHandler(
    filename <- function() {
        paste("[Chapter4] HGLMs; from method to algorithm", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Chapter4] HGLMs; from method to algorithm.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture5 <- downloadHandler(
    filename <- function() {
        paste("[Chapter5] HGLMs; modeling", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Chapter5] HGLMs; modeling.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture6 <- downloadHandler(
    filename <- function() {
        paste("[Chapter6] DHGLMs", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Chapter6] DHGLMs.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture7 <- downloadHandler(
    filename <- function() {
        paste("[Chapter7] MDHGLMs", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Chapter7] MDHGLMs.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture8 <- downloadHandler(
    filename <- function() {
        paste("[Chapter8] Survival Analysis", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Chapter8] Survival Analysis.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture9 <- downloadHandler(
    filename <- function() {
        paste("[Chapter9] Joint Models", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Chapter9] Joint Models.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture10 <- downloadHandler(
    filename <- function() {
        paste("[Chapter10] Further Topics", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Chapter10] Further Topics.mp4", file)
    },
    contentType = "application/zip"
)


# Lecture (English) Download ####

output$downloadLectureNote_eng <- downloadHandler(
    filename <- function() {
        paste("Lecture Note", "zip", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/Lecture Note.zip", file)
    },
    contentType = "application/zip"
)

output$downloadLecture0_eng <- downloadHandler(
    filename <- function() {
        paste("[Eng_Chapter0] Introduction", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Eng_Chapter0] Introduction.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture1_eng <- downloadHandler(
    filename <- function() {
        paste("[Eng_Chapter1] Linear Models", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Eng_Chapter1] Linear Models.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture2_eng <- downloadHandler(
    filename <- function() {
        paste("[Eng_Chapter2] Generalized Linear Models", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Eng_Chapter2] Generalized Linear Models.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture3_eng <- downloadHandler(
    filename <- function() {
        paste("[Eng_Chapter3] Inference for Models with Unobservables", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Eng_Chapter3] Inference for Models with Unobservables.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture4_eng <- downloadHandler(
    filename <- function() {
        paste("[Eng_Chapter4-5] HGLMs", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Eng_Chapter4-5] HGLMs.mp4", file)
    },
    contentType = "application/zip"
)

#output$downloadLecture5_eng <- downloadHandler(
#    filename <- function() {
#        paste("[Eng_Chapter5] HGLMs; modeling", "mp4", sep = ".")
#    },
#    content <- function(file) {
#        file.copy("/srv/shiny-server/files/[Chapter5] HGLMs; modeling.mp4", file)
#    },
#    contentType = "application/zip"
#)

output$downloadLecture6_eng <- downloadHandler(
    filename <- function() {
        paste("[Eng_Chapter6] DHGLMs", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Eng_Chapter6] DHGLMs.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture7_eng <- downloadHandler(
    filename <- function() {
        paste("[Eng_Chapter7] MDHGLMs", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Eng_Chapter7] MDHGLMs.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture8_eng <- downloadHandler(
    filename <- function() {
        paste("[Eng_Chapter8] Survival Analysis", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Eng_Chapter8] Survival Analysis.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture9_eng <- downloadHandler(
    filename <- function() {
        paste("[Eng_Chapter9] Joint Models", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Eng_Chapter9] Joint Models.mp4", file)
    },
    contentType = "application/zip"
)

output$downloadLecture10_eng <- downloadHandler(
    filename <- function() {
        paste("[Eng_Chapter10] Further Topics", "mp4", sep = ".")
    },
    content <- function(file) {
        file.copy("/srv/shiny-server/files/[Eng_Chapter10] Further Topics.mp4", file)
    },
    contentType = "application/zip"
)