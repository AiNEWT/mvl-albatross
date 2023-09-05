navbarMenu(
    "Materials",
    tabPanel(
        "Data-sets & Manual",
        h5("Download Dataset & Manual"),
        downloadButton("downloadDatasets1", label = "Download Data-sets"),
        br(), br(),
        downloadButton("downloadDatasets2", label = "Download Albatross-data"),
        br(), br(),
        downloadButton("downloadManual1", label = "Download Manual"),
        br(), br(),
        downloadButton("downloadManual2", label = "Download Manual(Minimal)"),
        br(), br(),
        downloadButton("downloadManual3", label = "Download Manual(Korean)")
    ),
    tabPanel(
        "Lecture (Korean)",
        h5("Download Lecture"),
        downloadButton("downloadLectureNote", label = "Lecture Note"),
        br(), br(),
        downloadButton("downloadLecture0", label = "[Chapter0] Introduction"),
        br(), br(),
        downloadButton("downloadLecture1", label = "[Chapter1] Linear Models"),
        br(), br(),
        downloadButton("downloadLecture2", label = "[Chapter2] Generalized Linear Models"),
        br(), br(),
        downloadButton("downloadLecture3", label = "[Chapter3] Inference for Models with Unobservables"),
        br(), br(),
        downloadButton("downloadLecture4", label = "[Chapter4] HGLMs; from Method to algorithm"),
        br(), br(),
        downloadButton("downloadLecture5", label = "[Chapter5] HGLMs; modeling"),
        br(), br(),
        downloadButton("downloadLecture6", label = "[Chapter6] DHGLMs"),
        br(), br(),
        downloadButton("downloadLecture7", label = "[Chapter7] MDHGLMs"),
        br(), br(),
        downloadButton("downloadLecture8", label = "[Chapter8] Survival Analysis"),
        br(), br(),
        downloadButton("downloadLecture9", label = "[Chapter9] Joint Models"),
        br(), br(),
        downloadButton("downloadLecture10", label = "[Chapter10] Further Topics"),
        br(), br()
    ),
    tabPanel(
        "Lecture (English)",
        h5("Download Lecture"),
        downloadButton("downloadLectureNote_eng", label = "Lecture Note"),
        br(), br(),
        downloadButton("downloadLecture0_eng", label = "[Chapter0] Introduction"),
        br(), br(),
        downloadButton("downloadLecture1_eng", label = "[Chapter1] Linear Models"),
        br(), br(),
        downloadButton("downloadLecture2_eng", label = "[Chapter2] Generalized Linear Models"),
        br(), br(),
        downloadButton("downloadLecture3_eng", label = "[Chapter3] Inference for Models with Unobservables"),
        br(), br(),
        downloadButton("downloadLecture4_eng", label = "[Chapter4-5] HGLMs"),
        br(), br(),
        #downloadButton("downloadLecture5_eng", label = "[Chapter5] HGLMs; modeling"),
        #br(), br(),
        downloadButton("downloadLecture6_eng", label = "[Chapter6] DHGLMs"),
        br(), br(),
        downloadButton("downloadLecture7_eng", label = "[Chapter7] MDHGLMs"),
        br(), br(),
        downloadButton("downloadLecture8_eng", label = "[Chapter8] Survival Analysis"),
        br(), br(),
        downloadButton("downloadLecture9_eng", label = "[Chapter9] Joint Models"),
        br(), br(),
        downloadButton("downloadLecture10_eng", label = "[Chapter10] Further Topics"),
        br(), br()
    ),
    tabPanel(
        "Version",
        h3("Version 1.2.0, 25/Feb/2022")
    ),
    tabPanel(
        "Albatross Forum",
        tags$ul(
            tags$li(
                tags$a(
                    href = "https://cms.pknu.ac.kr/albatross/view.do?no=11663", 
                    "Albatross Forum"
                )
            ),
       ),
       br(),
       hr()
    )

)