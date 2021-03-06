library(shiny)
library(TmCalculator)
library(stringr)
library(magrittr)
library(data.table)

ui <- fluidPage(
    titlePanel("DSN Probes"),
    
    sidebarLayout(
        sidebarPanel(
            textInput("dsnTaskTxt", "Task",
                      "AGATATTCCAACGTGCAAGTGGCTGGACCAGTG(G>A)ACAGAACTAGCTCAAAGGTATGTCCTAAATTAAATATAAGT"),
            sliderInput("probeMinMaxLength", "Probes Length", 
                        5, 40, c(15, 30)),
            sliderInput("probeMinMaxTm", "Probes Tm", 
                        30, 70, c(55, 65)),
            sliderInput("probeMinFlank", "Probe Min flanking Nucleotides", 
                        1, 10, 3),
            sliderInput("probeMaxOverlapTm", "Probes Overlap Max Tm", 
                        1, 60, 36),
            sliderInput("abProbesMaxDTm", "A-B Probes Max ∆Tm", 
                        0, 20, 1.2, step = 0.1),
            wellPanel(
                numericInput("oligoConc", "Oligo Conc, nM", 200),
                numericInput("naConc", "Na+ Conc, mM", 50),
                numericInput("mgConc", "Mg++ Conc, mM", 3),
                numericInput("dntpsConc", "dNTPs Conc, mM", 0.8)
            )
        ),

        mainPanel(
            tabsetPanel(
                tabPanel("Text Output",
                         verbatimTextOutput("txtResults")
                ),
                tabPanel("Table Output",
                         verbatimTextOutput("alnTxt"),
                         DT::dataTableOutput("tblResults"))
            )
        )
    )
)

server <- function(input, output) {
    MAX_DELTA_TM_FACTOR <- 0.8
    DEL_INS_PENALTY_TM <- 10
    complement <- function(ntseq) chartr("ATGCatgc", "TACGtacg", ntseq)
    reverse <- function(ntseq) stringi::stri_reverse(ntseq)
    
    probes <- reactive({
        dsnTask <- toupper(input$dsnTaskTxt)
        dsnTaskParsed <- as.vector(str_match(dsnTask, "([ATGC]+)\\(([ATGC]*)>([ATGC]*)\\)([ATGC]+)"))
        snpPosition <- str_locate(dsnTask, "\\([ATGC]*>[ATGC]*\\)")[1, "start"]
        
        flankSeq5 <- tolower(dsnTaskParsed[2])
        wtAllele <- toupper(dsnTaskParsed[3])
        wtAlleleL <- str_length(wtAllele)
        mAllele <- toupper(dsnTaskParsed[4])
        mAlleleL <- str_length(mAllele)
        flankSeq3 <- tolower(dsnTaskParsed[5])
        
        delInsType <- (wtAlleleL - mAlleleL) != 0
        
        wtSeq <- paste0(flankSeq5, wtAllele, flankSeq3)
        wtSeqL <- str_length(wtSeq)
        wtSeqRange <- 1:wtSeqL
        wtSeqRC <- reverse(complement(wtSeq))
        mSeq <- paste0(flankSeq5, mAllele, flankSeq3)
        mSeqRC <- reverse(complement(mSeq))
        
        # PROBE_MIN_MAX_LENGTH_NUC <- c(15, 30)
        # PROBE_MIN_MAX_TM <- c(48, 55)
        # PROBES_MAX_OVERLAP_TM <- 36
        # A_B_PROBES_MAX_DELTA_TM <- 1.2
        # PROBE_MIN_FLANK_NUC <- 3
        # OPTIMAL_TM <- 52
        # MAX_DELTA_TM_FACTOR <- 0.8
        
        
        probeLengthRange <- seq(input$probeMinMaxLength[1], input$probeMinMaxLength[2])
        
        calcProbes <- function(wtSeq, mSeq, probesOrientation) {
            if (probesOrientation == "B")
                snpPosition <- wtSeqL - snpPosition - wtAlleleL + 1
            probes <- 
                data.table(
                    probeStart = 
                        rep(seq(snpPosition - input$probeMinMaxLength[2] +
                                    input$probeMinFlank + wtAlleleL - 1,
                                snpPosition - input$probeMinFlank - 1),
                            each = length(probeLengthRange)
                        )
                )[
                    , probeStop := probeStart + probeLengthRange - 1, by = probeStart
                    ][
                        probeStop > snpPosition + input$probeMinFlank + (wtAlleleL - 1)
                        ][
                            , probeL := probeStop - probeStart + 1
                            ]
            
            probes <- 
                probes[
                    , probe := str_sub(wtSeq, probeStart, probeStop)
                    ][
                        , wtTm := Tm_NN(probe, 
                                        dnac1 = input$oligoConc,
                                        dnac2 = input$oligoConc,
                                        Na = input$naConc,
                                        Mg = input$mgConc,
                                        dNTPs = input$dntpsConc),
                        by = 1:NROW(probes)
                        ][
                            wtTm >= input$probeMinMaxTm[1] & wtTm <= input$probeMinMaxTm[2]
                            ][
                                , probe_m_comp := 
                                    complement(str_sub(mSeq, probeStart, probeStop))
                                ]
            if (delInsType)
                probes[, mTm := wtTm - DEL_INS_PENALTY_TM]
            else
                probes[, mTm := Tm_NN(probe, comSeq = probe_m_comp, 
                                      dnac1 = input$oligoConc,
                                      dnac2 = input$oligoConc,
                                      Na = input$naConc,
                                      Mg = input$mgConc,
                                      dNTPs = input$dntpsConc),
                       by = 1:length(probe)
                       ]
            probes[ , deltaTm := wtTm - mTm]
            colnames(probes) <- paste(colnames(probes), probesOrientation, sep = "_")
            probes
        }
        
        probesA <- calcProbes(wtSeq, mSeq, "A")
        probesB <- calcProbes(wtSeqRC, mSeqRC, "B")
        
        # crossjoin
        probes <- setkey(probesA[, c(k = 1, .SD)], k)[
            probesB[, c(k=1, .SD)], allow.cartesian = TRUE][, k:=NULL][
                , ABdeltaTm := abs(wtTm_A - wtTm_B)
                ][
                    ABdeltaTm < input$abProbesMaxDTm
                    ]
        probes <- 
            probes[
                , c("probeStart_B", "probeStop_B") :=
                    list(wtSeqL - probeStop_B + 1,
                         wtSeqL - probeStart_B + 1)][
                             , ABoverlap := {
                                 overlapRange <- 
                                     wtSeqRange[wtSeqRange %in% seq(probeStart_A, probeStop_A) &
                                                    wtSeqRange %in% seq(probeStart_B, probeStop_B)]
                                 str_sub(wtSeq,
                                         overlapRange[1],
                                         tail(overlapRange, 1))
                             },
                             by = 1:NROW(probes)]
        uniqueOverlaps <- unique(probes[, .(ABoverlap)])
        # for faster melting sort by start and stop.
        # Don't melt same start after ABoverlapTm > PROBES_MAX_OVERLAP_TM
        uniqueOverlaps <- uniqueOverlaps[
            , ABoverlapTm := Tm_NN(ABoverlap, 
                                   dnac1 = input$oligoConc,
                                   dnac2 = input$oligoConc,
                                   Na = input$naConc,
                                   Mg = input$mgConc,
                                   dNTPs = input$dntpsConc),
            by = 1:NROW(uniqueOverlaps)][
                ABoverlapTm < input$probeMaxOverlapTm
                ][
                    , ABoverlapL := str_length(ABoverlap)
                ]
        
        probes <- merge(probes, uniqueOverlaps, 
                        by.x = "ABoverlap", by.y = "ABoverlap", 
                        all.x = FALSE, all.y = FALSE)[
                            , c("deltaTmCoef", "probeLCoef") := 
                                list(deltaTm_A * deltaTm_B - (deltaTm_A - deltaTm_B) ^ 2,
                                     probeL_A * probeL_B - (probeL_A - probeL_B) ^ 2)
                            ][
                                deltaTmCoef > MAX_DELTA_TM_FACTOR * max(deltaTmCoef)
                                ]
        probes <- 
            probes[
                , probesTxt := 
                    sprintf(
                        paste("Probe A: %s, %i nt long, Tm = %.1f °C, ∆Tm(WT-M) = %.1f °C",
                              "Probe B: %s, %i nt long, Tm = %.1f °C, ∆Tm(WT-M) = %.1f °C",
                              "A-B Overlap: %s, %i nt long, Tm = %.1f °C\n",
                              "%s\n\n", sep = "\n"),
                        probe_A, probeL_A, wtTm_A, deltaTm_A,
                        probe_B, probeL_B, wtTm_B, deltaTm_B,
                        ABoverlap, ABoverlapL, ABoverlapTm,
                        paste(
                            str_pad(probe_A, width = probeStop_A, side = "left"),
                            wtSeq,
                            complement(wtSeq),
                            str_pad(reverse(probe_B), width = probeStop_B),
                            sep = "\n")
                    )
                ][
                    order(probeLCoef)
                    ]
        probes
    })

    output$txtResults <- renderText({
        req(probes())
        paste(probes()$probesTxt, collapse = "\n\n")
    })
    
    output$tblResults <- DT::renderDataTable({
        req(probes())
        dat <- probes()[
            , .("Probe A" = probe_A, 
                "Length" = probeL_A, 
                "Tm" = round(wtTm_A, 1), 
                "∆Tm(WT-M)" = round(deltaTm_A, 1),
                "Probe B" = probe_A, 
                "Length" = probeL_A, 
                "Tm" = round(wtTm_A, 1), 
                "∆Tm(WT-M)" = round(deltaTm_A, 1),
                "A-B Overlap" = ABoverlap,
                "Length" = ABoverlapL,
                "Tm" = round(ABoverlapTm, 1))]
        DT::datatable(
            data = dat,
            extensions = 'Buttons',
            options = list(paging = FALSE,
                           dom = 'Bfrtip',
                           buttons = c('copy', 'csv', 'excel')),
            selection = "single"
        )
    })
    
    output$alnTxt <- renderText({
        req(input$tblResults_rows_selected)
        probes()$probesTxt[input$tblResults_rows_selected]
    })
}

shinyApp(ui = ui, server = server)
