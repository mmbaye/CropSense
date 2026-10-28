## ─────────────────────────────────────────────────────────────────────────────
##  CropSense
##  Multi-crop: Sorghum | Groundnut Grain | Groundnut Flour
##  Spectral axis in wavelength (nm)  •  Instrument auto-detection
## ─────────────────────────────────────────────────────────────────────────────

library(shiny)
library(bslib)
library(bsicons)
library(readxl)
library(readr)
library(plotly)
library(dplyr)
library(pls)
library(ggplot2)
library(RColorBrewer)
library(DT)
library(kernlab)

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0) x else y

## ── 1. Load all datasets ──────────────────────────────────────────────────────
read_nir <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  df[, !names(df) %in% "wl_NA"]          # drop mis-labelled last column
}

message("Loading datasets …")
DATA <- list(
  sorghum          = read_nir("NIRwetlab_sorghum.csv"),
  groundnut_grain  = read_nir("NIRwetlab_groundnut_grain.csv"),
  groundnut_flour  = read_nir("NIRwetlab_groundnut_flour.csv")
)

## ── 2. Dataset configuration ─────────────────────────────────────────────────
CROP_CONFIG <- list(

  sorghum = list(
    label         = "🌾 Sorghum",
    meta_cols     = c("lab_num","sample_id","genotype","croptype",
                      "type_of_material","protein","starch","lab_number"),
    traits        = c("protein","starch"),
    trait_colors  = c(protein = "#3B82F6", starch = "#F97316"),
    trait_units   = c(protein = "DW (%)",  starch = "DW (%)"),
    id_col        = "lab_num",
    color_choices = c("None","protein","starch","genotype","type_of_material")
  ),

  groundnut_grain = list(
    label         = "🥜 Groundnut – Grain",
    meta_cols     = c("lab_num","sample_id","sample_env","genotype","croptype",
                      "l","a","b","moisture","proteine","oil","lab_number"),
    traits        = c("proteine","oil","moisture","l","a","b"),
    trait_colors  = c(proteine="#3B82F6", oil="#F97316", moisture="#22C55E",
                      l="#8B5CF6", a="#EC4899", b="#EAB308"),
    trait_units   = c(proteine="DW (%)", oil="DW (%)", moisture="%", l="", a="", b=""),
    id_col        = "lab_num",
    color_choices = c("None","proteine","oil","moisture","genotype","croptype")
  ),

  groundnut_flour = list(
    label         = "🥜 Groundnut – Flour",
    meta_cols     = c("lab_num","sample_id","sample_env","genotype","croptype",
                      "l","a","b","moisture","proteine","oil","lab_number"),
    traits        = c("proteine","oil","moisture","l","a","b"),
    trait_colors  = c(proteine="#3B82F6", oil="#F97316", moisture="#22C55E",
                      l="#8B5CF6", a="#EC4899", b="#EAB308"),
    trait_units   = c(proteine="%", oil="%", moisture="%", l="", a="", b=""),
    id_col        = "lab_num",
    color_choices = c("None","proteine","oil","moisture","genotype","croptype")
  )
)

## ── 3. Helpers ────────────────────────────────────────────────────────────────
get_wave_info <- function(df) {
  cols <- grep("^wl_[0-9]", names(df), value = TRUE)
  vals <- as.numeric(sub("^wl_", "", cols))
  list(cols = cols, vals = vals)
}

## Preprocessing
## Impute all-NA rows to 0 so lm()-based methods never see "0 non-NA cases"
impute_na_rows <- function(X) {
  X[is.na(X)] <- 0
  X
}

snv <- function(X) {
  t(apply(X, 1, function(x) {
    mn <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - mn) / s
  }))
}

msc <- function(X) {
  X   <- impute_na_rows(X)
  ref <- colMeans(X, na.rm = TRUE)
  t(apply(X, 1, function(x) {
    tryCatch({
      fit <- lm(x ~ ref)
      b0  <- coef(fit)[1]; b1 <- coef(fit)[2]
      if (is.na(b1) || b1 == 0) return(x)
      (x - b0) / b1
    }, error = function(e) x)
  }))
}

detrend_spec <- function(X) {
  X <- impute_na_rows(X)
  t(apply(X, 1, function(x) {
    tryCatch({
      fit <- lm(x ~ seq_along(x))
      x - fitted(fit)
    }, error = function(e) x)
  }))
}

preprocess <- function(X, method, sg_w = 11, sg_p = 2) {
  sg <- requireNamespace("prospectr", quietly = TRUE)
  Xp <- switch(method,
    "Raw"            = X,
    "SNV"            = snv(X),
    "MSC"            = msc(X),
    "Detrend"        = detrend_spec(X),
    "SNV+Detrend"    = detrend_spec(snv(X)),
    "SG-Smooth"      = if (sg) prospectr::savitzkyGolay(X, m=0, p=sg_p, w=sg_w) else X,
    "1st Deriv (SG)" = if (sg) prospectr::savitzkyGolay(X, m=1, p=sg_p, w=sg_w) else diff(t(X)) |> t(),
    "2nd Deriv (SG)" = if (sg) prospectr::savitzkyGolay(X, m=2, p=sg_p, w=sg_w) else diff(t(X), differences=2) |> t(),
    X
  )
  if (!is.null(Xp) && is.matrix(Xp) && nrow(Xp) == nrow(X)) rownames(Xp) <- rownames(X)
  Xp
}

## ── 4. Pre-train calibration models at startup ───────────────────────────────
message("Training calibration PLS-R models …")
CALIB_MODELS <- list()

for (cid in names(CROP_CONFIG)) {
  cfg  <- CROP_CONFIG[[cid]]
  df   <- DATA[[cid]]
  wi   <- get_wave_info(df)
  Xraw <- as.matrix(df[, wi$cols])
  ## Remove all-NA rows before any model training
  good_rows <- apply(Xraw, 1, function(r) any(!is.na(r)))
  Xraw      <- Xraw[good_rows, , drop = FALSE]
  df_clean  <- df[good_rows, ]

  CALIB_MODELS[[cid]] <- lapply(setNames(cfg$traits, cfg$traits), function(trait) {
    y_full <- df_clean[[trait]]; keep <- !is.na(y_full) & is.finite(y_full)
    if (sum(keep) < 20) return(NULL)
    Xc <- preprocess(Xraw, "SNV"); Xc[is.na(Xc)] <- 0; Xc <- Xc[keep, ]
    y  <- y_full[keep]
    tryCatch({
      m    <- plsr(y ~ Xc, ncomp = min(10, nrow(Xc)-1),
                   method = "kernelpls", validation = "CV", segments = 10)
      lv   <- which.min(RMSEP(m)$val[1,1,-1])
      yhat <- as.vector(m$validation$pred[,1,lv])
      list(model=m, lv=lv, trait=trait, n=sum(keep),
           r2cv   = round(cor(y, yhat, use="complete.obs")^2, 3),
           rmsecv = round(sqrt(mean((y-yhat)^2, na.rm=TRUE)), 3),
           y_range = range(y, na.rm=TRUE))
    }, error = function(e) NULL)
  })
  CALIB_MODELS[[cid]] <- Filter(Negate(is.null), CALIB_MODELS[[cid]])
  message(sprintf("  %-20s → %d models ready", cid, length(CALIB_MODELS[[cid]])))
}

PREP_CHOICES <- c("Raw","SNV","MSC","Detrend","SNV+Detrend",
                  "1st Deriv (SG)","2nd Deriv (SG)","SG-Smooth")

## ── 5. UI ─────────────────────────────────────────────────────────────────────
ui <- page_navbar(
  title = tags$span(
    tags$img(src="https://placehold.co/28x28/3B82F6/white?text=N",
             height="28px", style="border-radius:6px;margin-right:8px;vertical-align:middle"),
    "CropSense"
  ),
  theme = bs_theme(
    bootswatch  = "flatly", primary = "#3B82F6",
    base_font   = font_google("Inter"),
    `navbar-bg` = "#1e293b"
  ),
  fillable = TRUE,

  ## Global crop selector in navbar
  nav_item(
    selectInput("crop_id", NULL,
                choices  = setNames(names(CROP_CONFIG),
                                    sapply(CROP_CONFIG, `[[`, "label")),
                selected = "sorghum",
                width    = "240px")
  ),

  ## ── Tab 1: Spectra ──────────────────────────────────────────────────────────
  nav_panel("📈 Spectra",
    layout_sidebar(
      sidebar = sidebar(width = 285, bg = "#f8fafc",
        h6("Preprocessing", class="text-primary fw-bold"),
        selectInput("prep_method", "Method", choices = PREP_CHOICES),
        conditionalPanel(
          "input.prep_method.includes('SG') || input.prep_method.includes('Deriv')",
          sliderInput("sg_w", "SG Window (odd)", min=5, max=21, value=11, step=2),
          sliderInput("sg_p", "SG Polynomial",   min=1, max=5,  value=2)
        ),
        hr(),
        h6("Display", class="text-primary fw-bold"),
        selectInput("color_by", "Colour by", choices = c("None")),
        uiOutput("wave_range_ui"),
        numericInput("n_display", "Max spectra shown", value=100, min=10, max=900, step=10),
        hr(),
        h6("Summary", class="text-muted"),
        verbatimTextOutput("spec_summary")
      ),
      plotlyOutput("spec_plot", height = "460px"),
      card(card_header("Dataset overview"),
        layout_columns(col_widths = c(3,3,3,3),
          value_box("Samples",     textOutput("ov_samples"),  showcase=bs_icon("layers"),        theme="primary"),
          value_box("Channels",    textOutput("ov_channels"), showcase=bs_icon("activity"),      theme="info"),
          value_box("Range (nm)",  textOutput("ov_range"),    showcase=bs_icon("arrows-expand"), theme="success"),
          value_box("Traits",      textOutput("ov_traits"),   showcase=bs_icon("clipboard-data"),theme="warning")
        )
      )
    )
  ),

  ## ── Tab 2: PCA ──────────────────────────────────────────────────────────────
  nav_panel("🔵 PCA",
    layout_sidebar(
      sidebar = sidebar(width = 285, bg = "#f8fafc",
        h6("Options", class="text-primary fw-bold"),
        selectInput("pca_prep",  "Preprocessing", choices = PREP_CHOICES[-1]),
        sliderInput("pca_ncomp", "Nr of PCs", min=2, max=20, value=10),
        selectInput("pca_x",    "PC X-axis", choices=paste0("PC",1:20), selected="PC1"),
        selectInput("pca_y",    "PC Y-axis", choices=paste0("PC",1:20), selected="PC2"),
        selectInput("pca_color", "Colour by", choices=c("None")),
        hr(),
        h6("Variance explained", class="text-muted"),
        plotOutput("scree_plot", height="160px")
      ),
      layout_columns(col_widths = c(7,5),
        card(card_header("Scores plot"),  plotlyOutput("pca_scores",   height="400px")),
        card(card_header("Loadings plot"),plotlyOutput("pca_loadings", height="400px"))
      )
    )
  ),

  ## ── Tab 3: PLS-R ────────────────────────────────────────────────────────────
  nav_panel("📊 PLS-R",
    layout_sidebar(
      sidebar = sidebar(width = 285, bg = "#f8fafc",
        h6("Model", class="text-primary fw-bold"),
        selectInput("pls_trait", "Y trait",        choices = c("protein")),
        selectInput("pls_prep",  "Preprocessing",  choices = PREP_CHOICES[-1]),
        sliderInput("pls_ncomp", "Max LV", min=2, max=20, value=10),
        selectInput("pls_val",   "Validation",
                    choices = c("CV (10-fold)"="CV","None"="none")),
        actionButton("pls_run", "▶ Run PLS-R", class="btn-primary w-100 mt-2"),
        hr(),
        verbatimTextOutput("pls_opt_lv"),
        hr(),
        downloadButton("pls_download", "Download predictions")
      ),
      layout_columns(col_widths = c(6,6),
        card(card_header("RMSECV by LV"),         plotlyOutput("pls_cv_plot",   height="320px")),
        card(card_header("Measured vs Predicted"), plotlyOutput("pls_pred_plot", height="320px"))
      ),
      card(card_header("Model metrics"),
        layout_columns(col_widths = c(3,3,3,3),
          value_box("R² (CV)",  textOutput("pls_r2cv"),  showcase=bs_icon("graph-up"),  theme="primary"),
          value_box("RMSECV",   textOutput("pls_rmsecv"),showcase=bs_icon("bullseye"),  theme="info"),
          value_box("R² (Cal)", textOutput("pls_r2cal"), showcase=bs_icon("bar-chart"), theme="success"),
          value_box("RMSEC",    textOutput("pls_rmsec"), showcase=bs_icon("rulers"),    theme="warning")
        )
      )
    )
  ),

  ## ── Tab 4: Predict ──────────────────────────────────────────────────────────
  nav_panel("🔮 Predict",
    layout_sidebar(
      sidebar = sidebar(width = 305, bg = "#f8fafc",

        h6("① Upload new NIR spectra", class="text-primary fw-bold"),
        fileInput("pred_file", NULL,
                  accept      = c(".xlsx",".csv"),
                  buttonLabel = "Browse file…",
                  placeholder = "No file selected"),
        div(style = "font-size:0.78rem;color:#64748b;",
          "Formats: .xlsx or .csv", tags$br(),
          "Spectral columns: ", tags$b("wl_400"), "…", tags$b("wl_2499.5"), tags$br(),
          "FOSS ≈ 400–2500 nm  |  Bruker ≈ 800–2500 nm", tags$br(),
          "Optional ID column: ", tags$b("NIRS_Code"), " or ", tags$b("lab_num")
        ),
        uiOutput("pred_instrument_badge"),

        hr(),
        h6("② Preprocessing", class="text-primary fw-bold"),
        selectInput("pred_prep", "Method",
                    choices  = PREP_CHOICES[-1],
                    selected = "SNV"),
        helpText("Use SNV (same as calibration) for best consistency.",
                 style="font-size:0.77rem;color:#64748b;"),

        hr(),
        h6("③ Run", class="text-primary fw-bold"),
        actionButton("pred_run", "▶  Predict traits",
                     class="btn-primary w-100", icon=icon("flask")),

        hr(),
        h6("④ Export", class="text-muted fw-bold"),
        downloadButton("pred_download", "⬇  Download predictions (.csv)",
                       class="btn-outline-success w-100"),

        hr(),
        h6("Calibration model quality", class="text-muted fw-bold"),
        DTOutput("pred_cal_table")
      ),

      uiOutput("pred_banner"),
      uiOutput("pred_cards"),
      layout_columns(col_widths = c(7,5),
        card(card_header("Uploaded spectra (preprocessed)"),
             plotlyOutput("pred_spec_plot", height="350px")),
        card(card_header("Predicted value distributions"),
             plotlyOutput("pred_dist_plot", height="350px"))
      ),
      card(card_header("Prediction results"), DTOutput("pred_table"))
    )
  ),

  ## ── Tab 5: Data ─────────────────────────────────────────────────────────────
  nav_panel("📋 Data",
    card(card_header("Dataset – metadata & trait values"),
         DTOutput("data_table"))
  ),

  ## ── Tab 6 : GP Regression + Uncertainty ────────────────────────────────────
  nav_panel("🎯 GP + Uncertainty",
    layout_sidebar(
      sidebar = sidebar(width = 300, bg = "#f8fafc",

        h6("① Configuration", class = "text-primary fw-bold"),
        selectInput("gp_trait", "Trait Y",         choices = c("protein")),
        selectInput("gp_prep",  "Preprocessing",
                    choices = c("SNV","MSC","Detrend","SNV+Detrend","Raw"),
                    selected = "SNV"),
        sliderInput("gp_nlv",   "PLS Components (dim. reduction)",
                    min = 3, max = 15, value = 8, step = 1),
        helpText("GP on PLS scores: fast and efficient.",
                 style = "font-size:0.77rem;color:#64748b;"),

        hr(),
        h6("② RBF Kernel", class = "text-primary fw-bold"),
        sliderInput("gp_ls",    "Correlation length (ℓ)",
                    min = 0.1, max = 5.0, value = 1.0, step = 0.1),
        sliderInput("gp_noise", "Noise σ² (×0.01)",
                    min = 1, max = 50, value = 5, step = 1),
        helpText("Large ℓ = smooth predictions. High noise = less confident model.",
                 style = "font-size:0.77rem;color:#64748b;"),
        hr(),
        checkboxInput("gp_auto_opt", "🔧 Auto-optimize hyperparameters",
                      value = TRUE),
        helpText("Finds optimal ℓ and noise automatically via marginal likelihood.",
                 style = "font-size:0.77rem;color:#22C55E;font-weight:600;"),

        hr(),
        h6("③ Alert Threshold", class = "text-primary fw-bold"),
        sliderInput("gp_ci_thresh", "95% CI Threshold (trait unit)",
                    min = 0.1, max = 5.0, value = 1.5, step = 0.1),
        helpText("Samples whose CI exceeds this threshold → to confirm in wet lab.",
                 style = "font-size:0.77rem;color:#64748b;"),

        hr(),
        actionButton("gp_run", "▶  Run GP Regression",
                     class = "btn-primary w-100", icon = icon("circle-nodes")),
        hr(),
        verbatimTextOutput("gp_opt_info"),
        hr(),
        downloadButton("gp_download", "⬇  Export results (.csv)",
                       class = "btn-outline-success w-100")
      ),

      ## Status banner
      uiOutput("gp_banner"),

      ## Metrics
      layout_columns(col_widths = c(3,3,3,3),
        value_box("R² GP (LOO)",    textOutput("gp_r2"),
                  showcase = bs_icon("graph-up"),    theme = "primary"),
        value_box("RMSECV GP",      textOutput("gp_rmse"),
                  showcase = bs_icon("bullseye"),    theme = "info"),
        value_box("Mean 95% CI",   textOutput("gp_ci_mean"),
                  showcase = bs_icon("rulers"),      theme = "success"),
        value_box("⚠ Alerts",      textOutput("gp_n_alert"),
                  showcase = bs_icon("exclamation-triangle"), theme = "warning")
      ),

      ## Scatter full width
      card(card_header("Measured vs Predicted  +  95% CI"),
           plotlyOutput("gp_scatter", height = "420px")),

      ## All samples table – alert rows highlighted in red
      card(card_header("📋  All samples — rows in red = to confirm in wet lab"),
           DTOutput("gp_alert_table"))
    )
  )
)

## ── 6. SERVER ─────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## ── Reactive accessors ──────────────────────────────────────────────────────
  cur_cfg    <- reactive(CROP_CONFIG[[req(input$crop_id)]])
  cur_df     <- reactive(DATA[[req(input$crop_id)]])
  cur_wi     <- reactive(get_wave_info(cur_df()))

  ## Shared good-row mask: drops the 8 all-NA ghost rows in groundnut data.
  ## ALL downstream reactives (Xraw, meta, df_clean) use the same filter
  ## so rows are always aligned.
  cur_good   <- reactive({
    wi <- cur_wi(); df <- cur_df()
    X  <- as.matrix(df[, wi$cols])
    apply(X, 1, function(r) any(!is.na(r)))
  })
  cur_df_clean <- reactive({ df <- cur_df(); df[cur_good(), ] })

  cur_Xraw   <- reactive({
    wi   <- cur_wi()
    df_c <- cur_df_clean()
    X    <- as.matrix(df_c[, wi$cols])
    X[is.na(X)] <- 0          # impute residual scattered NAs
    rownames(X) <- as.character(df_c[[cur_cfg()$id_col]])
    X
  })
  cur_meta   <- reactive({
    df_c <- cur_df_clean(); cfg <- cur_cfg()
    df_c[, cfg$meta_cols]
  })
  cur_models <- reactive(CALIB_MODELS[[req(input$crop_id)]])

  ## ── Update UI when crop changes ─────────────────────────────────────────────
  observeEvent(input$crop_id, {
    cfg <- cur_cfg(); wi <- cur_wi()
    updateSelectInput(session, "color_by",   choices = cfg$color_choices)
    updateSelectInput(session, "pca_color",  choices = cfg$color_choices)
    updateSelectInput(session, "pls_trait",  choices = cfg$traits)
    updateSelectInput(session, "gp_trait",   choices = cfg$traits)
    pred_res(NULL); new_spectra(NULL); gp_res(NULL)
  })

  ## ── Wave range slider (rendered dynamically) ────────────────────────────────
  output$wave_range_ui <- renderUI({
    wi <- cur_wi()
    mn <- floor(min(wi$vals)); mx <- ceiling(max(wi$vals))
    sliderInput("wave_range", "Wavelength range (nm)",
                min=mn, max=mx, value=c(mn, mx), step=1)
  })

  ## ── Overview value boxes ────────────────────────────────────────────────────
  output$ov_samples  <- renderText(nrow(cur_df_clean()))
  output$ov_channels <- renderText(length(cur_wi()$cols))
  output$ov_range    <- renderText({
    wi <- cur_wi()
    paste0(floor(min(wi$vals)), "–", ceiling(max(wi$vals)), " nm")
  })
  output$ov_traits   <- renderText(length(cur_cfg()$traits))

  ## ── Preprocessed matrix reactive ────────────────────────────────────────────
  X_prep <- reactive({
    preprocess(cur_Xraw(), input$prep_method,
               sg_w = input$sg_w %||% 11,
               sg_p = input$sg_p %||% 2)
  })

  ## ── SPECTRA PLOT ────────────────────────────────────────────────────────────
  output$spec_plot <- renderPlotly({
    X  <- X_prep()
    wi <- cur_wi()
    wn <- wi$vals[seq_len(ncol(X))]
    wr <- req(input$wave_range)
    idx <- which(wn >= wr[1] & wn <= wr[2])
    X  <- X[, idx]; wn <- wn[idx]
    n  <- min(input$n_display %||% 100, nrow(X))
    sel <- sample(nrow(X), n)
    X   <- X[sel, ]
    meta <- cur_meta()

    if (input$color_by == "None" || !input$color_by %in% names(meta)) {
      cols <- rep("#3B82F6", n)
    } else {
      v <- meta[[input$color_by]][sel]
      if (is.numeric(v)) {
        pal <- colorRampPalette(c("#3B82F6","#22C55E","#EAB308","#EF4444"))(100)
        rng <- range(v, na.rm=TRUE); dif <- diff(rng)
        vn  <- if (is.finite(dif) && dif > 0) round((v-rng[1])/dif*99)+1 else rep(50L,length(v))
        vn[is.na(vn)|vn<1|vn>100] <- 50L; cols <- pal[vn]
      } else {
        lvls <- unique(na.omit(v))
        pal  <- colorRampPalette(brewer.pal(min(9,max(3,length(lvls))),"Set1"))(length(lvls))
        idx2 <- match(v, lvls); idx2[is.na(idx2)] <- 1L; cols <- pal[idx2]
      }
    }

    p <- plot_ly(type="scatter", mode="lines") %>%
      layout(xaxis = list(title="Wavelength (nm)", gridcolor="#e2e8f0"),
             yaxis = list(title=paste0(input$prep_method," Absorbance"), gridcolor="#e2e8f0"),
             paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", showlegend=FALSE)
    for (i in seq_len(nrow(X)))
      p <- p %>% add_trace(x=wn, y=X[i,], line=list(color=cols[i], width=0.6),
                           text=rownames(X)[i], hoverinfo="text+y")
    p
  })

  output$spec_summary <- renderText({
    X <- X_prep()
    paste0("n = ", nrow(X), " spectra\nChannels: ", ncol(X),
           "\nMean abs: ", round(mean(X,na.rm=TRUE),4),
           "\nSD abs:   ", round(sd(X,na.rm=TRUE),4))
  })

  ## ── PCA ─────────────────────────────────────────────────────────────────────
  pca_res <- reactive({
    Xp <- preprocess(cur_Xraw(), input$pca_prep); Xp[is.na(Xp)] <- 0
    prcomp(Xp, center=TRUE, scale.=FALSE, rank.=input$pca_ncomp)
  })

  output$scree_plot <- renderPlot({
    pca <- pca_res()
    ve  <- (pca$sdev^2/sum(pca$sdev^2))[1:min(10,length(pca$sdev))]*100
    df  <- data.frame(PC=factor(paste0("PC",seq_along(ve)),
                                levels=paste0("PC",seq_along(ve))), Var=ve)
    ggplot(df,aes(PC,Var)) + geom_col(fill="#3B82F6",alpha=0.8) +
      geom_line(aes(group=1),color="#1e293b",linewidth=0.6) +
      geom_point(color="#1e293b",size=2) +
      labs(x=NULL, y="Variance (%)") + theme_minimal(base_size=11) +
      theme(axis.text.x=element_text(angle=45,hjust=1))
  })

  output$pca_scores <- renderPlotly({
    pca     <- pca_res()
    scores  <- as.data.frame(pca$x)
    var_exp <- round(pca$sdev^2/sum(pca$sdev^2)*100, 1)
    pc_x    <- if (input$pca_x %in% names(scores)) input$pca_x else "PC1"
    pc_y    <- if (input$pca_y %in% names(scores)) input$pca_y else "PC2"
    pci     <- as.integer(gsub("PC","",pc_x)); pcj <- as.integer(gsub("PC","",pc_y))
    ids     <- as.character(cur_df_clean()[[cur_cfg()$id_col]])
    bl <- list(
      xaxis = list(title=paste0(pc_x," (",var_exp[pci],"%)"),
                   gridcolor="#e2e8f0",zeroline=TRUE,zerolinecolor="#94a3b8"),
      yaxis = list(title=paste0(pc_y," (",var_exp[pcj],"%)"),
                   gridcolor="#e2e8f0",zeroline=TRUE,zerolinecolor="#94a3b8"),
      paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)"
    )
    meta <- cur_meta(); cv <- input$pca_color
    xv <- scores[[pc_x]]; yv <- scores[[pc_y]]
    if (cv=="None" || !cv %in% names(meta))
      return(plot_ly(x=xv,y=yv,text=ids,type="scatter",mode="markers",showlegend=FALSE,
                     marker=list(color="#3B82F6",size=7,opacity=0.75)) %>% layout(bl))
    cv_raw <- meta[[cv]]
    if (is.numeric(cv_raw)) {
      cv_c <- cv_raw; cv_c[is.na(cv_c)] <- mean(cv_c,na.rm=TRUE)
      return(plot_ly(x=xv,y=yv,text=ids,type="scatter",mode="markers",showlegend=FALSE,
                     marker=list(color=cv_c,colorscale="RdYlGn",size=7,opacity=0.80,
                                 colorbar=list(title=cv,thickness=14,len=0.7),
                                 showscale=TRUE)) %>% layout(bl))
    }
    cv_ch <- as.character(cv_raw); cv_ch[is.na(cv_ch)] <- "NA"
    lvls  <- sort(unique(cv_ch))
    pal   <- colorRampPalette(brewer.pal(min(9,max(3,length(lvls))),"Set1"))(length(lvls))
    p <- plot_ly(x=xv,y=yv,text=ids,type="scatter",mode="markers",showlegend=FALSE,
                 marker=list(color=pal[match(cv_ch,lvls)],size=7,opacity=0.80))
    for (k in seq_along(lvls))
      p <- p %>% add_trace(x=NA_real_,y=NA_real_,type="scatter",mode="markers",
                           name=lvls[k],showlegend=TRUE,hoverinfo="none",
                           marker=list(color=pal[k],size=9))
    p %>% layout(bl)
  })

  output$pca_loadings <- renderPlotly({
    pca  <- pca_res(); wi <- cur_wi()
    wvls <- wi$vals[seq_len(nrow(pca$rotation))]
    pc_x <- if (input$pca_x %in% colnames(pca$rotation)) input$pca_x else "PC1"
    pc_y <- if (input$pca_y %in% colnames(pca$rotation)) input$pca_y else "PC2"
    ld   <- data.frame(wn=wvls, LX=pca$rotation[,pc_x], LY=pca$rotation[,pc_y])
    plot_ly(ld,x=~wn,y=~LX,type="scatter",mode="lines",name=pc_x,
            line=list(color="#3B82F6")) %>%
      add_trace(y=~LY,name=pc_y,line=list(color="#EF4444")) %>%
      layout(xaxis=list(title="Wavelength (nm)",gridcolor="#e2e8f0"),
             yaxis=list(title="Loading",gridcolor="#e2e8f0"),
             paper_bgcolor="rgba(0,0,0,0)",plot_bgcolor="rgba(0,0,0,0)",
             legend=list(x=0.01,y=0.99))
  })

  ## ── PLS-R ───────────────────────────────────────────────────────────────────
  pls_model <- eventReactive(input$pls_run, {
    withProgress(message="Fitting PLS-R …", value=0.3, {
      Xp <- preprocess(cur_Xraw(), input$pls_prep, sg_w=11, sg_p=2)
      Xp[is.na(Xp)] <- 0
      df <- cur_df_clean(); y_full <- df[[input$pls_trait]]; keep <- !is.na(y_full)
      Xp <- Xp[keep, ]; y <- y_full[keep]
      val_method <- if (input$pls_val=="CV") "CV" else "none"
      incProgress(0.4)
      m <- plsr(y ~ Xp, ncomp = min(input$pls_ncomp, nrow(Xp)-1),
                method = "kernelpls", validation = val_method,
                segments = if (val_method=="CV") 10 else NULL)
      incProgress(0.3)
      list(model=m, y=y, keep=keep, Xp=Xp)
    })
  })

  pls_ncomp_opt <- reactive({
    req(pls_model()); m <- pls_model()$model
    if (!is.null(m$validation)) which.min(RMSEP(m)$val[1,1,-1]) else m$ncomp
  })

  output$pls_opt_lv <- renderText({ req(pls_model()); paste("Optimal LV:", pls_ncomp_opt()) })

  output$pls_cv_plot <- renderPlotly({
    req(pls_model()); m <- pls_model()$model
    if (is.null(m$validation)) return(plot_ly() %>% layout(title="No CV performed"))
    rmsep <- RMSEP(m)$val[1,1,]
    df    <- data.frame(LV=0:(length(rmsep)-1), RMSECV=rmsep)
    opt   <- pls_ncomp_opt()
    plot_ly(df,x=~LV,y=~RMSECV,type="scatter",mode="lines+markers",
            line=list(color="#3B82F6"),marker=list(color="#3B82F6")) %>%
      add_segments(x=opt,xend=opt,y=min(rmsep)*.97,yend=max(rmsep)*1.03,
                   line=list(color="#EF4444",dash="dash"),showlegend=FALSE) %>%
      layout(xaxis=list(title="Latent Variables",gridcolor="#e2e8f0"),
             yaxis=list(title="RMSECV",gridcolor="#e2e8f0"),
             paper_bgcolor="rgba(0,0,0,0)",plot_bgcolor="rgba(0,0,0,0)")
  })

  output$pls_pred_plot <- renderPlotly({
    req(pls_model()); obj <- pls_model(); m <- obj$model; nc <- pls_ncomp_opt(); y <- obj$y
    yhat <- if (!is.null(m$validation)) as.vector(m$validation$pred[,1,nc])
            else as.vector(fitted(m)[,1,nc])
    df   <- data.frame(Measured=y, Predicted=yhat)
    lim  <- range(c(df$Measured,df$Predicted), na.rm=TRUE)
    plot_ly(df,x=~Measured,y=~Predicted,type="scatter",mode="markers",
            marker=list(color="#3B82F6",size=6,opacity=0.7)) %>%
      add_segments(x=lim[1],xend=lim[2],y=lim[1],yend=lim[2],
                   line=list(color="#EF4444",dash="dash"),showlegend=FALSE) %>%
      layout(xaxis=list(title=paste("Measured",input$pls_trait),gridcolor="#e2e8f0"),
             yaxis=list(title="Predicted",gridcolor="#e2e8f0"),
             paper_bgcolor="rgba(0,0,0,0)",plot_bgcolor="rgba(0,0,0,0)")
  })

  pls_metrics <- reactive({
    req(pls_model()); obj <- pls_model(); m <- obj$model; y <- obj$y; nc <- pls_ncomp_opt()
    yhat_cal <- as.vector(fitted(m)[,1,nc])
    rmsec    <- sqrt(mean((y-yhat_cal)^2, na.rm=TRUE))
    r2cal    <- cor(y,yhat_cal,use="complete.obs")^2
    if (!is.null(m$validation)) {
      yhat_cv <- as.vector(m$validation$pred[,1,nc])
      rmsecv  <- sqrt(mean((y-yhat_cv)^2, na.rm=TRUE))
      r2cv    <- cor(y,yhat_cv,use="complete.obs")^2
    } else { rmsecv <- NA; r2cv <- NA }
    list(rmsec=rmsec, r2cal=r2cal, rmsecv=rmsecv, r2cv=r2cv)
  })

  output$pls_r2cv   <- renderText({ req(pls_metrics()); m <- pls_metrics(); if (is.na(m$r2cv))   "N/A" else round(m$r2cv,3) })
  output$pls_rmsecv <- renderText({ req(pls_metrics()); m <- pls_metrics(); if (is.na(m$rmsecv)) "N/A" else round(m$rmsecv,3) })
  output$pls_r2cal  <- renderText({ req(pls_metrics()); round(pls_metrics()$r2cal,3) })
  output$pls_rmsec  <- renderText({ req(pls_metrics()); round(pls_metrics()$rmsec,3) })

  output$pls_download <- downloadHandler(
    filename = function() paste0("PLS_", input$pls_trait, "_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(pls_model()); obj <- pls_model(); nc <- pls_ncomp_opt()
      yhat <- as.vector(fitted(obj$model)[,1,nc])
      df   <- cur_df_clean(); ids <- df[[cur_cfg()$id_col]][obj$keep]
      write.csv(data.frame(ID=ids, Measured=obj$y, Predicted=yhat), file, row.names=FALSE)
    })

  ## ── PREDICT TAB ─────────────────────────────────────────────────────────────
  output$pred_cal_table <- renderDT({
    mdls <- cur_models()
    df <- do.call(rbind, lapply(mdls, function(m)
      data.frame(Trait=m$trait, LV=m$lv, R2_CV=m$r2cv, RMSECV=m$rmsecv, N=m$n)))
    datatable(df, rownames=FALSE, options=list(dom="t", pageLength=10),
              class="table table-sm") %>%
      formatStyle("R2_CV",
                  backgroundColor=styleInterval(c(0.7,0.9),c("#fef2f2","#fefce8","#f0fdf4")))
  })

  new_spectra <- reactiveVal(NULL)
  pred_res    <- reactiveVal(NULL)

  ## Instrument badge
  output$pred_instrument_badge <- renderUI({
    ns <- new_spectra()
    if (is.null(ns)) return(NULL)
    col <- if (grepl("FOSS", ns$instrument)) "#22C55E" else "#3B82F6"
    div(class="mt-2 p-2 rounded",
        style=paste0("background:",col,"22;border-left:3px solid ",col,";font-size:0.78rem;"),
        tags$b(ns$instrument), tags$br(),
        paste0(ns$n_matched, "/", ns$n_calib, " channels matched"))
  })

  ## File upload
  observeEvent(input$pred_file, {
    req(input$pred_file)
    ext <- tools::file_ext(input$pred_file$name)
    withProgress(message = "Reading file …", value = 0.5, {
      df_new <- tryCatch(
        if (tolower(ext)=="xlsx") read_excel(input$pred_file$datapath)
        else read_csv(input$pred_file$datapath, show_col_types=FALSE),
        error = function(e) {
          showNotification(paste("Cannot read file:", e$message), type="error"); NULL }
      )
      req(!is.null(df_new))

      wi_cal  <- cur_wi()
      new_wc  <- grep("^wl_[0-9]", names(df_new), value=TRUE)

      if (length(new_wc) < 50) {
        showNotification("Too few spectral columns found. Check column names.", type="error")
        return()
      }

      # Instrument detection
      min_nm     <- min(as.numeric(sub("^wl_","",new_wc)), na.rm=TRUE)
      instrument <- if (min_nm <= 500) "FOSS (400–2500 nm)" else "Bruker (~800–2500 nm)"

      # Align to calibration columns
      matched <- intersect(new_wc, wi_cal$cols)
      if (length(matched) < 50) {
        showNotification(paste0("Only ", length(matched),
          " spectral columns matched calibration wavenumbers."),
          type="error", duration=10)
        return()
      }

      X_new <- matrix(0, nrow=nrow(df_new), ncol=length(wi_cal$cols),
                      dimnames=list(NULL, wi_cal$cols))
      X_new[, matched] <- as.matrix(df_new[, matched])

      ids <- if ("NIRS_Code" %in% names(df_new)) as.character(df_new$NIRS_Code)
             else if ("lab_num"  %in% names(df_new)) as.character(df_new$lab_num)
             else if ("Sample"   %in% names(df_new)) as.character(df_new$Sample)
             else paste0("S", seq_len(nrow(df_new)))
      rownames(X_new) <- ids

      new_spectra(list(X=X_new, ids=ids,
                       n_matched=length(matched), n_calib=length(wi_cal$cols),
                       instrument=instrument, fname=input$pred_file$name))
      pred_res(NULL)
      showNotification(
        paste0("✅ ", nrow(X_new), " spectra loaded  |  ", instrument,
               "  |  ", length(matched), " channels matched"),
        type="message", duration=5)
    })
  })

  ## Banner
  output$pred_banner <- renderUI({
    ns  <- new_spectra()
    lbl <- cur_cfg()$label
    if (is.null(ns))
      div(class="alert alert-info m-2",
          tags$b("Crop: "), lbl, tags$br(),
          tags$b("① "), "Upload a NIR file (.xlsx or .csv).", tags$br(),
          tags$b("② "), "Select preprocessing.", tags$br(),
          tags$b("③ "), "Click ", tags$b("▶ Predict"), ".")
    else
      div(class="alert alert-success m-2",
          tags$b(paste0("✅ ", ns$fname)), "  |  ",
          paste0(nrow(ns$X), " samples  |  "), ns$instrument, "  |  ",
          paste0(ns$n_matched, "/", ns$n_calib, " channels"),
          "  →  Click ", tags$b("▶ Predict"))
  })

  ## Run predictions
  observeEvent(input$pred_run, {
    ns <- new_spectra(); req(ns)
    mdls  <- cur_models(); df_cal <- cur_df_clean(); Xraw_cal <- cur_Xraw()

    withProgress(message = "Running PLS-R predictions …", value = 0.1, {
      incProgress(0.2, detail="Preprocessing …")
      X_prep_new <- tryCatch({
        Xp <- preprocess(ns$X, input$pred_prep, sg_w=11, sg_p=2)
        Xp[is.na(Xp)] <- 0; Xp
      }, error=function(e) {
        showNotification(paste("Preprocessing error:", e$message), type="error"); NULL })
      req(!is.null(X_prep_new))

      results <- data.frame(NIRS_Code=ns$ids, stringsAsFactors=FALSE)
      step    <- 0.6 / max(1, length(mdls))

      for (trait in names(mdls)) {
        mi <- mdls[[trait]]
        incProgress(step, detail=paste0("Predicting ", trait, " …"))

        y_full <- df_cal[[trait]]; keep <- !is.na(y_full)
        Xc <- tryCatch({
          Xcc <- preprocess(Xraw_cal, input$pred_prep, sg_w=11, sg_p=2)
          Xcc[is.na(Xcc)] <- 0; Xcc[keep, ]
        }, error=function(e) NULL)

        preds <- tryCatch({
          req(!is.null(Xc))
          lv_u  <- min(mi$lv, nrow(Xc)-1, ncol(Xc))
          m_new <- plsr(y_full[keep] ~ Xc, ncomp=lv_u, method="kernelpls")
          as.vector(predict(m_new, newdata=X_prep_new, ncomp=lv_u))
        }, error=function(e) {
          tryCatch(as.vector(predict(mi$model, newdata=X_prep_new, ncomp=mi$lv)),
                   error=function(e2) rep(NA_real_, nrow(X_prep_new)))
        })

        buf   <- diff(mi$y_range) * 0.15
        preds <- pmax(pmin(preds, mi$y_range[2]+buf), mi$y_range[1]-buf)
        results[[paste0(trait,"_pred")]] <- round(preds, 3)
      }

      pred_res(results)
      showNotification(paste0("✅ Predicted ", nrow(results), " samples for ",
                               length(mdls), " traits."), type="message", duration=4)
    })
  })

  ## Metric cards
  output$pred_cards <- renderUI({
    df <- pred_res(); req(df)
    cfg <- cur_cfg(); mdls <- cur_models()
    cards <- lapply(names(mdls), function(trait) {
      pc   <- paste0(trait,"_pred")
      col  <- cfg$trait_colors[[trait]] %||% "#3B82F6"
      un   <- cfg$trait_units[[trait]]  %||% ""
      vals <- df[[pc]]
      mu   <- round(mean(vals,na.rm=TRUE),2); sg <- round(sd(vals,na.rm=TRUE),2)
      rng  <- paste0(round(min(vals,na.rm=TRUE),1)," – ",round(max(vals,na.rm=TRUE),1)," ",un)
      mi   <- mdls[[trait]]
      value_box(
        title    = gsub("_"," ",trait),
        value    = paste0(mu," ± ",sg," ",un),
        showcase = bs_icon("graph-up-arrow"),
        p(paste0("n=",sum(!is.na(vals)),"  range: ",rng), style="font-size:0.73rem;"),
        p(paste0("Cal: R²=",mi$r2cv,"  RMSECV=",mi$rmsecv," ",un),
          style=paste0("font-size:0.70rem;color:",col,";")),
        theme="primary", style=paste0("border-left:4px solid ",col,";")
      )
    })
    ncw <- 12 %/% max(1, min(3, length(cards)))
    layout_columns(col_widths=rep(ncw, length(cards)), .list=cards)
  })

  ## Spectra preview
  output$pred_spec_plot <- renderPlotly({
    ns <- new_spectra(); req(ns)
    X_p <- tryCatch({
      Xp <- preprocess(ns$X, input$pred_prep, sg_w=11, sg_p=2)
      Xp[is.na(Xp)] <- 0; Xp }, error=function(e) ns$X)
    wn  <- cur_wi()$vals[seq_len(ncol(X_p))]
    pal <- colorRampPalette(c("#3B82F6","#22C55E","#F97316","#EF4444"))(max(1,nrow(X_p)))
    p   <- plot_ly(type="scatter", mode="lines") %>%
      layout(xaxis=list(title="Wavelength (nm)",gridcolor="#e2e8f0"),
             yaxis=list(title=paste0(input$pred_prep," Absorbance"),gridcolor="#e2e8f0"),
             paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", showlegend=FALSE)
    for (i in seq_len(min(nrow(X_p),200)))
      p <- p %>% add_trace(x=wn, y=X_p[i,], line=list(color=pal[i], width=0.7),
                           text=rownames(X_p)[i], hoverinfo="text+y")
    p
  })

  ## Distribution box plot
  output$pred_dist_plot <- renderPlotly({
    df <- pred_res(); req(df)
    cfg <- cur_cfg(); mdls <- cur_models()
    p   <- plot_ly(type="box")
    for (trait in names(mdls)) {
      pc  <- paste0(trait,"_pred")
      col <- cfg$trait_colors[[trait]] %||% "#3B82F6"
      un  <- cfg$trait_units[[trait]]  %||% ""
      p <- p %>% add_trace(y=df[[pc]],
                           name=paste0(gsub("_"," ",trait)," (",un,")"),
                           fillcolor=paste0(substr(col,1,7),"44"),
                           line=list(color=col), marker=list(color=col,size=3,opacity=0.5))
    }
    p %>% layout(yaxis=list(title="Predicted value",gridcolor="#e2e8f0"),
                 paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                 showlegend=FALSE)
  })

  ## Results table
  output$pred_table <- renderDT({
    df <- pred_res(); req(df)
    nm <- names(df); names(df) <- gsub("_pred"," (predicted)",nm)
    datatable(df, rownames=FALSE, filter="top",
              options=list(pageLength=20, scrollX=TRUE,
                           columnDefs=list(list(className="dt-center",targets="_all"))),
              class="table table-sm table-hover") %>%
      formatRound(grep("predicted",names(df),value=TRUE), digits=3)
  })

  output$pred_download <- downloadHandler(
    filename = function() paste0("NIR_pred_",input$crop_id,"_",Sys.Date(),".csv"),
    content  = function(file) { df <- pred_res(); req(df); write.csv(df, file, row.names=FALSE) })

  ## ── GP REGRESSION + INCERTITUDE ────────────────────────────────────────────

  ## GP core : noyau RBF + posterior ──────────────────────────────────────────
  gp_rbf_kernel <- function(A, B, ls = 1.0) {
    # Noyau RBF (Squared Exponential) : K(x,x') = exp(-||x-x'||² / 2ℓ²)
    d2 <- outer(rowSums(A^2), rowSums(B^2), "+") - 2 * (A %*% t(B))
    d2[d2 < 0] <- 0   # numerical stability
    exp(-d2 / (2 * ls^2))
  }

  gp_predict <- function(X_tr, y_tr, X_te,
                          ls = 1.0, sigma2 = 1.0, noise = 0.05) {
    # Training covariance matrix + regularization
    K_tt <- gp_rbf_kernel(X_tr, X_tr, ls) * sigma2
    diag(K_tt) <- diag(K_tt) + noise

    # Cholesky pour inversion stable
    L <- tryCatch(t(chol(K_tt)), error = function(e) {
      diag(K_tt) <- diag(K_tt) + 1e-4; t(chol(K_tt))
    })

    alpha <- backsolve(t(L), forwardsolve(L, y_tr))

    K_st <- gp_rbf_kernel(X_te, X_tr, ls) * sigma2
    mu   <- as.vector(K_st %*% alpha)

    # Variance posterieure
    v      <- forwardsolve(L, t(K_st))
    K_ss   <- sigma2   # diagonale uniquement (diagonal RBF kernel = sigma2)
    var_p  <- pmax(K_ss - colSums(v^2), 1e-8)
    list(mean = mu, variance = var_p, std = sqrt(var_p))
  }

  ## Optimize GP hyperparameters via log marginal likelihood ──────────────────
  gp_optimize_params <- function(X, y) {
    sigma2 <- var(y)
    neg_log_ml <- function(log_params) {
      ls    <- exp(log_params[1])
      noise <- exp(log_params[2])
      K <- gp_rbf_kernel(X, X, ls) * sigma2
      diag(K) <- diag(K) + noise
      tryCatch({
        L     <- t(chol(K))
        alpha <- backsolve(t(L), forwardsolve(L, y))
        ld    <- 2 * sum(log(pmax(diag(L), 1e-15)))
        as.numeric(0.5 * (t(y) %*% alpha + ld + length(y) * log(2*pi)))
      }, error = function(e) 1e10)
    }
    opt <- tryCatch(
      optim(c(0, -3), neg_log_ml, method = "L-BFGS-B",
            lower = c(-3, -8), upper = c(3, 0),
            control = list(maxit = 200, factr = 1e9)),
      error = function(e) list(par = c(0, -3))
    )
    list(ls = exp(opt$par[1]), noise = exp(opt$par[2]))
  }

  ## Leave-One-Out GP (cross-validation honnête) ─────────────────────────────
  gp_loo <- function(X, y, ls = 1.0, sigma2 = 1.0, noise = 0.05) {
    n      <- nrow(X)
    mu_loo <- numeric(n)
    sd_loo <- numeric(n)

    K  <- gp_rbf_kernel(X, X, ls) * sigma2
    diag(K) <- diag(K) + noise

    L  <- tryCatch(t(chol(K)), error = function(e) {
      diag(K) <<- diag(K) + 1e-3; t(chol(K))
    })
    K_inv <- chol2inv(t(L))

    # Formules LOO exactes (Rasmussen & Williams 2006, §5.4.2)
    # Avoids n regressions: exact O(n³) formula (Rasmussen & Williams 2006)
    alpha      <- as.vector(K_inv %*% y)
    K_inv_diag <- diag(K_inv)
    mu_loo     <- y - alpha / K_inv_diag
    var_loo    <- pmax(1 / K_inv_diag, 1e-8)
    sd_loo     <- sqrt(var_loo)

    list(mean = mu_loo, variance = var_loo, std = sd_loo)
  }

  ## Réactif résultat GP ────────────────────────────────────────────────────
  gp_res <- reactiveVal(NULL)

  ## Lancement GP ─────────────────────────────────────────────────────────────
  observeEvent(input$gp_run, {
    withProgress(message = "GP Regression…", value = 0.1, {

      incProgress(0.15, detail = "Preprocessing…")
      Xp <- tryCatch(
        { Xr <- preprocess(cur_Xraw(), input$gp_prep); Xr[is.na(Xr)] <- 0; Xr },
        error = function(e) { showNotification(paste("Erreur prep:", e$message), type="error"); NULL }
      )
      req(!is.null(Xp))

      df   <- cur_df_clean()
      y_full <- df[[input$gp_trait]]
      keep   <- !is.na(y_full) & is.finite(y_full)
      Xp <- Xp[keep, ]; y <- y_full[keep]
      ids <- as.character(df[[cur_cfg()$id_col]])[keep]

      incProgress(0.25, detail = "PLS reduction…")
      # Dimensionality reduction: PLS scores
      nlv <- min(input$gp_nlv, nrow(Xp) - 1, ncol(Xp))
      pls_red <- pls::plsr(y ~ Xp, ncomp = nlv, method = "kernelpls")
      scores  <- as.matrix(pls_red$scores)   # n × nlv

      # Normaliser les scores (important pour le noyau RBF)
      sc_mean <- colMeans(scores)
      sc_sd   <- apply(scores, 2, sd)
      sc_sd[sc_sd < 1e-10] <- 1
      scores_n <- scale(scores, center = sc_mean, scale = sc_sd)

      incProgress(0.25, detail = "Optimizing hyperparameters…")
      sigma2 <- var(y)

      # Auto-optimize hyperparameters via marginal likelihood
      # (overrides manual sliders when auto-optimize is on)
      if (input$gp_auto_opt) {
        withProgress(message = "Optimizing GP hyperparameters…", value = 0.5, {
          opt_params <- gp_optimize_params(scores_n, y)
          ls    <- opt_params$ls
          noise <- opt_params$noise
        })
      } else {
        ls    <- input$gp_ls
        noise <- input$gp_noise * 0.01
      }

      incProgress(0.15, detail = "LOO GP…")
      loo <- tryCatch(
        gp_loo(scores_n, y, ls = ls, sigma2 = sigma2, noise = noise),
        error = function(e) {
          showNotification(paste("Erreur GP:", e$message), type = "error"); NULL
        }
      )
      req(!is.null(loo))

      incProgress(0.25, detail = "Computing metrics…")
      ci95   <- 1.96 * loo$std
      r2_gp  <- cor(y, loo$mean, use = "complete.obs")^2
      rmse_gp <- sqrt(mean((y - loo$mean)^2, na.rm = TRUE))

      # Alert: CI 95% exceeds threshold
      alert <- ci95 > input$gp_ci_thresh

      result_df <- data.frame(
        ID            = ids,
        Measured      = round(y, 3),
        GP_predicted  = round(loo$mean, 3),
        GP_std        = round(loo$std,  3),
        CI95_lower    = round(loo$mean - ci95, 3),
        CI95_upper    = round(loo$mean + ci95, 3),
        CI95_width    = round(ci95 * 2, 3),
        Wetlab_alert  = alert,
        stringsAsFactors = FALSE
      )

      gp_res(list(
        df        = result_df,
        y         = y,
        ids       = ids,
        r2        = r2_gp,
        rmse      = rmse_gp,
        ci_mean   = mean(ci95 * 2),
        n_alert   = sum(alert),
        trait     = input$gp_trait,
        units     = cur_cfg()$trait_units[[input$gp_trait]] %||% "",
        ls_used   = ls,
        noise_used = noise
      ))

      showNotification(
        paste0("✅ GP done — R²=", round(r2_gp, 3),
               "  RMSECV=", round(rmse_gp, 3),
               "  ⚠ ", sum(alert), " alerts"),
        type = "message", duration = 5
      )
    })
  })

  ## Bannière ─────────────────────────────────────────────────────────────────
  output$gp_banner <- renderUI({
    if (is.null(gp_res()))
      div(class = "alert alert-info m-2",
          tags$b("How does GP work?"), tags$br(),
          "① Spectral preprocessing is applied (SNV recommended).", tags$br(),
          "② Spectra are projected onto", tags$b("n PLS components"), " (4200→n reduction).", tags$br(),
          "③ A", tags$b("Gaussian Process (RBF kernel)"), "is fitted on these scores.", tags$br(),
          "④ ", tags$b("Leave-One-Out validation"), " gives a prediction AND a ", tags$b("95% confidence interval"), " for each sample.", tags$br(),
          "⑤ Samples with CI > threshold are ", tags$b("flagged for wet lab confirmation."))
    else {
      gp <- gp_res()
      div(class = "alert alert-success m-2",
          tags$b(paste0("✅ GP — ", gp$trait)),
          paste0("  R²=", round(gp$r2,3),
                 "  RMSECV=", round(gp$rmse,3), " ", gp$units,
                 "  Mean CI95=", round(gp$ci_mean,3), " ", gp$units,
                 "  ⚠ ", gp$n_alert, " alerts"))
    }
  })

  ## Métriques ────────────────────────────────────────────────────────────────
  output$gp_r2      <- renderText({ req(gp_res()); round(gp_res()$r2,    3) })
  output$gp_rmse    <- renderText({ req(gp_res()); round(gp_res()$rmse,  3) })
  output$gp_ci_mean <- renderText({ req(gp_res()); round(gp_res()$ci_mean, 3) })
  output$gp_n_alert <- renderText({ req(gp_res()); gp_res()$n_alert })

  output$gp_opt_info <- renderText({
    req(gp_res())
    gp <- gp_res()
    paste0("Trait    : ", gp$trait, "\n",
           "R² LOO   : ", round(gp$r2,   3), "\n",
           "RMSECV   : ", round(gp$rmse, 3), " ", gp$units, "\n",
           "Mean CI95: ", round(gp$ci_mean, 3), " ", gp$units, "\n",
           "Alerts   : ", gp$n_alert, " / ", nrow(gp$df), "\n",
           "ℓ        : ", round(gp$ls_used, 3), "\n",
           "noise    : ", round(gp$noise_used, 4))
  })

  ## Scatter mesuré vs prédit + barres d'erreur IC95 ─────────────────────────
  output$gp_scatter <- renderPlotly({
    req(gp_res()); gp <- gp_res(); df <- gp$df
    un  <- gp$units
    lim <- range(c(df$CI95_lower, df$CI95_upper, df$Measured), na.rm = TRUE)

    # Couleur selon alerte
    cols <- ifelse(df$Alerte_wetlab, "#EF4444", "#3B82F6")
    symb <- ifelse(df$Alerte_wetlab, "diamond", "circle")

    plot_ly() %>%
      # Barres d'erreur (IC95)
      add_trace(
        x = df$Measured, y = df$GP_predicted,
        type = "scatter", mode = "markers",
        error_y = list(
          type       = "data",
          symmetric  = FALSE,
          array      = df$GP_predicted - df$CI95_lower,   # lower CI
          arrayminus = df$CI95_upper    - df$GP_predicted,  # upper CI
          color      = "#94A3B8",
          thickness  = 1.2, width = 4
        ),
        marker = list(color = cols, symbol = symb, size = 8, opacity = 0.85),
        text   = paste0(df$ID, "<br>Measured=", df$Measured,
                        "<br>Predicted=", df$GP_predicted,
                        "<br>CI95=[", df$CI95_lower, ", ", df$CI95_upper, "]",
                        "<br>", ifelse(df$Wetlab_alert, "⚠ ALERT", "")),
        hoverinfo = "text",
        showlegend = FALSE
      ) %>%
      # Ligne 1:1
      add_segments(x = lim[1], xend = lim[2], y = lim[1], yend = lim[2],
                   line = list(color = "#EF4444", dash = "dash", width = 1.5),
                   showlegend = FALSE) %>%
      # Manual legend (alert / ok)
      add_trace(x = NA_real_, y = NA_real_, type = "scatter", mode = "markers",
                name = "OK", marker = list(color = "#3B82F6", symbol = "circle", size = 10),
                showlegend = TRUE) %>%
      add_trace(x = NA_real_, y = NA_real_, type = "scatter", mode = "markers",
                name = "⚠ Wet lab alert",
                marker = list(color = "#EF4444", symbol = "diamond", size = 10),
                showlegend = TRUE) %>%
      layout(
        xaxis      = list(title = paste("Measured", un),  gridcolor = "#e2e8f0"),
        yaxis      = list(title = paste("GP Predicted", un), gridcolor = "#e2e8f0"),
        paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
        legend = list(x = 0.02, y = 0.98, bgcolor = "rgba(255,255,255,0.8)")
      )
  })

  ## Alert table – all samples, alert rows in red ──────────────────────────────
  output$gp_alert_table <- renderDT({
    req(gp_res()); gp <- gp_res()

    # All samples, sorted: alerts first then by decreasing CI width
    df      <- gp$df
    df_meta <- cur_df_clean()
    cfg     <- cur_cfg()

    # Add genotype column from metadata (always show it)
    gen_col <- intersect(c("genotype","Genotype","GENOTYPE"), names(df_meta))
    if (length(gen_col) > 0) {
      df$Genotype <- df_meta[[gen_col[1]]][match(df$ID, as.character(df_meta[[cfg$id_col]]))]
      # Reorder: ID, Genotype, then rest
      df <- df[, c("ID","Genotype", setdiff(names(df), c("ID","Genotype")))]
    }

    # Sort: alerts first, then by CI95 width descending
    df <- df[order(-df$Wetlab_alert, -df$CI95_width), , drop = FALSE]

    # Store alert flag for row coloring then remove the boolean col
    alert_flag      <- df$Wetlab_alert
    df$Wetlab_alert <- ifelse(df$Wetlab_alert, "⚠ Yes", "✓ OK")
    names(df)[names(df) == "Wetlab_alert"] <- "Wet lab?"

    datatable(
      df, rownames = FALSE, filter = "top",
      options = list(
        pageLength = 20, scrollX = TRUE,
        columnDefs = list(list(className = "dt-center", targets = "_all"))
      ),
      class = "table table-sm table-hover"
    ) %>%
      # Entire row red for alert samples
      formatStyle(
        "Wet lab?",
        target           = "row",
        backgroundColor  = styleEqual(c("⚠ Yes", "✓ OK"),
                                       c("#FEE2E2", "white")),
        color            = styleEqual(c("⚠ Yes", "✓ OK"),
                                       c("#991B1B", "#374151")),
        fontWeight       = styleEqual(c("⚠ Yes", "✓ OK"),
                                       c("bold", "normal"))
      ) %>%
      # CI95 width gradient
      formatStyle(
        "CI95_width",
        background = styleColorBar(range(df$CI95_width, na.rm=TRUE), "#BFDBFE"),
        backgroundSize = "100% 80%", backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      ) %>%
      formatRound(c("Measured","GP_predicted","GP_std","CI95_lower","CI95_upper","CI95_width"),
                  digits = 3)
  })

  ## Export ───────────────────────────────────────────────────────────────────
  output$gp_download <- downloadHandler(
    filename = function() paste0("GP_", input$gp_trait, "_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(gp_res()); write.csv(gp_res()$df, file, row.names = FALSE)
    })

  ## ── DATA TABLE ──────────────────────────────────────────────────────────────

  output$data_table <- renderDT({
    df  <- cur_df_clean(); cfg <- cur_cfg()
    df_meta <- df[, cfg$meta_cols]

    # Always bring genotype column first (after ID column)
    id_col  <- cfg$id_col
    gen_col <- intersect(c("genotype","Genotype","GENOTYPE"), names(df_meta))
    other   <- setdiff(names(df_meta), c(id_col, gen_col))
    col_order <- c(id_col, gen_col, other)
    col_order <- col_order[col_order %in% names(df_meta)]
    df_meta <- df_meta[, col_order]

    df_meta %>%
      mutate(across(where(is.numeric), ~round(.,3))) %>%
      datatable(
        options = list(
          pageLength  = 20, scrollX = TRUE,
          columnDefs  = list(list(className = "dt-center", targets = "_all")),
          # Pin genotype column visible
          fixedColumns = list(leftColumns = min(3, ncol(df_meta)))
        ),
        rownames = FALSE, filter = "top",
        class    = "table table-sm table-hover"
      ) %>%
      formatStyle(
        gen_col,
        fontWeight       = "bold",
        backgroundColor  = styleEqual(unique(df_meta[[gen_col[1]]]),
                             colorRampPalette(c("#DBEAFE","#FEF9C3","#DCFCE7","#FEE2E2","#EDE9FE"))(
                               length(unique(df_meta[[gen_col[1]]]))))
      )
  })
}

## ── 7. Run ────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
