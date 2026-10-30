suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(mboost)
  library(kernelshap)
  library(shapviz)
  library(ggplot2)
  library(patchwork)
})

app_dir <- if (exists("appDir", inherits = FALSE)) appDir else getwd()
rdata_path <- file.path(app_dir, "shiny.Rdata")
if (!file.exists(rdata_path)) {
  stop("Missing shiny.Rdata in the app directory.")
}
load(rdata_path)

if (!exists("THRESH_MODE") || !THRESH_MODE %in% c("youden", "fixed")) {
  THRESH_MODE <- "youden"
}
if (!exists("THRESH_FIXED")) THRESH_FIXED <- 0.5
if (!exists("THRESH_YOUDEN")) THRESH_YOUDEN <- THRESH_FIXED

pred_fun <- function(object, newdata) {
  as.numeric(predict(object, newdata = newdata, type = "response"))
}

resolve_thresh <- function(mode) {
  if (identical(mode, "fixed")) THRESH_FIXED else THRESH_YOUDEN
}

SHAP_RB <- c("#FF0051", "#008BFB")

theme_shap_clean <- function(base_size = 13) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line.x = element_line(colour = "black", linewidth = 0.45),
      axis.ticks.y = element_blank(),
      axis.text = element_text(colour = "black", size = base_size - 1),
      axis.title = element_text(colour = "black", size = base_size),
      plot.title = element_blank(),
      plot.margin = margin(8, 12, 8, 8)
    )
}

build_input_row <- function(input) {
  vals <- lapply(feats, function(f) as.numeric(input[[paste0("x_", f)]]))
  names(vals) <- feats
  as.data.frame(vals, stringsAsFactors = FALSE)
}

card_wrap <- function(p, header) {
  p +
    labs(title = header) +
    theme(
      plot.title = element_text(
        hjust = 0, face = "bold", size = 12, colour = "#212529",
        margin = margin(b = 8)
      ),
      plot.background = element_rect(fill = "white", colour = "#DEE2E6", linewidth = 0.6),
      plot.margin = margin(10, 12, 10, 12)
    )
}

build_web_result_plot <- function(new_x, p_hat, sv, thr, thresh_mode) {
  risk_lab <- if (isTRUE(p_hat >= thr)) "Higher risk" else "Lower risk"
  thr_lab <- if (identical(thresh_mode, "youden")) {
    sprintf("Youden (%.4f)", thr)
  } else {
    sprintf("Fixed (%.2f)", thr)
  }

  page_title <- ggplot() +
    annotate(
      "text", x = 0, y = 0.5, label = "Sepsis-3 Risk Predictor",
      hjust = 0, vjust = 0.5, size = 7, fontface = "bold", colour = "#1F4E79"
    ) +
    xlim(0, 1) + ylim(0, 1) + theme_void() +
    theme(plot.margin = margin(2, 4, 2, 4))

  input_lines <- c("Threshold", paste0("  *  ", thr_lab), "")
  for (f in feats) {
    lab <- feat_meta$label[feat_meta$name == f][1]
    val <- as.numeric(new_x[[f]][1])
    if (identical(feat_meta$type[feat_meta$name == f][1], "binary")) {
      val_lab <- if (isTRUE(abs(val - 1) < 1e-8)) "Yes" else "No"
    } else {
      val_lab <- format(round(val, 3), nsmall = 3, trim = TRUE)
    }
    input_lines <- c(input_lines, lab, paste0("  ", val_lab), "")
  }
  input_lines <- c(input_lines, "[ Predict ]", "[ Download result PDF ]")
  p_input <- ggplot() +
    annotate(
      "text", x = 0.02, y = 0.98, label = paste(input_lines, collapse = "\n"),
      hjust = 0, vjust = 1, size = 3.15, colour = "#212529", lineheight = 1.15
    ) +
    xlim(0, 1) + ylim(0, 1) + theme_void()
  p_input <- card_wrap(p_input, "Input")

  p_pred <- ggplot() +
    annotate(
      "text", x = 0.5, y = 0.62, label = risk_lab,
      hjust = 0.5, vjust = 0.5, size = 8, fontface = "bold", colour = "#212529"
    ) +
    annotate(
      "text", x = 0.5, y = 0.28, label = sprintf("%.1f%%", 100 * p_hat),
      hjust = 0.5, vjust = 0.5, size = 6.5, fontface = "bold", colour = "#1F4E79"
    ) +
    xlim(0, 1) + ylim(0, 1) + theme_void()
  p_pred <- card_wrap(p_pred, "Prediction")

  p_wf <- sv_waterfall(
    sv, row_id = 1L, max_display = length(feats), fill_colors = SHAP_RB
  ) + theme_shap_clean(12)
  p_wf <- card_wrap(p_wf, "SHAP waterfall")

  p_fc <- sv_force(
    sv, row_id = 1L, max_display = length(feats), fill_colors = SHAP_RB
  ) + theme_shap_clean(12) +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.06)))
  p_fc <- card_wrap(p_fc, "SHAP force")

  right <- (p_pred / p_wf / p_fc) + plot_layout(heights = c(1.1, 3.6, 2.0))
  body <- (p_input | right) + plot_layout(widths = c(1.05, 2.2))
  (page_title / body) + plot_layout(heights = c(0.28, 5.8))
}

ui <- page_fluid(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#1F4E79",
    secondary = "#5B8FA8",
    "font-scale" = 0.95
  ),
  tags$head(
    tags$style(HTML("
      html, body { height: auto !important; overflow-y: auto !important; }
      .container-fluid { max-width: 1200px; padding-bottom: 48px; }
      .card { overflow: visible !important; }
      .card-body { overflow: visible !important; }
    "))
  ),
  h2("Sepsis-3 Risk Predictor", style = "color:#1F4E79; font-weight:700; margin:16px 0 20px 0;"),
  layout_columns(
    col_widths = c(4, 8),
    gap = "1.25rem",
    card(
      full_screen = FALSE,
      card_header("Input"),
      radioButtons(
        "thresh_mode",
        "Threshold",
        choiceNames = c(
          sprintf("Youden (%.4f)", THRESH_YOUDEN),
          sprintf("Fixed (%.2f)", THRESH_FIXED)
        ),
        choiceValues = c("youden", "fixed"),
        selected = THRESH_MODE
      ),
      lapply(seq_len(nrow(feat_meta)), function(i) {
        f <- feat_meta$name[[i]]
        lab <- feat_meta$label[[i]]
        if (identical(feat_meta$type[[i]], "binary")) {
          selectInput(
            inputId = paste0("x_", f),
            label = lab,
            choices = c("No" = 0, "Yes" = 1),
            selected = feat_meta$default[[i]]
          )
        } else {
          numericInput(
            inputId = paste0("x_", f),
            label = lab,
            value = round(feat_meta$default[[i]], 3),
            min = feat_meta$min[[i]],
            max = feat_meta$max[[i]],
            step = feat_meta$step[[i]]
          )
        }
      }),
      actionButton(
        "predict",
        "Predict",
        class = "btn-primary btn-lg w-100",
        style = "margin-top:4px; font-weight:600;"
      ),
      downloadButton(
        "download_pdf",
        "Download result PDF",
        class = "btn-outline-primary w-100",
        style = "margin-top:10px; font-weight:600;"
      )
    ),
    div(
      card(
        full_screen = FALSE,
        card_header("Prediction"),
        div(
          style = "text-align:center; padding:28px 12px;",
          h2(textOutput("risk_class"), style = "font-weight:700; margin:0;"),
          h3(textOutput("risk_prob"), style = "margin-top:12px; color:#1F4E79;")
        )
      ),
      br(),
      card(
        full_screen = FALSE,
        card_header("SHAP waterfall"),
        plotOutput("shap_waterfall", height = "480px")
      ),
      br(),
      card(
        full_screen = FALSE,
        card_header("SHAP force"),
        plotOutput("shap_force", height = "280px")
      )
    )
  )
)

server <- function(input, output, session) {
  active_thresh <- reactive({
    resolve_thresh(input$thresh_mode)
  })

  run_prediction <- function(new_x) {
    if (anyNA(as.matrix(new_x))) {
      stop("Please complete all seven features.")
    }
    p_hat <- as.numeric(pred_fun(fit_gb, new_x))
    ks <- kernelshap(
      object = fit_gb,
      X = new_x,
      bg_X = bg_X,
      pred_fun = pred_fun,
      exact = TRUE,
      verbose = FALSE
    )
    sv <- shapviz(ks, X = new_x)
    list(new_x = new_x, p_hat = p_hat, sv = sv)
  }

  pred_res <- reactiveVal(NULL)
  observeEvent(input$predict, {
    pred_res(run_prediction(build_input_row(input)))
  }, ignoreInit = TRUE)

  output$risk_class <- renderText({
    res <- pred_res()
    req(res)
    thr <- active_thresh()
    if (isTRUE(res$p_hat >= thr)) "Higher risk" else "Lower risk"
  })

  output$risk_prob <- renderText({
    res <- pred_res()
    req(res)
    sprintf("%.1f%%", 100 * res$p_hat)
  })

  output$shap_waterfall <- renderPlot({
    res <- pred_res()
    req(res)
    sv_waterfall(
      res$sv,
      row_id = 1L,
      max_display = length(feats),
      fill_colors = SHAP_RB
    ) +
      theme_shap_clean(13)
  }, res = 110)

  output$shap_force <- renderPlot({
    res <- pred_res()
    req(res)
    sv_force(
      res$sv,
      row_id = 1L,
      max_display = length(feats),
      fill_colors = SHAP_RB
    ) +
      theme_shap_clean(13) +
      theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      scale_x_continuous(expand = expansion(mult = c(0.06, 0.06)))
  }, res = 110)

  output$download_pdf <- downloadHandler(
    filename = function() {
      sprintf("sepsis3_prediction_%s.pdf", format(Sys.time(), "%Y%m%d_%H%M%S"))
    },
    content = function(file) {
      res <- pred_res()
      req(res)
      thr <- active_thresh()
      p_rep <- build_web_result_plot(
        new_x = res$new_x,
        p_hat = res$p_hat,
        sv = res$sv,
        thr = thr,
        thresh_mode = input$thresh_mode
      )
      ggsave(
        filename = file,
        plot = p_rep,
        width = 11,
        height = 9.2,
        units = "in",
        device = pdf,
        dpi = 300,
        bg = "white"
      )
    }
  )
}

shinyApp(ui = ui, server = server)
