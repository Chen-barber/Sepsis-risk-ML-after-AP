options(shiny.launch.browser = TRUE)
app_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
port <- as.integer(Sys.getenv("SHINY_PORT", "8080"))
host <- Sys.getenv("SHINY_HOST", "127.0.0.1")
url <- sprintf("http://%s:%s", host, port)
message("Local calculator URL: ", url)
message("Press Esc or Ctrl+C to stop the app.")
shiny::runApp(
  appDir = app_dir,
  host = host,
  port = port,
  launch.browser = TRUE
)
