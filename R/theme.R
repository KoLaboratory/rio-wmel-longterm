# -----------------------------------------------------------------------------
# Shared plotting helpers and project paths
# -----------------------------------------------------------------------------
#
# The figures in this project are built on ggplot2::theme_bw() with per-figure
# text sizes. These helpers are optional conveniences; they do not change the
# appearance of the published figures unless you opt into them.
# -----------------------------------------------------------------------------

# Projected CRS used throughout (SIRGAS 2000 / UTM zone 23S, metres).
PROJ_CRS <- 31983

# Standard output directory for figures and tables.
OUT_DIR <- "outputs"

#' theme_bw() with a configurable base text size.
#'
#' A thin wrapper so a whole figure's text size can be set in one place.
#'
#' @param base_size Base font size (points).
#' @return A ggplot2 theme object.
theme_wmel <- function(base_size = 13) {
  ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = base_size))
}

#' Save a figure to the outputs directory in one or more formats.
#'
#' @param plot A ggplot (or patchwork) object.
#' @param name File stem (no extension), e.g. "fig2".
#' @param width,height Dimensions in inches.
#' @param formats Character vector of extensions, e.g. c("pdf", "png").
#' @param dir Output directory (defaults to OUT_DIR).
save_fig <- function(plot, name, width, height,
                     formats = c("pdf", "png"), dir = OUT_DIR) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  for (fmt in formats) {
    ggplot2::ggsave(
      filename = file.path(dir, paste0(name, ".", fmt)),
      plot = plot, width = width, height = height
    )
  }
  invisible(file.path(dir, paste0(name, ".", formats)))
}
